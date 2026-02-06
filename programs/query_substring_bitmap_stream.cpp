// query_substring_bitmap_stream_windowed_mix_complete.cpp
//
// Windowed shard mixing for prefix-diverse output, while keeping lexicographic scanning
// inside each shard (fast iterator). COMPLETE for construct_k > k0 exhaustion.
//
// Features:
// - Filters: GC% (min/max), substring, construct_k > k0 (expansion via left/right appends)
// - Pagination with cursor
// - Random shard permutation (deterministic with seed) when --random_access is set
// - NEW: --window W loads W shards concurrently and interleaves outputs across them,
//        producing much more prefix diversity than single-shard draining.
// - NEW: --reverse_complement (default off). When set and --substring is set,
//        matches substring OR reverse-complement(substring).
//
// Output:
//   __META__ <cursor> <hasMore 0/1> <returned_count> <kout>
//   then one k-mer per line.
//
// Cursor (BCW2):
//  magic 'B','C','W','2'
//  flags(u8): bit0=random_access
//  k0(u8), kout(u8), d(u8)
//  numShards(u32)
//  seed(u64)                 // permutation seed if random_access, else 0
//  next_perm_pos(u32)        // next perm position to load when a lane finishes
//  window(u16)               // W, must match request
//  burst(u16)                // burst, must match request (kept for sanity)
//  lane_count(u16)           // == window
//  then for each lane i in [0..W-1]:
//     active(u8)
//     if active:
//        perm_pos(u32)       // which shard (by perm position) this lane is scanning
//        mode(u8): 0=konly, 1=expand
//        if mode==0:
//           after(u64)       // last returned 18-mer id from this shard; UINT64_MAX means "not started"
//        else:
//           parent_anchor(u64)    // UINT64_MAX means "not started"
//           child_present(u8)     // 1 if mid-parent expansion
//           if child_present:
//              L(u8), left_idx(u64), right_idx(u64)
//
// Notes:
// - For kout==k0, we can still use GC-hist skipping at shard selection time.
// - For construct_k>k0, GC-hist skipping is disabled (post-filter at kout can resurrect candidates).
//
// Compile:
//   g++ -O3 -march=native -mtune=native -std=c++17 -pthread \
//     query_substring_bitmap_stream_windowed_mix_complete.cpp -lroaring \
//     -o query_substring_bitmap_stream_windowed_mix
//
// Example:
//   ./query_substring_bitmap_stream_windowed_mix \
//     --shards shards_18 --gc-hist gc_hist_shards_18.json \
//     --window 64 --burst 1 --limit 200 --threads 4 \
//     --random_access --ra_seed 12345 \
//     --substring CGCGCC --reverse_complement --gc-min 40 --gc-max 60 --construct_k 20
//

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include <sys/resource.h>

#include <roaring/roaring64.h>

using namespace std;
using Clock = std::chrono::steady_clock;
using Sec   = std::chrono::duration<double>;

static inline long peak_rss_kb() { rusage r; getrusage(RUSAGE_SELF, &r); return r.ru_maxrss; }

// ---------------- KBITv1 portable ----------------
struct KbitHeader {
  uint64_t total_bits=0, ones=0, k=0, seed=0, flags=0, payload_len=0;
};
static inline uint64_t read_le64(const unsigned char* p) {
  uint64_t v=0; for(int i=0;i<8;i++) v |= (uint64_t)p[i] << (8*i); return v;
}
static roaring64_bitmap_t* load_kbit_portable(const string& path, KbitHeader& h) {
  ifstream in(path, ios::binary);
  if (!in) { perror(("open " + path).c_str()); return nullptr; }

  unsigned char hdr[64];
  in.read((char*)hdr, 64);
  if (!in || memcmp(hdr, "KBITv1\0", 8) != 0) {
    cerr << "Invalid shard header: " << path << "\n";
    return nullptr;
  }
  h.total_bits  = read_le64(hdr + 8);
  h.ones        = read_le64(hdr + 16);
  h.k           = read_le64(hdr + 24);
  h.seed        = read_le64(hdr + 32);
  h.flags       = read_le64(hdr + 40);
  h.payload_len = read_le64(hdr + 48);

  if (h.flags != 2) {
    cerr << "Shard not portable flags=2: " << path << "\n";
    return nullptr;
  }

  vector<char> payload(h.payload_len);
  in.read(payload.data(), (streamsize)h.payload_len);
  if ((uint64_t)in.gcount() != h.payload_len) {
    cerr << "Truncated shard payload: " << path << "\n";
    return nullptr;
  }
  roaring64_bitmap_t* bm = roaring64_bitmap_portable_deserialize_safe(payload.data(), payload.size());
  if (!bm) cerr << "Deserialize failed for shard: " << path << "\n";
  return bm;
}

// ---------------- DNA helpers ----------------
static inline int base4_digit(char c){
  switch(c){
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return -1;
  }
}
static inline char base4_char(int d){ return "ACGT"[d&3]; }
static string decode_kmer(uint64_t val, int k){
  string s(k,' ');
  for(int i=k-1;i>=0;--i){ s[i]=base4_char((int)(val&3ULL)); val>>=2; }
  return s;
}

static inline char dna_comp(char c) {
  switch (c) {
    case 'A': return 'T';
    case 'a': return 't';
    case 'C': return 'G';
    case 'c': return 'g';
    case 'G': return 'C';
    case 'g': return 'c';
    case 'T': return 'A';
    case 't': return 'a';
    default:  return '?';
  }
}

static string revcomp_string(const string& s) {
  string out;
  out.resize(s.size());
  for (size_t i = 0; i < s.size(); ++i) {
    out[s.size() - 1 - i] = dna_comp(s[i]);
  }
  return out;
}

// ---------------- Filters ----------------
static inline bool passes_gc_percent(uint64_t v, int k, int gcMinPct, int gcMaxPct) {
  int gc = 0;
  for (int i = 0; i < k; ++i) {
    uint64_t d = v & 3ULL;
    gc += (d == 1ULL) | (d == 2ULL);
    v >>= 2;
  }
  const int lhs = gc * 100;
  const int lo = gcMinPct * k;
  const int hi = gcMaxPct * k;
  return lhs >= lo && lhs <= hi;
}

// Substring patterns for fast check via masks
struct Pattern { uint64_t mask, bits; };
static inline bool contains_sub(uint64_t v, const Pattern* pats, size_t np) {
  for(size_t i=0;i<np;i++){
    const auto& p=pats[i];
    if(((v ^ p.bits) & p.mask) == 0ULL) return true;
  }
  return false;
}

// ---------------- index.json parsing ----------------
static bool read_index(const string& dir, unsigned& numShards, vector<string>& files,
                       uint64_t& k_out, uint64_t& total_bits_out,
                       vector<uint64_t>& starts, vector<uint64_t>& ends) {
  ifstream in(dir + "/index.json");
  if (!in) return false;

  string line;
  numShards = 0;
  files.clear();
  k_out = 0;
  total_bits_out = 0;
  starts.clear();
  ends.clear();

  auto parse_u64_field = [&](const string& s, const string& key, uint64_t& out)->bool {
    auto pos = s.find(key);
    if (pos == string::npos) return false;
    pos = s.find(':', pos);
    if (pos == string::npos) return false;
    pos++;
    while (pos < s.size() && (s[pos] == ' ' || s[pos] == '\t')) pos++;
    size_t end = s.find_first_of(",}", pos);
    if (end == string::npos || end <= pos) return false;
    out = stoull(s.substr(pos, end - pos));
    return true;
  };

  while (getline(in, line)) {
    if (line.find("\"num_shards\"") != string::npos) {
      auto p = line.find(':');
      if (p != string::npos) numShards = (unsigned)stoul(line.substr(p+1));
    }
    if (line.find("\"total_bits\"") != string::npos) {
      auto p = line.find(':');
      if (p != string::npos) total_bits_out = (uint64_t)stoull(line.substr(p+1));
    }
    if (line.find("\"k\"") != string::npos && line.find("\"seed\"") == string::npos) {
      auto p = line.find(':');
      if (p != string::npos) k_out = (uint64_t)stoull(line.substr(p+1));
    }
    auto fpos = line.find("\"file\"");
    if (fpos != string::npos) {
      uint64_t sstart=0, send=0;
      bool has_start = parse_u64_field(line, "\"start\"", sstart);
      bool has_end = parse_u64_field(line, "\"end\"", send);

      auto colon = line.find(':', fpos);
      if (colon == string::npos) continue;
      auto s1 = line.find('"', colon);
      if (s1 == string::npos) continue;
      auto s2 = line.find('"', s1 + 1);
      if (s2 == string::npos) continue;
      files.push_back(line.substr(s1 + 1, s2 - (s1 + 1)));

      if (has_start && has_end) {
        starts.push_back(sstart);
        ends.push_back(send);
      }
    }
  }

  if (numShards == 0) numShards = (unsigned)files.size();
  if (files.size() != numShards) {
    files.clear();
    files.reserve(numShards);
    for (unsigned i = 0; i < numShards; ++i) {
      ostringstream os;
      os << "shard_" << setw(4) << setfill('0') << i << ".kbit";
      files.push_back(os.str());
    }
  }
  return numShards > 0 && !files.empty() && k_out > 0;
}

// ---------------- GC histogram JSON parsing (minimal) ----------------
static bool load_gc_hist_json(const string& path, int& k_out, vector<vector<uint64_t>>& hists_out) {
  ifstream in(path);
  if (!in) return false;
  string s((istreambuf_iterator<char>(in)), istreambuf_iterator<char>());

  auto skip_ws = [&](size_t& i){
    while (i < s.size()) {
      char c = s[i];
      if (c==' '||c=='\n'||c=='\r'||c=='\t') i++;
      else break;
    }
  };
  auto parse_int = [&](size_t& i, long long& out)->bool{
    skip_ws(i);
    bool neg=false;
    if (i < s.size() && s[i]=='-') { neg=true; i++; }
    if (i >= s.size() || s[i]<'0' || s[i]>'9') return false;
    long long v=0;
    while (i < s.size() && s[i]>='0' && s[i]<='9') { v = v*10 + (s[i]-'0'); i++; }
    out = neg ? -v : v;
    return true;
  };

  // "k"
  {
    auto pos = s.find("\"k\"");
    if (pos == string::npos) return false;
    size_t i = pos;
    if ((i = s.find(':', i)) == string::npos) return false;
    i++;
    long long kk=0;
    if (!parse_int(i, kk)) return false;
    k_out = (int)kk;
    if (k_out <= 0 || k_out > 32) return false;
  }

  // optional "num_shards"
  int num_shards = -1;
  {
    auto pos = s.find("\"num_shards\"");
    if (pos != string::npos) {
      size_t i = pos;
      if ((i = s.find(':', i)) != string::npos) {
        i++;
        long long ns=0;
        if (parse_int(i, ns)) num_shards = (int)ns;
      }
    }
  }

  if (num_shards > 0) {
    hists_out.assign((size_t)num_shards, vector<uint64_t>((size_t)k_out+1, 0));
  } else {
    hists_out.clear();
  }

  size_t i = 0;
  while (true) {
    auto sp = s.find("\"shard\"", i);
    if (sp == string::npos) break;
    size_t j = sp;
    if ((j = s.find(':', j)) == string::npos) break;
    j++;
    long long shard_id_ll=0;
    if (!parse_int(j, shard_id_ll)) { i = sp+7; continue; }
    int shard_id = (int)shard_id_ll;
    i = j;

    auto gh = s.find("\"gc_hist\"", i);
    if (gh == string::npos) break;
    size_t kpos = gh;
    if ((kpos = s.find('[', kpos)) == string::npos) break;
    kpos++;

    if (shard_id < 0) { i = gh+8; continue; }
    if ((size_t)shard_id >= hists_out.size()) {
      hists_out.resize((size_t)shard_id + 1, vector<uint64_t>((size_t)k_out+1, 0));
    }

    for (int b=0; b<=k_out; ++b) {
      long long v=0;
      if (!parse_int(kpos, v)) return false;
      hists_out[(size_t)shard_id][(size_t)b] = (uint64_t)v;
      skip_ws(kpos);
      if (b < k_out) {
        if (kpos < s.size() && s[kpos] == ',') kpos++;
      }
    }
    auto close = s.find(']', kpos);
    if (close == string::npos) break;
    i = close + 1;
  }
  return true;
}

// ---------------- splitmix + perm ----------------
static inline uint64_t splitmix64(uint64_t x) {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}
static vector<uint32_t> build_perm(uint32_t n, uint64_t seed) {
  vector<uint32_t> p(n);
  for (uint32_t i=0;i<n;i++) p[i]=i;
  uint64_t st = seed ? seed : 1;
  for (uint32_t i=n; i>1; --i) {
    st = splitmix64(st);
    uint32_t j = (uint32_t)(st % i);
    swap(p[i-1], p[j]);
  }
  return p;
}

// ---------------- Expansion helpers ----------------
static inline uint64_t pow4(int n) { return (n<=0) ? 1ULL : (1ULL << (2*n)); }
static inline uint64_t make_value(uint64_t parentB, int k0, int kout,
                                  int L, uint64_t left_idx, uint64_t right_idx) {
  int d = kout - k0;
  int R = d - L;
  return (left_idx << (2 * (k0 + R))) | (parentB << (2 * R)) | right_idx;
}
static inline void init_state_first(int d, uint8_t& L, uint64_t& left_idx, uint64_t& right_idx) {
  L = (uint8_t)d;
  left_idx = 0;
  right_idx = 0;
}
static inline bool advance_state(int d, uint8_t& L, uint64_t& left_idx, uint64_t& right_idx) {
  int Li = (int)L;
  int R = d - Li;

  uint64_t right_lim = pow4(R);
  right_idx++;
  if (right_idx < right_lim) return true;

  right_idx = 0;
  uint64_t left_lim = pow4(Li);
  left_idx++;
  if (left_idx < left_lim) return true;

  left_idx = 0;
  if (Li == 0) return false;
  L = (uint8_t)(Li - 1);
  return true;
}

// ---------------- Base64url + LE pack/unpack ----------------
static const char* B64URL = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_";
static string b64url_encode(const vector<uint8_t>& in) {
  string out;
  out.reserve((in.size() * 4 + 2) / 3);
  uint32_t buf = 0;
  int bits = 0;
  for (uint8_t c : in) {
    buf = (buf << 8) | c;
    bits += 8;
    while (bits >= 6) {
      bits -= 6;
      out.push_back(B64URL[(buf >> bits) & 63]);
    }
  }
  if (bits) out.push_back(B64URL[(buf << (6 - bits)) & 63]);
  return out;
}
static bool b64url_decode(const string& s, vector<uint8_t>& out) {
  int8_t T[256];
  memset(T, -1, sizeof(T));
  for (int i = 0; i < 64; ++i) T[(unsigned char)B64URL[i]] = (int8_t)i;

  uint32_t buf = 0;
  int bits = 0;
  out.clear();
  out.reserve((s.size() * 3) / 4);

  for (char ch : s) {
    int8_t v = T[(unsigned char)ch];
    if (v < 0) return false;
    buf = (buf << 6) | (uint32_t)v;
    bits += 6;
    if (bits >= 8) {
      bits -= 8;
      out.push_back((uint8_t)((buf >> bits) & 0xFF));
    }
  }
  return true;
}
static inline void push_u16_le(vector<uint8_t>& b, uint16_t x) {
  b.push_back((uint8_t)(x & 0xFF));
  b.push_back((uint8_t)((x >> 8) & 0xFF));
}
static inline bool read_u16_le(const vector<uint8_t>& b, size_t off, uint16_t& x) {
  if (off + 2 > b.size()) return false;
  x = (uint16_t)b[off] | ((uint16_t)b[off+1] << 8);
  return true;
}
static inline void push_u32_le(vector<uint8_t>& b, uint32_t x) {
  b.push_back((uint8_t)(x & 0xFF));
  b.push_back((uint8_t)((x >> 8) & 0xFF));
  b.push_back((uint8_t)((x >> 16) & 0xFF));
  b.push_back((uint8_t)((x >> 24) & 0xFF));
}
static inline bool read_u32_le(const vector<uint8_t>& b, size_t off, uint32_t& x) {
  if (off + 4 > b.size()) return false;
  x = (uint32_t)b[off] |
      ((uint32_t)b[off+1] << 8) |
      ((uint32_t)b[off+2] << 16) |
      ((uint32_t)b[off+3] << 24);
  return true;
}
static inline void push_u64_le(vector<uint8_t>& b, uint64_t x) {
  for (int i = 0; i < 8; ++i) b.push_back((uint8_t)((x >> (8 * i)) & 0xFF));
}
static inline bool read_u64_le(const vector<uint8_t>& b, size_t off, uint64_t& x) {
  if (off + 8 > b.size()) return false;
  uint64_t v = 0;
  for (int i = 0; i < 8; ++i) v |= (uint64_t)b[off+i] << (8 * i);
  x = v;
  return true;
}

// ---------------- Window cursor BCW2 ----------------
struct WindowCursor {
  bool present = false;
  uint8_t flags = 0; // bit0=random_access
  uint8_t k0=0, kout=0, d=0;
  uint32_t numShards=0;
  uint64_t seed=0;
  uint32_t next_perm_pos=0;
  uint16_t window=0;
  uint16_t burst=0;

  struct LaneState {
    bool active=false;
    uint32_t perm_pos=0;
    uint8_t mode=0; // 0=konly, 1=expand

    // konly:
    uint64_t after=UINT64_MAX;

    // expand:
    uint64_t parent_anchor=UINT64_MAX; // UINT64_MAX means not started
    bool child_present=false;
    uint8_t L=0;
    uint64_t left_idx=0;
    uint64_t right_idx=0;
  };
  vector<LaneState> lanes;
};

static string make_cursor_bcw2(const WindowCursor& c) {
  vector<uint8_t> b;
  b.reserve(64 + c.lanes.size() * 64);
  b.push_back('B'); b.push_back('C'); b.push_back('W'); b.push_back('2');
  b.push_back(c.flags);
  b.push_back(c.k0); b.push_back(c.kout); b.push_back(c.d);
  push_u32_le(b, c.numShards);
  push_u64_le(b, c.seed);
  push_u32_le(b, c.next_perm_pos);
  push_u16_le(b, c.window);
  push_u16_le(b, c.burst);
  push_u16_le(b, (uint16_t)c.lanes.size());

  for (const auto& ln : c.lanes) {
    b.push_back(ln.active ? 1 : 0);
    if (!ln.active) continue;

    push_u32_le(b, ln.perm_pos);
    b.push_back(ln.mode);

    if (ln.mode == 0) {
      push_u64_le(b, ln.after);
    } else {
      push_u64_le(b, ln.parent_anchor);
      b.push_back(ln.child_present ? 1 : 0);
      if (ln.child_present) {
        b.push_back(ln.L);
        push_u64_le(b, ln.left_idx);
        push_u64_le(b, ln.right_idx);
      }
    }
  }
  return b64url_encode(b);
}

static bool parse_cursor_bcw2(const string& token, WindowCursor& c) {
  vector<uint8_t> b;
  if (!b64url_decode(token, b)) return false;
  if (b.size() < 4 + 1 + 3 + 4 + 8 + 4 + 2 + 2 + 2) return false;
  if (!(b[0]=='B' && b[1]=='C' && b[2]=='W' && b[3]=='2')) return false;

  c.flags = b[4];
  c.k0 = b[5]; c.kout = b[6]; c.d = b[7];
  if (!read_u32_le(b, 8, c.numShards)) return false;
  if (!read_u64_le(b, 12, c.seed)) return false;
  if (!read_u32_le(b, 20, c.next_perm_pos)) return false;

  size_t off = 24;
  if (!read_u16_le(b, off, c.window)) return false; off += 2;
  if (!read_u16_le(b, off, c.burst)) return false; off += 2;

  uint16_t lane_count=0;
  if (!read_u16_le(b, off, lane_count)) return false; off += 2;

  c.lanes.clear();
  c.lanes.resize(lane_count);

  for (uint16_t i=0; i<lane_count; ++i) {
    if (off + 1 > b.size()) return false;
    c.lanes[i].active = (b[off++] != 0);
    if (!c.lanes[i].active) continue;

    if (off + 4 + 1 > b.size()) return false;
    if (!read_u32_le(b, off, c.lanes[i].perm_pos)) return false; off += 4;
    c.lanes[i].mode = b[off++];

    if (c.lanes[i].mode == 0) {
      if (!read_u64_le(b, off, c.lanes[i].after)) return false; off += 8;
    } else {
      if (!read_u64_le(b, off, c.lanes[i].parent_anchor)) return false; off += 8;
      if (off + 1 > b.size()) return false;
      c.lanes[i].child_present = (b[off++] != 0);
      if (c.lanes[i].child_present) {
        if (off + 1 > b.size()) return false;
        c.lanes[i].L = b[off++];
        if (!read_u64_le(b, off, c.lanes[i].left_idx)) return false; off += 8;
        if (!read_u64_le(b, off, c.lanes[i].right_idx)) return false; off += 8;
      } else {
        c.lanes[i].L = 0;
        c.lanes[i].left_idx = 0;
        c.lanes[i].right_idx = 0;
      }
    }
  }

  c.present = true;
  return true;
}

// ---------------- CLI ----------------
struct Args {
  string shardsDir;
  string gcHistPath;

  string substring;
  bool substring_set=false;

  bool reverse_complement=false; // NEW

  int gcMinPct=0, gcMaxPct=100;
  uint64_t limit=200;
  int threads=4;

  int construct_k=-1;

  bool cursor_set=false;
  string cursor_token;

  bool random_access=false;
  bool ra_seed_set=false;
  uint64_t ra_seed=0;

  uint16_t window=16;
  uint16_t burst=1;

  uint32_t refill_chunk=256;
};

static void usage(const char* prog) {
  cerr << "Usage: " << prog
       << " --shards <dir> --gc-hist <json>"
       << " [--construct_k X]"
       << " [--substring <DNA>]"
       << " [--reverse_complement]"
       << " [--gc-min 0..100] [--gc-max 0..100]"
       << " [--limit N]"
       << " [--threads N]"
       << " [--window W] [--burst B]"
       << " [--cursor <BCW2...>]"
       << " [--random_access [--ra_seed U64]]\n";
}

static bool parse_args(int argc, char** argv, Args& a) {
  for (int i=1;i<argc;i++) {
    string s(argv[i]);
    if (s=="--shards" && i+1<argc) a.shardsDir=argv[++i];
    else if (s=="--gc-hist" && i+1<argc) a.gcHistPath=argv[++i];
    else if (s=="--construct_k" && i+1<argc) a.construct_k=stoi(argv[++i]);
    else if (s=="--substring" && i+1<argc) { a.substring=argv[++i]; a.substring_set=!a.substring.empty(); }
    else if (s=="--reverse_complement") a.reverse_complement=true; // NEW
    else if (s=="--gc-min" && i+1<argc) a.gcMinPct=stoi(argv[++i]);
    else if (s=="--gc-max" && i+1<argc) a.gcMaxPct=stoi(argv[++i]);
    else if (s=="--limit" && i+1<argc) a.limit=stoull(argv[++i]);
    else if (s=="--threads" && i+1<argc) a.threads=max(1, stoi(argv[++i]));
    else if (s=="--window" && i+1<argc) a.window=(uint16_t)max(1, stoi(argv[++i]));
    else if (s=="--burst" && i+1<argc) a.burst=(uint16_t)max(1, stoi(argv[++i]));
    else if (s=="--cursor" && i+1<argc) { a.cursor_set=true; a.cursor_token=argv[++i]; }
    else if (s=="--random_access") a.random_access=true;
    else if (s=="--ra_seed" && i+1<argc) { a.ra_seed_set=true; a.ra_seed=(uint64_t)stoull(argv[++i]); }
    else if (s=="--refill_chunk" && i+1<argc) a.refill_chunk=(uint32_t)max(16, stoi(argv[++i]));
    else { cerr << "Unknown arg: " << s << "\n"; return false; }
  }

  if (a.shardsDir.empty() || a.gcHistPath.empty()) return false;
  if (a.gcMinPct < 0 || a.gcMinPct > 100 || a.gcMaxPct < 0 || a.gcMaxPct > 100 || a.gcMinPct > a.gcMaxPct) {
    cerr << "GC range must satisfy 0<=gc-min<=gc-max<=100\n";
    return false;
  }
  if (a.limit < 1) return false;
  return true;
}

// ---------------- Lane runtime ----------------
struct LaneRuntime {
  bool active=false;
  uint32_t perm_pos=0;
  unsigned shardIdx=0;
  string shardPath;

  roaring64_bitmap_t* bm=nullptr;
  KbitHeader hdr;

  // Mode 0 (konly)
  uint64_t after=UINT64_MAX;

  // Mode 1 (expand) COMPLETE state:
  uint64_t parent_anchor=UINT64_MAX; // UINT64_MAX => not started
  bool child_present=false;
  uint8_t L=0;
  uint64_t left_idx=0;
  uint64_t right_idx=0;

  vector<uint64_t> buf;
  size_t buf_pos=0;

  void clear_buf(){ buf.clear(); buf_pos=0; }

  void free_all() {
    if (bm) { roaring64_bitmap_free(bm); bm=nullptr; }
    active=false;
    clear_buf();
  }
};

static inline bool leaf_ok(uint64_t vX, int kout,
                           int gcMinPct, int gcMaxPct,
                           bool substring_set, const vector<Pattern>& patterns) {
  if (!passes_gc_percent(vX, kout, gcMinPct, gcMaxPct)) return false;
  if (substring_set && !contains_sub(vX, patterns.data(), patterns.size())) return false;
  return true;
}


// Fill lane buffer by scanning lexicographically in-shard
static void refill_lane(LaneRuntime& lane,
                        int k0, int kout,
                        int gcMinPct, int gcMaxPct,
                        bool substring_set, const vector<Pattern>& patterns,
                        uint32_t refill_target,
                        const vector<uint64_t>& shard_starts,
                        const vector<uint64_t>& shard_ends)
{
  lane.clear_buf();
  if (!lane.active || !lane.bm) return;

  if (kout == k0) {
    uint64_t shardIdx = lane.shardIdx;
    if (shardIdx >= shard_starts.size() || shardIdx >= shard_ends.size()) {
      lane.active = false;
      return;
    }
    uint64_t start = shard_starts[shardIdx];
    uint64_t end = shard_ends[shardIdx];
    uint64_t v = (lane.after == UINT64_MAX) ? start : lane.after + 1;

    for (; v < end && lane.buf.size() < refill_target; ++v) {
      if (roaring64_bitmap_contains(lane.bm, v)) continue;
      if (!leaf_ok(v, kout, gcMinPct, gcMaxPct, substring_set, patterns)) continue;
      lane.buf.push_back(v);
    }

    if (v == end) {
      lane.active = false;
    } else {
      lane.after = v - 1;
    }
    return;
  }

  // construct_k>k0 expansion stream
  const int d = kout - k0;

  uint64_t shardIdx = lane.shardIdx;
  if (shardIdx >= shard_starts.size() || shardIdx >= shard_ends.size()) {
    lane.active = false;
    return;
  }
  const uint64_t start = shard_starts[shardIdx];
  const uint64_t end = shard_ends[shardIdx];

  while (lane.buf.size() < refill_target) {
    uint64_t parentB = 0;
    if (lane.parent_anchor == UINT64_MAX) {
      parentB = start;
    } else if (lane.child_present) {
      parentB = lane.parent_anchor;
    } else {
      parentB = lane.parent_anchor + 1;
    }

    while (parentB < end && roaring64_bitmap_contains(lane.bm, parentB)) parentB++;
    if (parentB >= end) { lane.active=false; break; }

    uint8_t Lcur;
    uint64_t li, ri;

    if (lane.child_present && lane.parent_anchor == parentB) {
      Lcur = lane.L;
      li = lane.left_idx;
      ri = lane.right_idx;
    } else {
      init_state_first(d, Lcur, li, ri);
    }

    bool exhausted_parent = false;
    while (!exhausted_parent && lane.buf.size() < refill_target) {
      uint64_t vX = make_value(parentB, k0, kout, (int)Lcur, li, ri);
      if (leaf_ok(vX, kout, gcMinPct, gcMaxPct, substring_set, patterns)) {
        lane.buf.push_back(vX);
      }
      if (!advance_state(d, Lcur, li, ri)) exhausted_parent = true;
    }

    if (lane.buf.size() >= refill_target) {
      lane.parent_anchor = parentB;
      lane.child_present = true;
      lane.L = Lcur;
      lane.left_idx = li;
      lane.right_idx = ri;
      break;
    }

    lane.parent_anchor = parentB;
    lane.child_present = false;
    lane.L = 0; lane.left_idx = 0; lane.right_idx = 0;
  }
}

// ---------------- Main ----------------
int main(int argc, char** argv) {
  ios::sync_with_stdio(false);
  cin.tie(nullptr);

  Args args;
  if (!parse_args(argc, argv, args)) { usage(argv[0]); return 1; }

  // We support k0 in {16,17,18} shard sets on-disk.
  // Expansion (construct_k > k0) is ONLY allowed when k0==18.
  // If the caller requests construct_k > 18, we force base shards to 18-mers.
  const int requested_kout = (args.construct_k > 0) ? args.construct_k : -1;

  unsigned numShards=0;
  vector<string> shardFiles;
  uint64_t k_from_index=0;
  uint64_t total_bits_index=0;
  vector<uint64_t> shard_starts;
  vector<uint64_t> shard_ends;
  if (!read_index(args.shardsDir, numShards, shardFiles, k_from_index,
                  total_bits_index, shard_starts, shard_ends)) {
    cerr << "Failed to read " << args.shardsDir << "/index.json\n";
    return 1;
  }

  const int k_index = (int)k_from_index;
  if (k_index <= 0 || k_index > 32) {
    cerr << "Error: invalid k in index.json: " << k_index << "\n";
    return 1;
  }

  // Decide the effective base-k (k0) and output-k (kout).
  int k0 = k_index;
  int kout = (requested_kout > 0) ? requested_kout : k0;
  if (kout > 32) { cerr << "Error: construct_k>32 not supported in uint64 encoding\n"; return 1; }

  // Enforce: only 18-mer base allows expansion.
  // For kout>18, we require k0==18. If the user pointed us at shards_16/shards_17, fail loudly.
  // (The web server will route to shards_18 automatically for kout>18.)
  if (kout > 18 && k0 != 18) {
    cerr << "Error: construct_k>18 expansion is only supported from k=18 base shards. "
         << "Got base k=" << k0 << ".\n";
    return 1;
  }

  // Enforce: no expansion for k0<18 at all.
  if (k0 < 18 && kout != k0) {
    cerr << "Error: expansion is disabled for k=" << k0 << ". "
         << "Use construct_k=" << k0 << ".\n";
    return 1;
  }
  int k_from_hist=0;
  vector<vector<uint64_t>> gc_hists;
  auto t_hist0 = Clock::now();
  if (!load_gc_hist_json(args.gcHistPath, k_from_hist, gc_hists)) {
    cerr << "Failed to load gc histogram json: " << args.gcHistPath << "\n";
    return 1;
  }
  auto t_hist1 = Clock::now();

  if (k_from_hist != k0) {
    cerr << "GC hist k (" << k_from_hist << ") != index k (" << k0 << ")\n";
    return 1;
  }

  if (shard_starts.size() != numShards || shard_ends.size() != numShards) {
    uint64_t total_bits = total_bits_index ? total_bits_index : (1ULL << (2 * k0));
    const uint64_t width = (total_bits + (uint64_t)numShards - 1) / (uint64_t)numShards;
    shard_starts.assign(numShards, 0);
    shard_ends.assign(numShards, 0);
    for (unsigned i=0;i<numShards;i++) {
      uint64_t start = (uint64_t)i * width;
      uint64_t end = min<uint64_t>(total_bits, start + width);
      shard_starts[i] = start;
      shard_ends[i] = end;
    }
  }

  // substring patterns for kout
  vector<Pattern> patterns;

  auto append_patterns_for = [&](const string& sub) -> bool {
    const int m = (int)sub.size();
    if (m == 0) return true;
    if (m > kout) { cerr << "substring longer than output k\n"; return false; }

    uint64_t sub_bits=0;
    for (char c : sub) {
      int dd = base4_digit(c);
      if (dd < 0) { cerr << "Invalid base in substring\n"; return false; }
      sub_bits = (sub_bits<<2) | (uint64_t)dd;
    }

    uint64_t base_mask = (m >= 32) ? ~0ULL : ((1ULL << (2*m)) - 1ULL);
    for (int pos=0; pos<=kout-m; ++pos) {
      int shift = 2 * (kout - m - pos);
      patterns.push_back({ base_mask << shift, sub_bits << shift });
    }
    return true;
  };

  if (args.substring_set) {
    patterns.reserve((size_t)(2 * (kout - (int)args.substring.size() + 1)));

    if (!append_patterns_for(args.substring)) return 1;

    if (args.reverse_complement) {
      string rc = revcomp_string(args.substring);
      // If substring is palindromic (rc == substring), avoid duplicating the same patterns.
      if (rc != args.substring) {
        if (!append_patterns_for(rc)) return 1;
      }
    }
  }

  // permutation seed
  uint64_t seed = 0;
  if (args.random_access) {
    if (args.ra_seed_set) seed = args.ra_seed;
    else {
      std::random_device rd;
      uint64_t s1 = ((uint64_t)rd() << 32) ^ (uint64_t)rd();
      uint64_t s2 = ((uint64_t)rd() << 32) ^ (uint64_t)rd();
      seed = splitmix64(s1 ^ (s2 << 1));
    }
    if (seed == 0) seed = 1;
  }

  vector<uint32_t> perm;
  if (args.random_access) perm = build_perm((uint32_t)numShards, seed);
  else { perm.resize(numShards); for (uint32_t i=0;i<numShards;i++) perm[i]=i; }

  auto shard_from_permpos = [&](uint32_t perm_pos)->unsigned { return (unsigned)perm[perm_pos]; };

  // Cursor init
  uint32_t next_perm_pos = 0;
  vector<WindowCursor::LaneState> lane_states;

  if (args.cursor_set) {
    WindowCursor in;
    if (!parse_cursor_bcw2(args.cursor_token, in)) {
      cerr << "Error: expected BCW2 cursor\n";
      return 1;
    }
    if (in.numShards != (uint32_t)numShards) { cerr << "Error: cursor mismatch numShards\n"; return 1; }
    if (in.k0 != (uint8_t)k0 || in.kout != (uint8_t)kout) { cerr << "Error: cursor mismatch k\n"; return 1; }
    if (in.window != args.window) { cerr << "Error: cursor window mismatch\n"; return 1; }
    if (in.burst != args.burst) { cerr << "Error: cursor burst mismatch\n"; return 1; }

    bool cursor_ra = (in.flags & 0x1) != 0;
    if (cursor_ra != args.random_access) { cerr << "Error: cursor random_access mismatch\n"; return 1; }

    if (args.random_access) {
      // cursor wins seed if user passed another
      seed = in.seed;
      if (seed == 0) seed = 1;
      perm = build_perm((uint32_t)numShards, seed);
    }

    next_perm_pos = in.next_perm_pos;
    lane_states = in.lanes;
  } else {
    next_perm_pos = 0;
    lane_states.resize(args.window);
    for (auto& ls : lane_states) ls.active=false;
  }

  // Runtime lanes
  vector<LaneRuntime> lanes(args.window);

  auto load_lane_from_state = [&](int i, const WindowCursor::LaneState& st)->bool {
    lanes[i].free_all();
    lanes[i].active = false;
    if (!st.active) return false;
    if (st.perm_pos >= numShards) return false;

    unsigned shardIdx = shard_from_permpos(st.perm_pos);

    lanes[i].perm_pos = st.perm_pos;
    lanes[i].shardIdx = shardIdx;
    lanes[i].shardPath = args.shardsDir + "/" + shardFiles[shardIdx];

    KbitHeader hdr;
    roaring64_bitmap_t* bm = load_kbit_portable(lanes[i].shardPath, hdr);
    if (!bm) return false;

    lanes[i].bm = bm;
    lanes[i].hdr = hdr;
    lanes[i].active = true;

    lanes[i].clear_buf();

    if (kout == k0) {
      lanes[i].after = st.after;
    } else {
      lanes[i].parent_anchor = st.parent_anchor;
      lanes[i].child_present = st.child_present;
      lanes[i].L = st.L;
      lanes[i].left_idx = st.left_idx;
      lanes[i].right_idx = st.right_idx;
    }
    return true;
  };

  // Load lanes from cursor state
  for (int i=0;i<(int)args.window;i++) {
    if (i < (int)lane_states.size() && lane_states[i].active) {
      if (!load_lane_from_state(i, lane_states[i])) lane_states[i].active = false;
    }
  }

  auto try_fill_empty_lane = [&](int i)->bool {
    if (lanes[i].active) return true;

    while (next_perm_pos < numShards) {
      uint32_t ppos = next_perm_pos++;
      unsigned shardIdx = shard_from_permpos(ppos);

      lanes[i].perm_pos = ppos;
      lanes[i].shardIdx = shardIdx;
      lanes[i].shardPath = args.shardsDir + "/" + shardFiles[shardIdx];

      KbitHeader hdr;
      roaring64_bitmap_t* bm = load_kbit_portable(lanes[i].shardPath, hdr);
      if (!bm) return false;

      lanes[i].bm = bm;
      lanes[i].hdr = hdr;
      lanes[i].clear_buf();
      lanes[i].active = true;

      if (kout == k0) {
        lanes[i].after = UINT64_MAX;
      } else {
        lanes[i].parent_anchor = UINT64_MAX;
        lanes[i].child_present = false;
        lanes[i].L = 0; lanes[i].left_idx = 0; lanes[i].right_idx = 0;
      }
      return true;
    }
    return false;
  };

  // Fill any empty lanes initially
  for (int i=0;i<(int)args.window;i++) if (!lanes[i].active) (void)try_fill_empty_lane(i);

  const uint64_t need = args.limit + 1;
  vector<uint64_t> out_vals;
  out_vals.reserve((size_t)need);

  double scan_sec_total=0.0;
  uint64_t shards_loaded=0;
  for (auto& ln : lanes) if (ln.bm) shards_loaded++;

  auto t_scan0 = Clock::now();

  while (out_vals.size() < need) {
    bool any_active=false;
    for (auto& ln : lanes) if (ln.active) { any_active=true; break; }
    if (!any_active) break;

    // Parallel refill of empty buffers
    {
      atomic<int> idx(0);
      int T = min(args.threads, (int)args.window);
      vector<thread> pool;
      pool.reserve((size_t)T);

      for (int t=0;t<T;t++) {
        pool.emplace_back([&](){
          while (true) {
            int i = idx.fetch_add(1);
            if (i >= (int)args.window) break;
            if (!lanes[i].active) continue;
            if (lanes[i].buf_pos < lanes[i].buf.size()) continue;

            refill_lane(lanes[i], k0, kout, args.gcMinPct, args.gcMaxPct,
                        args.substring_set, patterns, args.refill_chunk,
                        shard_starts, shard_ends);

            if (!lanes[i].active) {
              lanes[i].free_all();
              (void)try_fill_empty_lane(i);
              if (lanes[i].bm) shards_loaded++;
            }
          }
        });
      }
      for (auto& th : pool) th.join();
    }

    // Round-robin emission
    bool emitted_any=false;
    for (int i=0;i<(int)args.window && out_vals.size() < need; ++i) {
      if (!lanes[i].active) continue;

      uint16_t took=0;
      while (took < args.burst && out_vals.size() < need) {
        if (lanes[i].buf_pos >= lanes[i].buf.size()) break;
        uint64_t v = lanes[i].buf[lanes[i].buf_pos++];

        if (kout == k0) lanes[i].after = v;
        // expand mode state is already updated during refill (anchor/child state)

        out_vals.push_back(v);
        took++;
        emitted_any=true;
      }
    }

    if (!emitted_any) {
      bool still=false;
      for (auto& ln : lanes) if (ln.active) { still=true; break; }
      if (!still) break;
    }
  }

  auto t_scan1 = Clock::now();
  scan_sec_total = chrono::duration_cast<Sec>(t_scan1 - t_scan0).count();

  bool hasMore = false;
  if (out_vals.size() > args.limit) { out_vals.resize((size_t)args.limit); hasMore=true; }
  else {
    for (auto& ln : lanes) {
      if (!ln.active) continue;
      if (ln.buf_pos < ln.buf.size()) { hasMore=true; break; }
      if (kout > k0) { hasMore=true; break; }
    }
    if (!hasMore && next_perm_pos < numShards) hasMore=true;
  }

  // Build next cursor
  string cursorStr;
  if (hasMore) {
    WindowCursor outc;
    outc.flags = (args.random_access ? 0x1 : 0);
    outc.k0=(uint8_t)k0; outc.kout=(uint8_t)kout; outc.d=(uint8_t)(kout - k0);
    outc.numShards=(uint32_t)numShards;
    outc.seed=seed;
    outc.next_perm_pos=next_perm_pos;
    outc.window=args.window;
    outc.burst=args.burst;
    outc.lanes.resize(args.window);

    for (int i=0;i<(int)args.window;i++) {
      auto& st = outc.lanes[i];
      if (!lanes[i].active) { st.active=false; continue; }
      st.active=true;
      st.perm_pos=lanes[i].perm_pos;
      if (kout == k0) {
        st.mode=0;
        st.after=lanes[i].after;
      } else {
        st.mode=1;
        st.parent_anchor=lanes[i].parent_anchor;
        st.child_present=lanes[i].child_present;
        st.L=lanes[i].L;
        st.left_idx=lanes[i].left_idx;
        st.right_idx=lanes[i].right_idx;
      }
    }

    cursorStr = make_cursor_bcw2(outc);
  } else {
    cursorStr = "";
  }

  // Emit
  cout << "__META__\t" << cursorStr << "\t" << (hasMore ? "1" : "0") << "\t" << out_vals.size() << "\t" << kout << "\n";
  for (uint64_t v : out_vals) cout << decode_kmer(v, kout) << "\n";

  // Cleanup
  for (auto& ln : lanes) ln.free_all();

  long pk = peak_rss_kb();
  cerr << fixed << setprecision(6);
  cerr << "[INFO] Shards dir          : " << args.shardsDir << "\n";
  cerr << "[INFO] GC hist             : " << args.gcHistPath << "\n";
  cerr << "[INFO] Threads             : " << args.threads << "\n";
  cerr << "[INFO] Limit               : " << args.limit << "\n";
  cerr << "[INFO] window / burst      : " << args.window << " / " << args.burst << "\n";
  cerr << "[INFO] refill_chunk        : " << args.refill_chunk << "\n";
  cerr << "[INFO] k0 / kout           : " << k0 << " / " << kout << "\n";
  cerr << "[INFO] Random access       : " << (args.random_access ? "yes" : "no") << "\n";
  if (args.random_access) cerr << "[INFO] RA seed             : " << seed << "\n";
  cerr << "[INFO] GC% range           : " << args.gcMinPct << "-" << args.gcMaxPct << "\n";
  cerr << "[INFO] Substring           : " << (args.substring_set ? args.substring : "(none)") << "\n";
  cerr << "[INFO] Reverse complement  : " << (args.reverse_complement ? "yes" : "no") << "\n";
  cerr << "[INFO] Returned            : " << out_vals.size() << "\n";
  cerr << "[INFO] Has more            : " << (hasMore ? "yes" : "no") << "\n";
  cerr << "[INFO] Next cursor         : " << (cursorStr.empty() ? "(none)" : cursorStr) << "\n";
  cerr << "[INFO] Shards loaded        : " << shards_loaded << "\n";
  cerr << "[INFO] GC hist load time    : " << chrono::duration_cast<Sec>(t_hist1 - t_hist0).count() << " s\n";
  cerr << "[INFO] Scan time            : " << scan_sec_total << " s\n";
  cerr << "[INFO] Peak RSS             : " << pk << " KB (" << (pk/1024.0) << " MB)\n";

  return 0;
}
