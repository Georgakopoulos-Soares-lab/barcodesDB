// query_kmer_bitmap_single_fast.cpp
// Faster single-threaded random-access k-mer existence queries for a KBITv1 roaring64 payload.

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include <roaring/roaring64.h>
#include "motif_filters.hpp"

struct Args {
  // Sharded mode
  std::string shards;

  // Legacy single-bitmap mode (optional)
  std::string bitmap;

  // Optional explicit k (when caller wants to hard-require a specific length)
  int k = -1;

  std::string kmers;  // optional
  std::string out;    // optional

  int threads = 4;
  int min_hamming_distance = 0;

  // Motif filter options (all off by default)
  MotifFilterOptions motif_opts;
};

static void usage(const char* prog) {
  std::cerr << "Usage: " << prog
            << " --shards <dir> [--k 16|17|18] [--kmers <file>] [--out <file>]"
            << " [--threads N]\n"
            << "  Motif filters (all off by default):\n"
            << "    --motif-mode <off|flag|exclude>\n"
            << "    --filter-homopolymers  --max-homopolymer <int> (default 4)\n"
            << "    --filter-low-complexity --min-shannon-entropy <float> (default 1.5)\n"
            << "    --filter-dinucleotide-repeats\n"
            << "    --filter-trinucleotide-repeats\n"
            << "    --filter-tetranucleotide-repeats\n"
            << "    --filter-restriction-sites\n"
            << "    --filter-functional-motifs\n";
}

static bool parse_args(int argc, char** argv, Args& a) {
  for (int i = 1; i < argc; ++i) {
    std::string s(argv[i]);
    if (s == "--shards" && i + 1 < argc) a.shards = argv[++i];
    else if (s == "--bitmap" && i + 1 < argc) a.bitmap = argv[++i];
    else if (s == "--k" && i + 1 < argc) a.k = std::atoi(argv[++i]);
    else if (s == "--kmers" && i + 1 < argc) a.kmers = argv[++i];
    else if (s == "--out" && i + 1 < argc) a.out = argv[++i];
    else if (s == "--min-hamming-distance" && i + 1 < argc) a.min_hamming_distance = std::max(0, std::atoi(argv[++i]));
    // Motif filter parameters
    else if (s == "--motif-mode" && i + 1 < argc) a.motif_opts.motif_mode = argv[++i];
    else if (s == "--filter-homopolymers") a.motif_opts.filter_homopolymers = true;
    else if (s == "--max-homopolymer" && i + 1 < argc) a.motif_opts.max_homopolymer = std::max(1, std::atoi(argv[++i]));
    else if (s == "--filter-low-complexity") a.motif_opts.filter_low_complexity = true;
    else if (s == "--min-shannon-entropy" && i + 1 < argc) a.motif_opts.min_shannon_entropy = std::atof(argv[++i]);
    else if (s == "--filter-dinucleotide-repeats") a.motif_opts.filter_dinucleotide_repeats = true;
    else if (s == "--filter-trinucleotide-repeats") a.motif_opts.filter_trinucleotide_repeats = true;
    else if (s == "--filter-tetranucleotide-repeats") a.motif_opts.filter_tetranucleotide_repeats = true;
    else if (s == "--filter-restriction-sites") a.motif_opts.filter_restriction_sites = true;
    else if (s == "--filter-functional-motifs") a.motif_opts.filter_functional_motifs = true;
    else { std::cerr << "Unknown/invalid arg: " << s << "\n"; return false; }
  }

  if (a.shards.empty() && a.bitmap.empty()) {
    std::cerr << "Error: --shards is required (or --bitmap for legacy mode)\n";
    return false;
  }
  if (a.k != -1 && !(a.k == 16 || a.k == 17 || a.k == 18)) {
    std::cerr << "Error: --k must be 16, 17, or 18\n";
    return false;
  }
  return true;
}
static inline uint64_t read_le64_u(const unsigned char* p) {
  uint64_t v = 0;
  for (int i = 0; i < 8; ++i) v |= (uint64_t)p[i] << (8 * i);
  return v;
}

struct Header {
  uint64_t total_bits = 0;
  uint64_t ones = 0;
  uint64_t k = 0;
  uint64_t seed = 0;
  uint64_t flags = 0;
  uint64_t payload_len = 0;
};

static bool read_header(std::ifstream& f, Header& h) {
  unsigned char hdr[64];
  f.read(reinterpret_cast<char*>(hdr), 64);
  if (!f) return false;
  if (std::memcmp(hdr, "KBITv1\0", 8) != 0) {
    std::cerr << "Error: bad magic (not a KBITv1 file)\n";
    return false;
  }
  h.total_bits  = read_le64_u(hdr + 8);
  h.ones        = read_le64_u(hdr + 16);
  h.k           = read_le64_u(hdr + 24);
  h.seed        = read_le64_u(hdr + 32);
  h.flags       = read_le64_u(hdr + 40);
  h.payload_len = read_le64_u(hdr + 48);
  return true;
}

static roaring64_bitmap_t* load_kbit_portable(const std::string& path, Header& H, std::vector<char>& buf) {
  std::ifstream in(path, std::ios::binary);
  if (!in) { std::perror(("open shard: " + path).c_str()); return nullptr; }
  if (!read_header(in, H)) return nullptr;
  if (H.flags != 2) {
    std::cerr << "Error: expected roaring payload (flags=2) in " << path << "\n";
    return nullptr;
  }
  buf.resize(H.payload_len);
  in.read(buf.data(), (std::streamsize)buf.size());
  if ((uint64_t)in.gcount() != H.payload_len) {
    std::cerr << "Error: truncated payload in " << path << "\n";
    return nullptr;
  }
  roaring64_bitmap_t* rbm = roaring64_bitmap_portable_deserialize_safe(buf.data(), buf.size());
  if (!rbm) { std::cerr << "Error: deserialization failed for " << path << "\n"; return nullptr; }
  return rbm;
}

static roaring64_bitmap_t* load_bitmap_file(const std::string& path, Header& H, std::vector<char>& buf) {
  std::ifstream in(path, std::ios::binary);
  if (!in) { std::perror(("open bitmap: " + path).c_str()); return nullptr; }
  if (!read_header(in, H)) return nullptr;
  if (H.flags != 2) {
    std::cerr << "Error: expected roaring payload (flags=2) in " << path << "\n";
    return nullptr;
  }
  buf.resize(H.payload_len);
  in.read(buf.data(), (std::streamsize)buf.size());
  if ((uint64_t)in.gcount() != H.payload_len) {
    std::cerr << "Error: truncated payload in " << path << "\n";
    return nullptr;
  }
  roaring64_bitmap_t* rbm =
      roaring64_bitmap_portable_deserialize_safe(buf.data(), buf.size());
  if (!rbm) { std::cerr << "Error: deserialization failed for " << path << "\n"; return nullptr; }
  return rbm;
}

struct ShardInfo {
  uint64_t start = 0;
  uint64_t end = 0;
  std::string file;
};

static bool read_index_shards(const std::string& dir, unsigned& numShards, uint64_t& k_out,
                              std::vector<ShardInfo>& shards) {
  std::ifstream in(dir + "/index.json");
  if (!in) return false;

  std::string line;
  numShards = 0;
  k_out = 0;
  shards.clear();

  auto parse_u64_field = [&](const std::string& s, const std::string& key, uint64_t& out)->bool {
    auto pos = s.find(key);
    if (pos == std::string::npos) return false;
    pos = s.find(':', pos);
    if (pos == std::string::npos) return false;
    pos++;
    while (pos < s.size() && (s[pos] == ' ' || s[pos] == '\t')) pos++;
    size_t end = s.find_first_of(",}", pos);
    if (end == std::string::npos || end <= pos) return false;
    out = std::stoull(s.substr(pos, end - pos));
    return true;
  };

  while (std::getline(in, line)) {
    if (line.find("\"num_shards\"") != std::string::npos) {
      auto p = line.find(':');
      if (p != std::string::npos) numShards = (unsigned)std::stoul(line.substr(p+1));
    }
    if (line.find("\"k\"") != std::string::npos && line.find("\"seed\"") == std::string::npos) {
      auto p = line.find(':');
      if (p != std::string::npos) k_out = (uint64_t)std::stoull(line.substr(p+1));
    }

    auto fpos = line.find("\"file\"");
    if (fpos != std::string::npos) {
      uint64_t start=0, end=0;
      bool has_start = parse_u64_field(line, "\"start\"", start);
      bool has_end = parse_u64_field(line, "\"end\"", end);

      auto colon = line.find(':', fpos);
      if (colon == std::string::npos) continue;
      auto s1 = line.find('"', colon);
      if (s1 == std::string::npos) continue;
      auto s2 = line.find('"', s1 + 1);
      if (s2 == std::string::npos) continue;

      ShardInfo si;
      si.file = line.substr(s1 + 1, s2 - (s1 + 1));
      si.start = has_start ? start : 0;
      si.end = has_end ? end : 0;
      shards.push_back(si);
    }
  }

  if (numShards == 0) numShards = (unsigned)shards.size();
  if (numShards == 0 || k_out == 0) return false;
  if (shards.size() != numShards) {
    shards.clear();
    shards.reserve(numShards);
    for (unsigned i = 0; i < numShards; ++i) {
      ShardInfo si;
      char tmp[64];
      std::snprintf(tmp, sizeof(tmp), "shard_%04u.kbit", i);
      si.file = tmp;
      si.start = 0;
      si.end = 0;
      shards.push_back(si);
    }
  }
  return true;
}

static int find_shard(const std::vector<ShardInfo>& shards, uint64_t idx) {
  if (shards.empty()) return -1;
  size_t lo = 0;
  size_t hi = shards.size();
  while (lo < hi) {
    size_t mid = lo + (hi - lo) / 2;
    const auto& s = shards[mid];
    if (idx < s.start) hi = mid;
    else if (idx >= s.end) lo = mid + 1;
    else return (int)mid;
  }
  return -1;
}

// Fast base->2bit LUT. 0xFF means invalid.
static inline const uint8_t* base_lut() {
  static uint8_t lut[256];
  static bool inited = false;
  if (!inited) {
    std::fill(std::begin(lut), std::end(lut), (uint8_t)0xFF);
    lut[(unsigned)'A'] = lut[(unsigned)'a'] = 0;
    lut[(unsigned)'C'] = lut[(unsigned)'c'] = 1;
    lut[(unsigned)'G'] = lut[(unsigned)'g'] = 2;
    lut[(unsigned)'T'] = lut[(unsigned)'t'] = 3;
    inited = true;
  }
  return lut;
}

// Encode exactly k chars from sv; returns false if length!=k or invalid chars.
static inline bool encode_kmer_sv(std::string_view sv, int k, uint64_t& out) {
  if ((int)sv.size() != k) return false;
  const uint8_t* lut = base_lut();
  uint64_t v = 0;
  const unsigned char* p = (const unsigned char*)sv.data();
  for (int i = 0; i < k; ++i) {
    uint8_t d = lut[p[i]];
    if (d == 0xFF) return false;
    v = (v << 2) | (uint64_t)d;
  }
  out = v;
  return true;
}

// Very fast line reader (handles \n and optional \r).
struct FastLineReader {
  static constexpr size_t BUFSZ = 1 << 20; // 1 MiB
  FILE* f = nullptr;
  std::vector<char> buf;
  size_t pos = 0, len = 0;
  std::string line;

  explicit FastLineReader(FILE* fp) : f(fp), buf(BUFSZ) { line.reserve(256); }

  bool refill() {
    len = std::fread(buf.data(), 1, buf.size(), f);
    pos = 0;
    return len > 0;
  }

  // Reads next line into `line` (without trailing \r or \n). Returns false at EOF.
  bool next() {
    line.clear();

    while (true) {
      if (pos >= len) {
        if (!refill()) return !line.empty(); // last line without newline
      }

      // Scan for newline in current buffer
      size_t start = pos;
      while (pos < len && buf[pos] != '\n') ++pos;

      line.append(buf.data() + start, pos - start);

      if (pos < len && buf[pos] == '\n') { // consume newline
        ++pos;
        break;
      }
      // else: buffer ended without newline, continue refilling
    }

    // Trim trailing '\r' if present (Windows line endings)
    if (!line.empty() && line.back() == '\r') line.pop_back();
    return true;
  }
};

static int hamming_distance(const std::string& a, const std::string& b) {
  if (a.size() != b.size()) return -1;
  int dist = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    if (a[i] != b[i]) ++dist;
  }
  return dist;
}

int main(int argc, char** argv) {
  std::ios::sync_with_stdio(false);

  Args args;
  if (!parse_args(argc, argv, args)) { usage(argv[0]); return 1; }

  Header H;
  std::vector<char> buf;
  roaring64_bitmap_t* rbm = nullptr;

  unsigned numShards = 0;
  uint64_t k_from_index = 0;
  std::vector<ShardInfo> shards;

  if (!args.shards.empty()) {
    if (!read_index_shards(args.shards, numShards, k_from_index, shards)) {
      std::cerr << "Error: failed to read shards index: " << args.shards << "/index.json\n";
      return 2;
    }
  } else {
    rbm = load_bitmap_file(args.bitmap, H, buf);
    if (!rbm) return 2;
  }

  int k_fixed = (args.shards.empty()) ? (int)H.k : (int)k_from_index;
  if (args.k != -1 && k_fixed != args.k) {
    std::cerr << "Error: bitmap header k=" << k_fixed << " does not match requested --k " << args.k << "\n";
    if (rbm) roaring64_bitmap_free(rbm);
    return 2;
  }
  if (!(k_fixed == 16 || k_fixed == 17 || k_fixed == 18)) {
    std::cerr << "Error: unsupported bitmap k=" << k_fixed << " (expected 16/17/18)\n";
    if (rbm) roaring64_bitmap_free(rbm);
    return 2;
  }

  // --- Input FILE*
  FILE* fin = nullptr;
  if (!args.kmers.empty()) {
    fin = std::fopen(args.kmers.c_str(), "rb");
    if (!fin) { std::perror("open kmers"); roaring64_bitmap_free(rbm); return 1; }
  } else {
    // stdin
    fin = stdin;
    std::fprintf(stderr, "Enter k-mers (one per line, Ctrl+D to finish):\n");
  }

  // --- Output FILE*
  FILE* fout = nullptr;
  if (!args.out.empty()) {
    fout = std::fopen(args.out.c_str(), "wb");
    if (!fout) { std::perror("open out"); if (fin!=stdin) std::fclose(fin); roaring64_bitmap_free(rbm); return 1; }
  } else {
    fout = stdout;
  }

  // Bigger stdio buffers (often helps a lot)
  static char outbuf[1 << 20];
  std::setvbuf(fout, outbuf, _IOFBF, sizeof(outbuf));

  FastLineReader r(fin);

  std::vector<std::string> kmers;
  std::vector<uint64_t> kmer_vals;
  kmers.reserve(1 << 20);
  kmer_vals.reserve(1 << 20);

  while (r.next()) {
    if (r.line.empty()) continue;
    std::string_view sv(r.line);
    uint64_t idx = 0;
    bool ok = encode_kmer_sv(sv, k_fixed, idx);
    if (!ok) {
      std::cerr << "Error: encountered k-mer of length " << sv.size()
                << " but this query expects k=" << k_fixed << "\n";
      if (fout != stdout) std::fclose(fout);
      if (fin != stdin) std::fclose(fin);
      if (rbm) roaring64_bitmap_free(rbm);
      return 3;
    }
    kmers.emplace_back(sv);
    kmer_vals.push_back(idx);
  }

  std::vector<char> hits(kmers.size(), '0');

  if (!args.shards.empty()) {
    if (shards.empty()) {
      std::cerr << "Error: no shards listed in index.json\n";
      if (fout != stdout) std::fclose(fout);
      if (fin != stdin) std::fclose(fin);
      return 2;
    }

    for (auto& s : shards) {
      if (s.end <= s.start) {
        std::cerr << "Error: shard ranges missing in index.json (start/end)\n";
        if (fout != stdout) std::fclose(fout);
        if (fin != stdin) std::fclose(fin);
        return 2;
      }
    }

    std::vector<std::vector<size_t>> shard_to_indices(shards.size());
    for (size_t i = 0; i < kmer_vals.size(); ++i) {
      int sid = find_shard(shards, kmer_vals[i]);
      if (sid < 0) {
        std::cerr << "Error: k-mer index out of shard ranges\n";
        if (fout != stdout) std::fclose(fout);
        if (fin != stdin) std::fclose(fin);
        return 2;
      }
      shard_to_indices[(size_t)sid].push_back(i);
    }

    std::atomic<size_t> next_shard(0);
    int thread_count = std::min<int>(args.threads, (int)shards.size());
    std::vector<std::thread> pool;
    pool.reserve((size_t)thread_count);

    for (int t = 0; t < thread_count; ++t) {
      pool.emplace_back([&]() {
        std::vector<char> local_buf;
        Header local_h;
        while (true) {
          size_t sid = next_shard.fetch_add(1);
          if (sid >= shards.size()) break;
          if (shard_to_indices[sid].empty()) continue;

          const std::string shard_path = args.shards + "/" + shards[sid].file;
          roaring64_bitmap_t* sbm = load_kbit_portable(shard_path, local_h, local_buf);
          if (!sbm) continue;

          for (size_t idx_pos : shard_to_indices[sid]) {
            uint64_t val = kmer_vals[idx_pos];
            hits[idx_pos] = roaring64_bitmap_contains(sbm, val) ? '1' : '0';
          }
          roaring64_bitmap_free(sbm);
        }
      });
    }
    for (auto& th : pool) th.join();
  } else {
    for (size_t i = 0; i < kmer_vals.size(); ++i) {
      hits[i] = roaring64_bitmap_contains(rbm, kmer_vals[i]) ? '1' : '0';
    }
  }

  // Hamming distance computation (only when requested)
  std::vector<int> nearest_hamming(kmers.size(), -1);
  std::vector<char> passes_hamming(kmers.size(), '0');
  if (args.min_hamming_distance > 0) {
    // Collect indices of found k-mers
    std::vector<size_t> found_idx;
    for (size_t i = 0; i < hits.size(); ++i) {
      if (hits[i] == '1') found_idx.push_back(i);
    }
    // For each found k-mer, find minimum Hamming distance to any other found k-mer
    for (size_t fi = 0; fi < found_idx.size(); ++fi) {
      size_t i = found_idx[fi];
      int min_dist = k_fixed + 1; // larger than max possible distance
      for (size_t fj = 0; fj < found_idx.size(); ++fj) {
        if (fi == fj) continue;
        size_t j = found_idx[fj];
        int d = hamming_distance(kmers[i], kmers[j]);
        if (d >= 0 && d < min_dist) min_dist = d;
      }
      nearest_hamming[i] = (min_dist <= k_fixed) ? min_dist : -1;
      passes_hamming[i] = (min_dist >= args.min_hamming_distance) ? '1' : '0';
    }
  }

  // Reuse output line buffer (max 18 + tab + 3 + tab + 1 + \n + motif metadata)
  std::string out_line;
  out_line.reserve(128);
  bool motif_active = (args.motif_opts.motif_mode != "off");
  for (size_t i = 0; i < kmers.size(); ++i) {
    out_line.assign(kmers[i]);
    out_line.push_back('\t');
    out_line.push_back(hits[i]);
    if (args.min_hamming_distance > 0) {
      out_line.push_back('\t');
      if (nearest_hamming[i] >= 0) {
        out_line.append(std::to_string(nearest_hamming[i]));
      }
      out_line.push_back('\t');
      out_line.push_back(passes_hamming[i]);
    }
    // Motif filter metadata (appended as additional tab-delimited columns)
    if (motif_active) {
      MotifFilterResult mf = evaluate_motif_filters(kmers[i], args.motif_opts);
      out_line.push_back('\t');
      out_line.push_back(mf.passes ? '1' : '0');                     // motif_passes
      out_line.push_back('\t');
      out_line.append(std::to_string(mf.hits.size()));                // motif_hit_count
      out_line.push_back('\t');
      out_line.append(motif_hits_to_string(mf.hits));                 // motif_hits detail
    }
    out_line.push_back('\n');
    std::fwrite(out_line.data(), 1, out_line.size(), fout);
  }

  if (fout != stdout) std::fclose(fout);
  if (fin != stdin) std::fclose(fin);
  if (rbm) roaring64_bitmap_free(rbm);
  return 0;
}
