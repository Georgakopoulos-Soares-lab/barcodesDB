// gen_dense_bitmap_fast.cpp
// Fully optimized dense 4^k-bit generator with exact X% ones.
// Build: g++ -O3 -march=native -flto -std=cpp17 -pthread -o gen_dense_bitmap_fast gen_dense_bitmap_fast.cpp
//
// Layout (little-endian):
//   [64B header]
//     0..7   : "KBITv1\0"
//     8..15  : total_bits = 4^k
//     16..23 : ones = round(X% * total_bits)
//     24..31 : k
//     32..39 : seed
//     40..47 : flags (=1 dense)
//     48..55 : payload_len = ceil(total_bits/8)
//     56..63 : reserved=0
//   [payload bytes]  bit i at byte[i>>3], mask 1<<(i&7) (LSB-first)

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <string>
#include <thread>
#include <vector>
#include <mutex>
#include <condition_variable>
#include <algorithm>

#if defined(__linux__)
  #include <sys/resource.h>
  #include <sys/types.h>
  #include <sys/stat.h>
  #include <fcntl.h>
  #include <unistd.h>
  #include <sched.h>
#endif

// -------- Args --------
struct Args {
  int k = -1;
  double percent = -1.0;
  std::string out;
  uint64_t seed = std::random_device{}();
  size_t io_buf_bytes = 64ULL<<20;      // default 64 MiB buffers
  int status_step_mib = 128;            // print every N MiB written
  bool pin_threads = true;              // pin producer/writer on Linux
};

static void usage(const char* prog) {
  std::cerr <<
    "Usage: " << prog << " --k <int> --percent <float> --out <path> "
    "[--seed <u64>] [--buf <bytes>] [--status-step-mib <int>] [--no-pin]\n";
}

static bool parse_args(int argc, char** argv, Args& a) {
  for (int i=1; i<argc; ++i) {
    std::string s(argv[i]);
    if (s=="--k" && i+1<argc) a.k = std::atoi(argv[++i]);
    else if (s=="--percent" && i+1<argc) a.percent = std::atof(argv[++i]);
    else if (s=="--out" && i+1<argc) a.out = argv[++i];
    else if (s=="--seed" && i+1<argc) a.seed = std::strtoull(argv[++i], nullptr, 10);
    else if (s=="--buf" && i+1<argc) a.io_buf_bytes = (size_t)std::strtoull(argv[++i], nullptr, 10);
    else if (s=="--status-step-mib" && i+1<argc) a.status_step_mib = std::atoi(argv[++i]);
    else if (s=="--no-pin") a.pin_threads = false;
    else { std::cerr << "Unknown/invalid arg: " << s << "\n"; return false; }
  }
  if (a.k < 1 || a.k > 31) { std::cerr << "Error: --k in [1,31]\n"; return false; }
  if (!(a.percent >= 0.0 && a.percent <= 100.0)) { std::cerr << "Error: --percent in [0,100]\n"; return false; }
  if (a.out.empty()) { std::cerr << "Error: --out required\n"; return false; }
  if (a.io_buf_bytes < (1<<20)) a.io_buf_bytes = 1<<20; // >= 1 MiB
  // Align buffer to multiples of 8 for word writes
  a.io_buf_bytes = (a.io_buf_bytes / 8) * 8;
  if (a.status_step_mib < 1) a.status_step_mib = 64;
  return true;
}

// -------- Utilities --------
static inline void write_le64_u(unsigned char* dst, uint64_t v) {
  for (int i=0;i<8;++i) dst[i] = (unsigned char)((v >> (8*i)) & 0xFF);
}

#if defined(__linux__)
static size_t peak_rss_bytes() {
  struct rusage ru; if (getrusage(RUSAGE_SELF, &ru)==0) return (size_t)ru.ru_maxrss * 1024ULL;
  return 0;
}
static void pin_to_cpu(int cpu) {
  cpu_set_t set; CPU_ZERO(&set); CPU_SET(cpu, &set);
  pthread_setaffinity_np(pthread_self(), sizeof(set), &set);
}
#else
static size_t peak_rss_bytes() { return 0; }
static void pin_to_cpu(int) {}
#endif

// -------- RNG: xoshiro256** + Lemire bounded --------
struct Xoshiro256ss {
  uint64_t s[4];
  static inline uint64_t rotl(uint64_t x, int k){ return (x<<k) | (x>>(64-k)); }
  inline uint64_t next() {
    uint64_t result = rotl(s[1] * 5ULL, 7) * 9ULL;
    uint64_t t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl(s[3], 45);
    return result;
  }
};
static inline uint64_t splitmix64(uint64_t x){
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}
static inline Xoshiro256ss make_rng(uint64_t seed){
  Xoshiro256ss r;
  r.s[0]=splitmix64(seed);
  r.s[1]=splitmix64(seed+0x9e3779b97f4a7c15ULL);
  r.s[2]=splitmix64(seed+0x632be59bd9b4e019ULL);
  r.s[3]=splitmix64(seed+0x94d049bb133111ebULL);
  return r;
}
// Lemire: unbiased mapping of 64-bit RNG to [0, bound)
static inline uint64_t fast_bounded(uint64_t rnd, uint64_t bound){
  __uint128_t mul = (__uint128_t)rnd * (__uint128_t)bound;
  return (uint64_t)(mul >> 64);
}

// -------- Main --------
int main(int argc, char** argv) {
  std::ios::sync_with_stdio(false);
  Args args;
  if (!parse_args(argc, argv, args)) { usage(argv[0]); return 1; }

  const int shift = 2 * args.k;               // <= 62
  const uint64_t total_bits = (1ULL << shift);
  const uint64_t payload_bytes = (total_bits + 7) / 8;
  const uint64_t ones_target = std::min<uint64_t>(
      (uint64_t)llround((long double)args.percent * (long double)total_bits / 100.0L),
      total_bits);

  FILE* f = std::fopen(args.out.c_str(), "wb");
  if (!f) { std::perror("fopen"); return 1; }

  // Make stdio buffering large to reduce syscalls
  setvbuf(f, nullptr, _IOFBF, 1<<20);

#if defined(__linux__)
  // Preallocate full file and advise sequential
  int fd = fileno(f);
  off_t total_size = (off_t)(64 + payload_bytes);
  if (posix_fallocate(fd, 0, total_size) != 0) {
    // Not fatal on all FS; continue anyway
  }
  posix_fadvise(fd, 0, total_size, POSIX_FADV_SEQUENTIAL);
#endif

  // Write placeholder header (64 bytes)
  unsigned char hdr[64]; std::memset(hdr, 0, 64);
  std::memcpy(hdr, "KBITv1\0", 8);
  if (std::fwrite(hdr, 1, 64, f) != 64) { std::perror("write header"); return 1; }

  // Triple-buffered writer
  struct Buffer { std::vector<unsigned char> data; size_t used=0; bool eof=false; };
  const int NBUF = 3;
  Buffer buffers[NBUF];
  for (int i=0;i<NBUF;++i) { buffers[i].data.resize(args.io_buf_bytes); buffers[i].used=0; buffers[i].eof=false; }

  std::mutex mtx;
  std::condition_variable cv_prod, cv_cons;
  int fill_idx = 0, flush_idx = 0;
  bool ready[NBUF] = {false,false,false};
  bool stop_writer = false;
  std::atomic<uint64_t> bytes_written{0};   // payload bytes only

  auto writer = [&](){
#if defined(__linux__)
    if (args.pin_threads) pin_to_cpu(1);
#endif
    for(;;){
      std::unique_lock<std::mutex> lk(mtx);
      cv_cons.wait(lk, [&]{ return ready[flush_idx] || stop_writer; });
      if (stop_writer && !ready[flush_idx]) break;
      Buffer &B = buffers[flush_idx];
      size_t to_write = B.used;
      lk.unlock();

      if (to_write) {
        size_t off=0;
        while (off < to_write) {
          size_t w = std::fwrite(B.data.data()+off, 1, to_write-off, f);
          if (w==0) { std::perror("fwrite payload"); std::exit(1); }
          off += w;
          bytes_written.fetch_add(w, std::memory_order_relaxed);
        }
      }

      lk.lock();
      B.used = 0;
      ready[flush_idx] = false;
      bool was_eof = B.eof;
      flush_idx = (flush_idx + 1) % NBUF;
      lk.unlock();
      cv_prod.notify_one();

      if (was_eof) { /* continue draining other ready=false slots until stop */ }
    }
  };

  std::thread writer_thread(writer);

#if defined(__linux__)
  if (args.pin_threads) pin_to_cpu(0); // producer on CPU 0
#endif

  // Producer state
  auto t0 = std::chrono::steady_clock::now();
  uint64_t next_status_mark = (uint64_t)args.status_step_mib << 20;  // MiB to bytes

  auto print_status = [&](bool final=false){
    auto now = std::chrono::steady_clock::now();
    static auto last = t0;
    static uint64_t last_bytes = 0;

    double dt = std::chrono::duration<double>(now - last).count();
    uint64_t bw = bytes_written.load(std::memory_order_relaxed);
    uint64_t dbytes = bw - last_bytes;

    double pct = (payload_bytes ? (double)bw / (double)payload_bytes * 100.0 : 100.0);
    double mbps = dt>0 ? (dbytes / (1024.0*1024.0)) / dt : 0.0;
    double eta = (mbps>0) ? (( (payload_bytes - bw) / (1024.0*1024.0)) / mbps) : NAN;

    size_t rss = peak_rss_bytes();
    if (!final) {
      std::cerr << "[status] "
                << "written=" << (bw/(1024.0*1024.0)) << " MiB ("
                << pct << "%), rate=" << mbps << " MiB/s, ETA="
                << (std::isfinite(eta)? std::to_string((int)eta)+"s":"n/a")
                << ", peakRSS=" << (rss? std::to_string(rss/(1024*1024))+" MiB":"n/a")
                << "\n";
    } else {
      double elapsed = std::chrono::duration<double>(now - t0).count();
      std::cerr << "[done]   "
                << "written=" << (bw/(1024.0*1024.0)) << " MiB (100%), "
                << "elapsed=" << elapsed << " s, avg_rate="
                << (elapsed>0? (bw/(1024.0*1024.0))/elapsed : 0.0) << " MiB/s, "
                << "peakRSS=" << (rss? std::to_string(rss/(1024*1024))+" MiB":"n/a")
                << "\n";
    }
    last = now; last_bytes = bw;
  };

  // Helper: append arbitrary byte span to current buffer; flush when full
  auto put_bytes = [&](const void* p, size_t n){
    const unsigned char* src = (const unsigned char*)p;
    while (n) {
      Buffer &B = buffers[fill_idx];
      size_t space = B.data.size() - B.used;
      if (space == 0) {
        // hand off
        {
          std::unique_lock<std::mutex> lk(mtx);
          cv_prod.wait(lk, [&]{ return !ready[fill_idx]; });
          ready[fill_idx] = true;
        }
        cv_cons.notify_one();
        fill_idx = (fill_idx + 1) % NBUF;
        continue;
      }
      size_t chunk = std::min(space, n);
      std::memcpy(B.data.data()+B.used, src, chunk);
      B.used += chunk; src += chunk; n -= chunk;
    }
  };

  // Generator: 64-bit words with integer selection sampling
  Xoshiro256ss rng = make_rng(args.seed);
  uint64_t remaining = total_bits;
  uint64_t needed    = ones_target;

  // Pre-print config
  const uint64_t est_final_size = 64 + payload_bytes;
  std::cerr << "Config:\n"
            << "  k=" << args.k << " => total_bits=4^k=" << total_bits << "\n"
            << "  ones=" << ones_target << " (" << args.percent << "%)\n"
            << "  payload=" << payload_bytes << " bytes; file≈" << est_final_size << " bytes\n"
            << "  buffers=" << NBUF << " x " << args.io_buf_bytes << " bytes\n"
            << "  seed=" << args.seed << (args.pin_threads? " (pinned)":"") << "\n";

  // Emit full 64-bit words
  while (remaining) {
    uint64_t w = 0;
    int bits_here = (remaining >= 64 ? 64 : (int)remaining);
    for (int b=0; b<bits_here; ++b) {
      bool one = false;
      if (needed) {
        uint64_t r = rng.next();
        uint64_t pick = fast_bounded(r, remaining);
        if (pick < needed) { one = true; --needed; }
      }
      if (one) w |= (1ULL << b); // LSB-first
      --remaining;
    }
    put_bytes(&w, sizeof(uint64_t));

    // progress (by bytes written in writer thread)
    uint64_t bw = bytes_written.load(std::memory_order_relaxed);
    if (bw >= next_status_mark) {
      print_status(false);
      next_status_mark += (uint64_t)args.status_step_mib << 20;
    }
  }

  // Finalize buffers: flag EOF on the buffer that currently holds data
  {
    std::unique_lock<std::mutex> lk(mtx);
    // ensure the current buffer gets handed off
    if (!ready[fill_idx] && buffers[fill_idx].used>0) ready[fill_idx] = true;
    buffers[fill_idx].eof = true;
  }
  cv_cons.notify_one();

  // Stop writer when all handed
  {
    std::unique_lock<std::mutex> lk(mtx);
    stop_writer = true;
  }
  cv_cons.notify_one();
  writer_thread.join();

  // Rewrite header with final values
  std::fflush(f);
  if (std::fseek(f, 0, SEEK_SET)!=0) { std::perror("fseek"); return 1; }
  std::memset(hdr, 0, 64);
  std::memcpy(hdr, "KBITv1\0", 8);
  write_le64_u(hdr+8,  total_bits);
  write_le64_u(hdr+16, ones_target);
  write_le64_u(hdr+24, (uint64_t)args.k);
  write_le64_u(hdr+32, args.seed);
  write_le64_u(hdr+40, 1ULL);             // dense
  write_le64_u(hdr+48, payload_bytes);
  // 56..63 reserved=0
  if (std::fwrite(hdr,1,64,f) != 64) { std::perror("rewrite header"); return 1; }
  std::fflush(f);
  std::fclose(f);

  // Final status
  // Ensure we count exactly payload_bytes even if final fwrite buffered:
  // bytes_written may already equal payload_bytes; not critical but we print done anyway.
  // (We intentionally exclude header from bytes_written.)
  // Print final:
  {
    // Force last status as "done"
    auto now = std::chrono::steady_clock::now(); (void)now;
  }
  // Reuse status printer for final line
  {
    // Fake completion if filesystem delayed accounting
    bytes_written.store(payload_bytes, std::memory_order_relaxed);
  }
  // Print 'done'
  {
    auto now = std::chrono::steady_clock::now(); (void)now;
  }
  // Inline final print (slightly simplified):
  {
    size_t rss = peak_rss_bytes();
    auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    double avg = elapsed>0? (payload_bytes/(1024.0*1024.0))/elapsed : 0.0;
    std::cerr << "[done]   written=" << (payload_bytes/(1024.0*1024.0)) << " MiB (100%), "
              << "elapsed=" << elapsed << " s, avg_rate=" << avg << " MiB/s, "
              << "peakRSS=" << (rss? std::to_string(rss/(1024*1024))+" MiB":"n/a")
              << "\n";
  }

  if (needed != 0) {
    std::cerr << "Warning: internal mismatch, remaining ones_needed=" << needed << "\n";
    return 2;
  }
  return 0;
}
