#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cinttypes>
#include <cmath>
#include <string>
#include <vector>
#include <roaring/roaring64.h>

// -----------------------------------------------------------------------
// KBITv1 header (64 bytes, little-endian)
// -----------------------------------------------------------------------
#pragma pack(push, 1)
struct KbitHeader {
    char     magic[8];      // 0-7:   "KBITv1\0"
    uint64_t total_bits;    // 8-15:  4^k
    uint64_t ones;          // 16-23: count of set bits
    uint64_t k;             // 24-31: k-mer length
    uint64_t seed;          // 32-39: RNG seed
    uint64_t flags;         // 40-47: 1=dense, 2=roaring
    uint64_t payload_len;   // 48-55: payload bytes
    char     reserved[8];   // 56-63: zeros
};
#pragma pack(pop)

static_assert(sizeof(KbitHeader) == 64, "KbitHeader must be exactly 64 bytes");

// -----------------------------------------------------------------------
// K-mer encoding: base-4, A=0, C=1, G=2, T=3, leftmost = most significant
// -----------------------------------------------------------------------
static uint64_t encode_kmer(const std::string &kmer) {
    uint64_t val = 0;
    for (char c : kmer) {
        val <<= 2;
        switch (c) {
            case 'A': case 'a': val |= 0; break;
            case 'C': case 'c': val |= 1; break;
            case 'G': case 'g': val |= 2; break;
            case 'T': case 't': val |= 3; break;
            default:
                fprintf(stderr, "ERROR: invalid nucleotide '%c' in '%s'\n", c, kmer.c_str());
                return UINT64_MAX;
        }
    }
    return val;
}

static std::string decode_kmer(uint64_t val, int k) {
    std::string s(k, ' ');
    for (int i = k - 1; i >= 0; --i) {
        switch (val & 3) {
            case 0: s[i] = 'A'; break;
            case 1: s[i] = 'C'; break;
            case 2: s[i] = 'G'; break;
            case 3: s[i] = 'T'; break;
        }
        val >>= 2;
    }
    return s;
}

// -----------------------------------------------------------------------
// Write a single KBITv1 file using a roaring64 bitmap
// -----------------------------------------------------------------------
static void write_kbit(const std::string &path, int k,
                       roaring64_bitmap_t *bitmap,
                       uint64_t seed) {
    uint64_t ones = roaring64_bitmap_get_cardinality(bitmap);
    uint64_t total_bits = (uint64_t)1 << (2 * k);   // 4^k, safe for k <= 31

    // Serialize roaring bitmap to portable byte array
    size_t payload_len = roaring64_bitmap_portable_size_in_bytes(bitmap);
    std::vector<char> serialized(payload_len, 0);
    size_t written = roaring64_bitmap_portable_serialize(bitmap, serialized.data());
    if (written != payload_len) {
        fprintf(stderr, "ERROR: serialized size mismatch (%zu vs %zu)\n", written, payload_len);
        exit(1);
    }

    // Populate header
    KbitHeader hdr;
    std::memset(&hdr, 0, sizeof(hdr));
    std::memcpy(hdr.magic, "KBITv1\0", 8);
    hdr.total_bits  = total_bits;
    hdr.ones        = ones;
    hdr.k           = static_cast<uint64_t>(k);
    hdr.seed        = seed;
    hdr.flags       = 2;          // roaring
    hdr.payload_len = payload_len;

    // Write file (all fields in little-endian)
    FILE *f = std::fopen(path.c_str(), "wb");
    if (!f) {
        std::perror("fopen");
        std::exit(1);
    }

    auto write_le64 = [&](uint64_t v) {
        std::fwrite(&v, sizeof(v), 1, f);
    };

    std::fwrite(hdr.magic,     8, 1, f);
    write_le64(hdr.total_bits);
    write_le64(hdr.ones);
    write_le64(hdr.k);
    write_le64(hdr.seed);
    write_le64(hdr.flags);
    write_le64(hdr.payload_len);
    std::fwrite(hdr.reserved,  8, 1, f);

    std::fwrite(serialized.data(), 1, payload_len, f);
    std::fclose(f);

    std::printf("  -> %s  (k=%d, ones=%" PRIu64 ", payload=%zu B)\n",
                path.c_str(), k, ones, payload_len);
}

// -----------------------------------------------------------------------
// Print and optionally add a k-mer to a bitmap
// -----------------------------------------------------------------------
static void report_kmer(const std::string &kmer, int k, bool is_present,
                        roaring64_bitmap_t *bitmap) {
    uint64_t code = encode_kmer(kmer);
    if (code >= (uint64_t)1 << (2 * k)) {
        // K-mer is longer than k — cannot exist in this bitmap
        std::printf("  [k=%d] %s %s -> %" PRIu64 "  (out of range for k=%d)\n",
                    k, is_present ? "PRESENT" : "ABSENT",
                    kmer.c_str(), code, k);
        return;
    }
    if (is_present) {
        roaring64_bitmap_add(bitmap, code);
    }
    std::printf("  [k=%d] %s %s -> %" PRIu64 "\n",
                k, is_present ? "PRESENT" : "ABSENT",
                kmer.c_str(), code);
}

// -----------------------------------------------------------------------
// Print a JSON truth entry to stdout
// -----------------------------------------------------------------------
static void print_truth_json(int k, const std::string &fixture_path,
                             const std::vector<std::string> &present,
                             const std::vector<std::string> &absent) {
    auto esc = [](const std::string &s) -> std::string {
        std::string out;
        out.reserve(s.size() + 2);
        out += '"';
        for (char c : s) out += c;
        out += '"';
        return out;
    };

    auto join = [&](const std::vector<std::string> &v) -> std::string {
        std::string out;
        for (size_t i = 0; i < v.size(); ++i) {
            if (i) out += ", ";
            out += esc(v[i]);
        }
        return out;
    };

    std::printf("  {\n");
    std::printf("    \"k\": %d,\n", k);
    std::printf("    \"fixture_path\": %s,\n", esc(fixture_path).c_str());
    std::printf("    \"present_kmers\": [%s],\n", join(present).c_str());
    std::printf("    \"absent_kmers\": [%s],\n", join(absent).c_str());
    std::printf("    \"sequence_encoding\": {\"A\": 0, \"C\": 1, \"G\": 2, \"T\": 3}\n");
    std::printf("  }");
}

// =======================================================================
int main() {
    const uint64_t seed = 42;

    // Ensure output directory exists
    std::system("mkdir -p qc/tmp");

    std::printf("\n=== Generating synthetic KBITv1 fixtures ===\n\n");
    std::printf("Truth JSON:\n[\n");

    // -------------------------------------------------------------------
    // k = 16
    // -------------------------------------------------------------------
    {
        std::printf("\n--- k=16 ---\n");
        int k = 16;
        roaring64_bitmap_t *bm = roaring64_bitmap_create();

        std::vector<std::string> present = {
            "ACGTACGTACGTACGT",
            "TTTTCCCCAAAAGGGG",
            "GATTACAGATTACAGA",
            "CCCCGGGGAAAATTTT",
            "AAAAACCCCCCCCCCCC"
        };

        std::vector<std::string> absent = {
            "AAAAAAAAAAAAAAAA",
            "CCCCCCCCCCCCCCCC",
            "GGGGGGGGGGGGGGGG",
            "TATATATATATATATA",
            "GCTCCCTGTAAGACCCCA"
        };

        for (const auto &kmer : present) report_kmer(kmer, k, true,  bm);
        for (const auto &kmer : absent)  report_kmer(kmer, k, false, bm);

        std::string path = "qc/tmp/k16_fixture.kbit";
        write_kbit(path, k, bm, seed);
        roaring64_bitmap_free(bm);

        print_truth_json(k, path, present, absent);
        std::printf(",\n");
    }

    // -------------------------------------------------------------------
    // k = 17
    // -------------------------------------------------------------------
    {
        std::printf("\n--- k=17 ---\n");
        int k = 17;
        roaring64_bitmap_t *bm = roaring64_bitmap_create();

        std::vector<std::string> present = {
            "ACGTACGTACGTACGTA",
            "TTTTCCCCAAAAGGGGA",
            "GATTACAGATTACAGAT"
        };

        std::vector<std::string> absent = {
            "AAAAAAAAAAAAAAAAA",
            "CCCCCCCCCCCCCCCCC",
            "GGGGGGGGGGGGGGGGG"
        };

        for (const auto &kmer : present) report_kmer(kmer, k, true,  bm);
        for (const auto &kmer : absent)  report_kmer(kmer, k, false, bm);

        std::string path = "qc/tmp/k17_fixture.kbit";
        write_kbit(path, k, bm, seed);
        roaring64_bitmap_free(bm);

        print_truth_json(k, path, present, absent);
        std::printf(",\n");
    }

    // -------------------------------------------------------------------
    // k = 18
    // -------------------------------------------------------------------
    {
        std::printf("\n--- k=18 ---\n");
        int k = 18;
        roaring64_bitmap_t *bm = roaring64_bitmap_create();

        std::vector<std::string> present = {
            "ACGTACGTACGTACGTAC",
            "TTTTCCCCAAAAGGGGAA",
            "GATTACAGATTACAGATT"
        };

        std::vector<std::string> absent = {
            "AAAAAAAAAAAAAAAAAA",
            "CCCCCCCCCCCCCCCCCC",
            "GGGGGGGGGGGGGGGGGG",
            "GCTCCCTGTAAGACCCCA",
            "AAGCTTGTTGCTCGTCGG"
        };

        for (const auto &kmer : present) report_kmer(kmer, k, true,  bm);
        for (const auto &kmer : absent)  report_kmer(kmer, k, false, bm);

        std::string path = "qc/tmp/k18_fixture.kbit";
        write_kbit(path, k, bm, seed);
        roaring64_bitmap_free(bm);

        print_truth_json(k, path, present, absent);
        std::printf("\n");
    }

    std::printf("]\n\n=== Done ===\n");

    return 0;
}
