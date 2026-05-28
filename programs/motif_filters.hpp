// motif_filters.hpp — Shared sequence-level motif filter for barcodesDB
//
// Provides optional, lightweight checks on DNA barcode sequences for short
// motifs that may be undesirable in some experimental contexts:
//   - homopolymer runs
//   - low-complexity / low-entropy sequences
//   - dinucleotide repeats
//   - common restriction enzyme recognition sites
//   - short functional-like motifs (polyA signals, TATA-like, splice-like, etc.)
//
// These filters operate directly on candidate barcode strings and can either
// flag or exclude sequences depending on user settings. They are intended as
// practical safeguards, NOT as organism-specific functional annotation.
//
// All filters are OFF by default. No heavy dependencies.

#ifndef MOTIF_FILTERS_HPP
#define MOTIF_FILTERS_HPP

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

// ---------------- Data structures ----------------

/// A single motif hit detected in a sequence.
struct MotifHit {
    std::string category;      // e.g. "homopolymer", "restriction_site", "functional_motif"
    std::string motif;         // e.g. "AAAAA", "EcoRI:GAATTC", "polyA_signal:AATAAA"
    size_t      position = 0;  // 0-based start index in the original sequence
    std::string description;   // human-readable, e.g. "Homopolymer run exceeds threshold"
};

/// User-configurable filter options.
struct MotifFilterOptions {
    // Homopolymer filter
    bool filter_homopolymers = false;
    int  max_homopolymer     = 4;  // runs strictly longer than this fail

    // Low-complexity filter (Shannon entropy over A/C/G/T)
    bool   filter_low_complexity  = false;
    double min_shannon_entropy    = 1.5;  // fail if H < threshold

    // Dinucleotide repeat filter
    bool filter_dinucleotide_repeats = false;

    // Trinucleotide repeat filter (e.g., CAGCAGCAG, ATCATCATC)
    bool filter_trinucleotide_repeats = false;

    // Tetranucleotide repeat filter (e.g., GATAGATAGATA, CACTCACTCACT)
    bool filter_tetranucleotide_repeats = false;

    // Restriction-site filter
    bool filter_restriction_sites = false;

    // Functional-motif filter
    bool filter_functional_motifs = false;

    // Global mode: "off" (no-op), "flag" (report but don't exclude),
    //              "exclude" (skip sequences that fail)
    std::string motif_mode = "off";
};

/// Result of evaluating motif filters on a single sequence.
struct MotifFilterResult {
    bool passes = true;           // false if any enabled filter fails (exclude mode)
    std::vector<MotifHit> hits;   // all detected hits (even if passes==true in flag mode)
};

// ---------------- Helper: uppercase ----------------
inline std::string uppercase_dna(std::string_view seq) {
    std::string out;
    out.reserve(seq.size());
    for (char c : seq) out.push_back((char)std::toupper((unsigned char)c));
    return out;
}

// ---------------- Homopolymer runs ----------------
/// Finds runs longer than max_run. Returns hits for each run that exceeds.
inline std::vector<MotifHit> find_homopolymer_runs(std::string_view seq, int max_run) {
    std::vector<MotifHit> hits;
    if (seq.empty() || max_run < 1) return hits;

    size_t run_start = 0;
    char   run_base  = seq[0];
    for (size_t i = 1; i <= seq.size(); ++i) {
        if (i == seq.size() || seq[i] != run_base) {
            size_t run_len = i - run_start;
            if ((int)run_len > max_run) {
                MotifHit h;
                h.category    = "homopolymer";
                h.motif       = std::string(seq.substr(run_start, run_len));
                h.position    = run_start;
                h.description = "Homopolymer run exceeds threshold (" +
                                std::to_string(run_len) + " > " + std::to_string(max_run) + ")";
                hits.push_back(std::move(h));
            }
            if (i < seq.size()) {
                run_start = i;
                run_base  = seq[i];
            }
        }
    }
    return hits;
}

// ---------------- Shannon entropy ----------------
/// Compute Shannon entropy H = -sum(p_i * log2(p_i)) over A/C/G/T only.
/// Returns 0.0 for empty string.
inline double shannon_entropy(std::string_view seq) {
    if (seq.empty()) return 0.0;
    int counts[4] = {0, 0, 0, 0};
    int total = 0;
    for (char c : seq) {
        switch (c) {
            case 'A': case 'a': counts[0]++; total++; break;
            case 'C': case 'c': counts[1]++; total++; break;
            case 'G': case 'g': counts[2]++; total++; break;
            case 'T': case 't': counts[3]++; total++; break;
            default: break; // skip N and others (preserve current invalid-char behavior)
        }
    }
    if (total == 0) return 0.0;

    double H = 0.0;
    for (int i = 0; i < 4; ++i) {
        if (counts[i] == 0) continue;
        double p = (double)counts[i] / (double)total;
        H -= p * std::log2(p);
    }
    return H;
}

/// Check whether sequence entropy is below threshold.
inline bool has_low_entropy(std::string_view seq, double min_H) {
    return shannon_entropy(seq) < min_H;
}

// ---------------- Dinucleotide repeats ----------------
/// Detect dinucleotide repeat tracts of length >= 6 (3+ consecutive repeats).
inline std::vector<MotifHit> find_dinucleotide_repeats(std::string_view seq) {
    std::vector<MotifHit> hits;
    if (seq.size() < 6) return hits;

    for (size_t i = 0; i + 5 < seq.size(); ++i) {
        char b1 = seq[i];
        char b2 = seq[i + 1];
        // Skip if not A/C/G/T
        if (b1 != 'A' && b1 != 'C' && b1 != 'G' && b1 != 'T') continue;
        if (b2 != 'A' && b2 != 'C' && b2 != 'G' && b2 != 'T') continue;
        if (b1 == b2) continue; // homopolymer, not dinucleotide

        // Check that the pattern (b1,b2) repeats at least 3 times consecutively
        if (seq[i + 2] == b1 && seq[i + 3] == b2 &&
            seq[i + 4] == b1 && seq[i + 5] == b2) {
            // Extend to find full tract length
            size_t j = i + 6;
            while (j + 1 < seq.size() && seq[j] == b1 && seq[j + 1] == b2)
                j += 2;
            size_t tract_len = j - i;

            MotifHit h;
            h.category    = "dinucleotide_repeat";
            h.motif       = std::string(seq.substr(i, tract_len));
            h.position    = i;
            h.description = "Dinucleotide repeat tract (" + std::to_string(tract_len) + " bp)";
            hits.push_back(std::move(h));
            i = j - 1; // skip past this tract
        }
    }
    return hits;
}

// ---------------- Trinucleotide repeats ----------------
/// Detect trinucleotide repeat tracts of length >= 9 (3+ consecutive repeats).
/// Important for triplet expansion disorders and microsatellite instability.
inline std::vector<MotifHit> find_trinucleotide_repeats(std::string_view seq) {
    std::vector<MotifHit> hits;
    if (seq.size() < 9) return hits;

    for (size_t i = 0; i + 8 < seq.size(); ++i) {
        char b1 = seq[i];
        char b2 = seq[i + 1];
        char b3 = seq[i + 2];
        if (b1 != 'A' && b1 != 'C' && b1 != 'G' && b1 != 'T') continue;
        if (b2 != 'A' && b2 != 'C' && b2 != 'G' && b2 != 'T') continue;
        if (b3 != 'A' && b3 != 'C' && b3 != 'G' && b3 != 'T') continue;
        // Skip homopolymers and dinucleotide repeats (handled by other filters)
        if (b1 == b2 && b2 == b3) continue;
        if (b1 == b2 || b2 == b3) continue;

        // Check that the pattern (b1,b2,b3) repeats at least 3 times consecutively
        if (seq[i + 3] == b1 && seq[i + 4] == b2 && seq[i + 5] == b3 &&
            seq[i + 6] == b1 && seq[i + 7] == b2 && seq[i + 8] == b3) {
            // Extend to find full tract length
            size_t j = i + 9;
            while (j + 2 < seq.size() && seq[j] == b1 && seq[j + 1] == b2 && seq[j + 2] == b3)
                j += 3;
            size_t tract_len = j - i;

            MotifHit h;
            h.category    = "trinucleotide_repeat";
            h.motif       = std::string(seq.substr(i, tract_len));
            h.position    = i;
            h.description = "Trinucleotide repeat tract (" + std::to_string(tract_len) + " bp)";
            hits.push_back(std::move(h));
            i = j - 1;
        }
    }
    return hits;
}

// ---------------- Tetranucleotide repeats ----------------
/// Detect tetranucleotide repeat tracts of length >= 12 (3+ consecutive repeats).
/// Common in microsatellite markers and forensic STR loci.
inline std::vector<MotifHit> find_tetranucleotide_repeats(std::string_view seq) {
    std::vector<MotifHit> hits;
    if (seq.size() < 12) return hits;

    for (size_t i = 0; i + 11 < seq.size(); ++i) {
        char b1 = seq[i];
        char b2 = seq[i + 1];
        char b3 = seq[i + 2];
        char b4 = seq[i + 3];
        if (b1 != 'A' && b1 != 'C' && b1 != 'G' && b1 != 'T') continue;
        if (b2 != 'A' && b2 != 'C' && b2 != 'G' && b2 != 'T') continue;
        if (b3 != 'A' && b3 != 'C' && b3 != 'G' && b3 != 'T') continue;
        if (b4 != 'A' && b4 != 'C' && b4 != 'G' && b4 != 'T') continue;
        // Skip homopolymer/di/tri-nucleotide repeats
        if (b1 == b2 || b2 == b3 || b3 == b4) continue;

        // Check that the pattern (b1,b2,b3,b4) repeats at least 3 times consecutively
        if (seq[i + 4] == b1 && seq[i + 5] == b2 && seq[i + 6] == b3 && seq[i + 7] == b4 &&
            seq[i + 8] == b1 && seq[i + 9] == b2 && seq[i + 10] == b3 && seq[i + 11] == b4) {
            // Extend to find full tract length
            size_t j = i + 12;
            while (j + 3 < seq.size() && seq[j] == b1 && seq[j + 1] == b2 && seq[j + 2] == b3 && seq[j + 3] == b4)
                j += 4;
            size_t tract_len = j - i;

            MotifHit h;
            h.category    = "tetranucleotide_repeat";
            h.motif       = std::string(seq.substr(i, tract_len));
            h.position    = i;
            h.description = "Tetranucleotide repeat tract (" + std::to_string(tract_len) + " bp)";
            hits.push_back(std::move(h));
            i = j - 1;
        }
    }
    return hits;
}

// ---------------- Literal motif search ----------------
/// Find all occurrences of literal motifs in seq. Case-sensitive (expect uppercase).
inline std::vector<MotifHit> find_literal_motifs(std::string_view seq,
                                                  const std::vector<std::pair<std::string, std::string>>& motif_list,
                                                  const std::string& category) {
    std::vector<MotifHit> hits;
    for (const auto& [name, pattern] : motif_list) {
        if (pattern.empty()) continue;
        size_t pos = 0;
        while ((pos = seq.find(pattern, pos)) != std::string::npos) {
            MotifHit h;
            h.category    = category;
            h.motif       = name + ":" + pattern;
            h.position    = pos;
            h.description = "";
            hits.push_back(std::move(h));
            pos++; // overlapping matches
        }
    }
    return hits;
}

// ---------------- Restriction enzyme sites ----------------
/// Built-in list of common restriction enzyme recognition motifs.
/// Both forward and reverse-complement forms are checked.
inline const std::vector<std::pair<std::string, std::string>>& restriction_site_list() {
    // name, motif (forward strand)
    static const std::vector<std::pair<std::string, std::string>> sites = {
        // 6+ bp cutters
        {"EcoRI",  "GAATTC"},
        {"BamHI",  "GGATCC"},
        {"HindIII","AAGCTT"},
        {"NotI",   "GCGGCCGC"},
        {"XhoI",   "CTCGAG"},
        {"XbaI",   "TCTAGA"},
        {"SpeI",   "ACTAGT"},
        {"NheI",   "GCTAGC"},
        {"PstI",   "CTGCAG"},
        {"KpnI",   "GGTACC"},
        {"SacI",   "GAGCTC"},
        {"SalI",   "GTCGAC"},
        {"SmaI",   "CCCGGG"},
        {"MluI",   "ACGCGT"},
        {"AgeI",   "ACCGGT"},
        {"BglII",  "AGATCT"},
        {"AatII",  "GACGTC"},
        {"AccI",   "GTATAC"},
        {"AflII",  "CTTAAG"},
        {"AscI",   "GGCGCGCC"},
        {"AvrII",  "CCTAGG"},
        {"BclI",   "TGATCA"},
        {"BsiWI",  "CGTACG"},
        {"BspHI",  "TCATGA"},
        {"BsrGI",  "TGTACA"},
        {"BssHII", "GCGCGC"},
        {"ClaI",   "ATCGAT"},
        {"DraI",   "TTTAAA"},
        {"EcoRV",  "GATATC"},
        {"FseI",   "GGCCGGCC"},
        {"HpaI",   "GTTAAC"},
        {"MfeI",   "CAATTG"},
        {"NcoI",   "CCATGG"},
        {"NdeI",   "CATATG"},
        {"NruI",   "TCGCGA"},
        {"PacI",   "TTAATTAA"},
        {"PmeI",   "GTTTAAAC"},
        {"PvuII",  "CAGCTG"},
        {"SbfI",   "CCTGCAGG"},
        {"ScaI",   "AGTACT"},
        {"SnaBI",  "TACGTA"},
        {"SspI",   "AATATT"},
        {"StuI",   "AGGCCT"},
        {"SwaI",   "ATTTAAAT"},
        // 4-5 bp cutters (common in cloning workflows)
        {"AluI",   "AGCT"},
        {"HaeIII", "GGCC"},
        {"MseI",   "TTAA"},
        {"MspI",   "CCGG"},
        {"RsaI",   "GTAC"},
        {"TaqI",   "TCGA"},
    };
    return sites;
}

/// Compute reverse complement of a DNA string.
inline std::string reverse_complement(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (auto it = s.rbegin(); it != s.rend(); ++it) {
        switch (*it) {
            case 'A': out.push_back('T'); break;
            case 'T': out.push_back('A'); break;
            case 'C': out.push_back('G'); break;
            case 'G': out.push_back('C'); break;
            default:  out.push_back('?'); break;
        }
    }
    return out;
}

/// Check restriction sites (forward + reverse complement if not palindromic).
inline std::vector<MotifHit> find_restriction_sites(std::string_view seq) {
    std::vector<MotifHit> hits;

    for (const auto& [name, motif] : restriction_site_list()) {
        // Forward
        size_t pos = 0;
        while ((pos = seq.find(motif, pos)) != std::string::npos) {
            MotifHit h;
            h.category    = "restriction_site";
            h.motif       = name + ":" + motif;
            h.position    = pos;
            h.description = "Common restriction enzyme recognition site";
            hits.push_back(std::move(h));
            pos++;
        }
        // Reverse complement (only needed when motif is not self-complementary/palindromic)
        std::string rc = reverse_complement(motif);
        if (rc == motif) continue; // palindromic, already covered by forward search
        pos = 0;
        while ((pos = seq.find(rc, pos)) != std::string::npos) {
            MotifHit h;
            h.category    = "restriction_site";
            h.motif       = name + "(rc):" + rc;
            h.position    = pos;
            h.description = "Common restriction enzyme recognition site (reverse complement)";
            hits.push_back(std::move(h));
            pos++;
        }
    }
    return hits;
}

// ---------------- Functional motifs ----------------
/// Conservative built-in list of short functional-like motifs.
/// These are simple sequence patterns, NOT organism-specific annotations.
inline const std::vector<std::pair<std::string, std::string>>& functional_motif_list() {
    static const std::vector<std::pair<std::string, std::string>> motifs = {
        // --- Core promoter elements ---
        {"TATA_box",                "TATAAA"},
        {"GC_box",                  "GGGCGG"},
        {"CCAAT_box",               "CCAAT"},
        {"CRE",                     "TGACGTCA"},
        {"E_box_canonical",         "CACGTG"},

        // --- Transcription factor binding sites ---
        {"AP1_site",                "TGACTCA"},
        {"NFkB_site",               "GGGACTTTCC"},
        {"STAT3_site",              "TTCCGGGAA"},
        {"OCT4_site",               "ATGCAAAT"},
        {"E2F_site",                "TTTCCCGC"},
        {"Ets_binding_core",        "CCGGAAGT"},

        // --- mRNA processing signals ---
        {"polyA_signal",            "AATAAA"},
        {"alt_polyA_signal",        "ATTAAA"},
        {"polyA_variant_AATACA",    "AATACA"},
        {"polyA_variant_GATAAA",    "GATAAA"},
        {"polyA_variant_AATAGA",    "AATAGA"},
        {"polyA_variant_AATGAA",    "AATGAA"},
        {"ARE_mRNA_stability",      "ATTTAT"},
        {"splice_donor_like",       "GTAGT"},
        {"splice_acceptor_like",    "CAGG"},
        {"Kozak_like_core",         "ACCATGG"},

        // --- Prokaryotic signals ---
        {"Shine_Dalgarno",          "AGGAGG"},
        {"Pribnow_box",             "TATAAT"},
        {"minus35_sigma70",         "TTGACA"},

        // --- Very short motifs (conservatively included) ---
        {"CpG",                     "CG"},
    };
    return motifs;
}

// ---------------- Main evaluation function ----------------
/// Evaluate all enabled motif filters on a DNA sequence.
/// seq should be uppercase A/C/G/T, but the function handles case insensitively.
inline MotifFilterResult evaluate_motif_filters(std::string_view seq_orig,
                                                  const MotifFilterOptions& opts) {
    MotifFilterResult result;
    result.passes = true;

    if (opts.motif_mode == "off") return result; // no-op

    std::string seq = uppercase_dna(seq_orig);

    bool exclude = (opts.motif_mode == "exclude");

    // 1. Homopolymer filter
    if (opts.filter_homopolymers) {
        auto runs = find_homopolymer_runs(seq, opts.max_homopolymer);
        for (auto& h : runs) {
            result.hits.push_back(std::move(h));
            if (exclude) result.passes = false;
        }
    }

    // 2. Low-complexity filter
    if (opts.filter_low_complexity) {
        if (has_low_entropy(seq, opts.min_shannon_entropy)) {
            MotifHit h;
            h.category    = "low_complexity";
            h.motif       = "entropy";
            h.position    = 0;
            h.description = "Sequence entropy below threshold (" +
                            std::to_string(shannon_entropy(seq)) +
                            " < " + std::to_string(opts.min_shannon_entropy) + ")";
            result.hits.push_back(std::move(h));
            if (exclude) result.passes = false;
        }
    }

    // 3. Dinucleotide repeat filter
    if (opts.filter_dinucleotide_repeats) {
        auto reps = find_dinucleotide_repeats(seq);
        for (auto& h : reps) {
            result.hits.push_back(std::move(h));
            if (exclude) result.passes = false;
        }
    }

    // 3b. Trinucleotide repeat filter
    if (opts.filter_trinucleotide_repeats) {
        auto reps = find_trinucleotide_repeats(seq);
        for (auto& h : reps) {
            result.hits.push_back(std::move(h));
            if (exclude) result.passes = false;
        }
    }

    // 3c. Tetranucleotide repeat filter
    if (opts.filter_tetranucleotide_repeats) {
        auto reps = find_tetranucleotide_repeats(seq);
        for (auto& h : reps) {
            result.hits.push_back(std::move(h));
            if (exclude) result.passes = false;
        }
    }

    // 4. Restriction-site filter
    if (opts.filter_restriction_sites) {
        auto sites = find_restriction_sites(seq);
        for (auto& h : sites) {
            result.hits.push_back(std::move(h));
            if (exclude) result.passes = false;
        }
    }

    // 5. Functional-motif filter
    if (opts.filter_functional_motifs) {
        auto motifs = find_literal_motifs(seq, functional_motif_list(), "functional_motif");
        for (auto& h : motifs) {
            result.hits.push_back(std::move(h));
            if (exclude) result.passes = false;
        }
    }

    return result;
}

/// Convert motif hits to a compact JSON-like string for embedding in TSV output.
/// Format: [category:motif@pos|...] or empty if none.
inline std::string motif_hits_to_string(const std::vector<MotifHit>& hits) {
    if (hits.empty()) return "";
    std::string out;
    for (size_t i = 0; i < hits.size(); ++i) {
        if (i > 0) out += "|";
        out += hits[i].category + ":" + hits[i].motif + "@" + std::to_string(hits[i].position);
    }
    return out;
}

#endif // MOTIF_FILTERS_HPP
