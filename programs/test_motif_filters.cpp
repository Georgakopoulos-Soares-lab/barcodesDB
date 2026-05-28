// test_motif_filters.cpp — Quick test program for motif_filters.hpp
//
// Compile:
//   g++ -std=c++17 -O2 -o test_motif_filters test_motif_filters.cpp
//
// Run:
//   ./test_motif_filters

#include "motif_filters.hpp"
#include <cassert>
#include <iostream>
#include <string>

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name, expr) do { \
    tests_run++; \
    if (!(expr)) { std::cerr << "FAIL: " << name << " (" << __LINE__ << ")\n"; } \
    else { tests_passed++; std::cout << "  PASS: " << name << "\n"; } \
} while(0)

int main() {
    std::cout << "=== motif_filters.hpp tests ===\n\n";

    // --- Homopolymer ---
    {
        std::cout << "-- Homopolymer filter --\n";
        MotifFilterOptions opts;
        opts.motif_mode = "exclude";
        opts.filter_homopolymers = true;
        opts.max_homopolymer = 4;

        // AAAAACGT fails with max_homopolymer=4 (run of 5 A's)
        auto res1 = evaluate_motif_filters("AAAAACGT", opts);
        TEST("AAAAACGT fails homopolymer", res1.passes == false && res1.hits.size() >= 1);

        // AAAACGT passes with max_homopolymer=4
        auto res2 = evaluate_motif_filters("AAAACGT", opts);
        TEST("AAAACGT passes homopolymer", res2.passes == true && res2.hits.empty());

        // AAAAA homopolymer on its own
        auto res3 = evaluate_motif_filters("AAAAA", opts);
        TEST("AAAAA fails homopolymer", res3.passes == false);

        // All different bases
        auto res4 = evaluate_motif_filters("ACGTACGTACGT", opts);
        TEST("ACGTACGTACGT passes homopolymer", res4.passes == true);
    }

    // --- Restriction sites ---
    {
        std::cout << "-- Restriction site filter --\n";
        MotifFilterOptions opts;
        opts.motif_mode = "exclude";
        opts.filter_restriction_sites = true;

        auto res1 = evaluate_motif_filters("GAATTC", opts);
        TEST("GAATTC triggers EcoRI", res1.passes == false && res1.hits.size() >= 1);

        bool has_ecori = false;
        for (auto& h : res1.hits) if (h.motif.find("EcoRI") != std::string::npos) has_ecori = true;
        TEST("EcoRI hit present", has_ecori);

        auto res2 = evaluate_motif_filters("GGATCC", opts);
        TEST("GGATCC triggers BamHI", res2.passes == false);

        // Sequence without listed motifs (harder now with 50+ enzymes — use verified-clean sequence)
        auto res3 = evaluate_motif_filters("AATACTGAATACTGAATC", opts);
        TEST("Clean sequence passes restriction filter", res3.passes == true);
    }

    // --- Low complexity / Shannon entropy ---
    {
        std::cout << "-- Low-complexity filter --\n";
        MotifFilterOptions opts;
        opts.motif_mode = "exclude";
        opts.filter_low_complexity = true;
        opts.min_shannon_entropy = 1.5;

        // AAAAAAAAAAAAAAAA has very low entropy
        double h1 = shannon_entropy("AAAAAAAAAAAAAAAA");
        TEST("AAAA... has low entropy", h1 < 0.5);

        // ACGTACGTACGT has balanced entropy (~2.0)
        double h2 = shannon_entropy("ACGTACGTACGT");
        TEST("ACGT... has high entropy", h2 > 1.5);

        auto res1 = evaluate_motif_filters("AAAAAAAAAAAAAAAA", opts);
        TEST("AAAA... fails low-complexity filter", res1.passes == false);

        auto res2 = evaluate_motif_filters("ACGTACGTACGTACGT", opts);
        TEST("ACGT... passes low-complexity filter", res2.passes == true);
    }

    // --- Dinucleotide repeats ---
    {
        std::cout << "-- Dinucleotide repeat filter --\n";
        MotifFilterOptions opts;
        opts.motif_mode = "exclude";
        opts.filter_dinucleotide_repeats = true;

        auto res1 = evaluate_motif_filters("ATATATATCG", opts);
        TEST("ATATAT triggers dinucleotide repeat", res1.passes == false);

        auto res2 = evaluate_motif_filters("ACGTACGTAC", opts);
        TEST("ACGTAC does not trigger dinucleotide repeat", res2.passes == true);
    }

    // --- Functional motifs ---
    {
        std::cout << "-- Functional motif filter --\n";
        MotifFilterOptions opts;
        opts.motif_mode = "flag";
        opts.filter_functional_motifs = true;

        auto res1 = evaluate_motif_filters("AATAAA", opts);
        TEST("AATAAA triggers polyA signal", res1.passes == true && res1.hits.size() >= 1);
        bool has_polya = false;
        for (auto& h : res1.hits) if (h.motif.find("polyA_signal") != std::string::npos) has_polya = true;
        TEST("polyA_signal hit present in flag mode", has_polya);

        auto res2 = evaluate_motif_filters("TATAAA", opts);
        bool has_tata = false;
        for (auto& h : res2.hits) if (h.motif.find("TATA_box") != std::string::npos) has_tata = true;
        TEST("TATAAA triggers TATA-like motif", has_tata);
    }

    // --- Mode: off ---
    {
        std::cout << "-- Mode: off --\n";
        MotifFilterOptions opts;
        opts.motif_mode = "off";
        opts.filter_homopolymers = true;
        opts.max_homopolymer = 4;

        auto res = evaluate_motif_filters("AAAAA", opts);
        TEST("off mode: passes and no hits", res.passes == true && res.hits.empty());
    }

    // --- Mode: flag vs exclude ---
    {
        std::cout << "-- Mode: flag vs exclude --\n";
        MotifFilterOptions opts_flag;
        opts_flag.motif_mode = "flag";
        opts_flag.filter_homopolymers = true;
        opts_flag.max_homopolymer = 4;
        auto res_flag = evaluate_motif_filters("AAAAA", opts_flag);
        TEST("flag mode: passes even with homopolymer fail", res_flag.passes == true && !res_flag.hits.empty());

        MotifFilterOptions opts_excl;
        opts_excl.motif_mode = "exclude";
        opts_excl.filter_homopolymers = true;
        opts_excl.max_homopolymer = 4;
        auto res_excl = evaluate_motif_filters("AAAAA", opts_excl);
        TEST("exclude mode: fails with homopolymer fail", res_excl.passes == false && !res_excl.hits.empty());
    }

    // --- motif_hits_to_string ---
    {
        std::cout << "-- motif_hits_to_string --\n";
        std::vector<MotifHit> hits;
        MotifHit h;
        h.category = "homopolymer";
        h.motif = "AAAAA";
        h.position = 0;
        hits.push_back(h);
        std::string s = motif_hits_to_string(hits);
        TEST("motif_hits_to_string non-empty", !s.empty() && s.find("AAAAA") != std::string::npos);
    }

    // --- uppercase_dna ---
    {
        std::cout << "-- uppercase_dna --\n";
        std::string s = uppercase_dna("acgtACGT");
        TEST("uppercase_dna", s == "ACGTACGT");
    }

    std::cout << "\n=== Results: " << tests_passed << "/" << tests_run << " passed ===\n";
    return (tests_passed == tests_run) ? 0 : 1;
}
