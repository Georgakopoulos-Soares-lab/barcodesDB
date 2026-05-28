#!/usr/bin/env bash
# test_motif_filters.sh — Quick smoke tests for motif filter logic
# Run from the website/programs directory
set -euo pipefail

echo "=== Motif Filter Unit Tests ==="
echo ""

# Build a minimal test binary that exercises motif_filters.hpp directly
cat > /tmp/test_motif.cpp << 'CPPEOF'
#include "motif_filters.hpp"
#include <iostream>
#include <cassert>
#include <string>

static int tests = 0;
static int passed = 0;
static int failed = 0;

#define TEST(name) do { tests++; std::cout << "  TEST " << tests << ": " << name << " ... "; } while(0)
#define PASS() do { passed++; std::cout << "PASS\n"; } while(0)
#define FAIL(msg) do { failed++; std::cout << "FAIL: " << msg << "\n"; } while(0)

int main() {

  // --- Homopolymer ---
  TEST("AAAAACGT fails with max_homopolymer=4");
  {
    auto hits = find_homopolymer_runs("AAAAACGT", 4);
    if (!hits.empty()) PASS(); else FAIL("expected homopolymer hit AAAAA");
  }

  TEST("AAAACGT passes with max_homopolymer=4");
  {
    auto hits = find_homopolymer_runs("AAAACGT", 4);
    if (hits.empty()) PASS(); else FAIL("expected no hits for length-4 run");
  }

  TEST("CCCCC fails with max_homopolymer=4");
  {
    auto hits = find_homopolymer_runs("CCCCC", 4);
    if (!hits.empty()) PASS(); else FAIL("expected homopolymer hit CCCCC");
  }

  TEST("ACGT passes homopolymer check");
  {
    auto hits = find_homopolymer_runs("ACGT", 4);
    if (hits.empty()) PASS(); else FAIL("expected no hits");
  }

  // --- Low complexity ---
  TEST("low complexity: AAAAAAAACCCCCCCC has lower entropy than ACGTACGTACGT");
  {
    double h1 = shannon_entropy("AAAAAAAAACCCCCCCC");
    double h2 = shannon_entropy("ACGTACGTACGT");
    if (h1 < h2) PASS(); else FAIL("expected lower entropy for biased sequence");
  }

  TEST("entropy threshold: all-A fails with min_entropy=1.5");
  {
    bool low = has_low_entropy("AAAAAAAAAAAAAAAAAA", 1.5);
    if (low) PASS(); else FAIL("all-A should be below entropy threshold");
  }

  TEST("entropy threshold: balanced passes with min_entropy=1.5");
  {
    bool low = has_low_entropy("ACGTACGTACGTACGTAC", 1.5);
    if (!low) PASS(); else FAIL("balanced sequence should pass entropy threshold");
  }

  // --- Dinucleotide repeats ---
  TEST("ATATAT triggers dinucleotide repeat");
  {
    auto hits = find_dinucleotide_repeats("ATATAT");
    if (!hits.empty()) PASS(); else FAIL("expected dinucleotide repeat hit for ATATAT");
  }

  TEST("ACGTAC does NOT trigger dinucleotide repeat");
  {
    auto hits = find_dinucleotide_repeats("ACGTAC");
    if (hits.empty()) PASS(); else FAIL("expected no dinucleotide repeat for ACGTAC");
  }

  TEST("ACACACAC triggers dinucleotide repeat");
  {
    auto hits = find_dinucleotide_repeats("ACACACAC");
    if (!hits.empty()) PASS(); else FAIL("expected dinucleotide repeat hit for ACACACAC");
  }

  // --- Restriction sites ---
  TEST("EcoRI site GAATTC detected");
  {
    auto hits = find_restriction_sites("ACGTGAATTCACGT");
    bool found = false;
    for (auto& h : hits) if (h.motif.find("EcoRI") != std::string::npos) found = true;
    if (found) PASS(); else FAIL("expected EcoRI hit");
  }

  TEST("No false positive restriction site");
  {
    auto hits = find_restriction_sites("AACTGAATACTGAATC");
    if (hits.empty()) PASS(); else FAIL("expected no restriction hits");
  }

  // --- Functional motifs ---
  TEST("AATAAA triggers polyA signal");
  {
    auto hits = find_literal_motifs("CAATAAAG", functional_motif_list(), "functional_motif");
    bool found = false;
    for (auto& h : hits) if (h.motif.find("polyA_signal") != std::string::npos) found = true;
    if (found) PASS(); else FAIL("expected polyA_signal hit");
  }

  TEST("TATAAA triggers TATA box-like motif");
  {
    auto hits = find_literal_motifs("ATATAAAC", functional_motif_list(), "functional_motif");
    bool found = false;
    for (auto& h : hits) if (h.motif.find("TATA_box") != std::string::npos) found = true;
    if (found) PASS(); else FAIL("expected TATA_box hit");
  }

  TEST("Clean sequence has no functional motif hits");
  {
    auto hits = find_literal_motifs("AAAAAAAAAAAAAAAAAA", functional_motif_list(), "functional_motif");
    if (hits.empty()) PASS(); else FAIL("expected no functional motif hits for all-A");
  }

  // --- evaluate_motif_filters ---
  TEST("evaluate_motif_filters: mode=off returns empty");
  {
    MotifFilterOptions opts;
    opts.motif_mode = "off";
    opts.filter_homopolymers = true;
    auto result = evaluate_motif_filters("AAAAAA", opts);
    if (result.passes && result.hits.empty()) PASS(); else FAIL("off mode should skip all checks");
  }

  TEST("evaluate_motif_filters: mode=flag reports but passes");
  {
    MotifFilterOptions opts;
    opts.motif_mode = "flag";
    opts.filter_homopolymers = true;
    opts.max_homopolymer = 4;
    auto result = evaluate_motif_filters("AAAAAA", opts);
    if (result.passes && !result.hits.empty()) PASS(); else FAIL("flag mode should report but not fail");
  }

  TEST("evaluate_motif_filters: mode=exclude fails homopolymer");
  {
    MotifFilterOptions opts;
    opts.motif_mode = "exclude";
    opts.filter_homopolymers = true;
    opts.max_homopolymer = 4;
    auto result = evaluate_motif_filters("AAAAAA", opts);
    if (!result.passes && !result.hits.empty()) PASS(); else FAIL("exclude mode should fail on homopolymer");
  }

  // --- motif_hits_to_string ---
  TEST("motif_hits_to_string formats correctly");
  {
    std::vector<MotifHit> hits;
    MotifHit h;
    h.category = "homopolymer";
    h.motif = "AAAAA";
    h.position = 0;
    hits.push_back(h);
    std::string s = motif_hits_to_string(hits);
    if (s == "homopolymer:AAAAA@0") PASS(); else FAIL("got: " + s);
  }

  std::cout << "\n=== Results: " << passed << "/" << tests << " passed";
  if (failed > 0) std::cout << ", " << failed << " FAILED";
  std::cout << " ===\n";
  return failed > 0 ? 1 : 0;
}
CPPEOF

echo "Compiling test binary..."
cd /home1/10899/kimopro/WORK/barcodes/website/programs
g++ -std=c++17 -I. /tmp/test_motif.cpp -o /tmp/test_motif 2>&1

echo ""
echo "Running tests..."
/tmp/test_motif

echo ""
echo "=== Integration: query_kmer_bitmap with --motif-mode=flag ==="
echo "Testing with a known k-mer containing homopolymer..."
if [ -x "./query_kmer_bitmap_test" ]; then
  echo "AAAAAACCCCCCCCCC" > /tmp/test_kmers.txt
  ./query_kmer_bitmap_test --shards /home1/10899/kimopro/WORK/barcodes/shards_16 --k 16 --kmers /tmp/test_kmers.txt --motif-mode flag --filter-homopolymers --max-homopolymer 4 2>&1 | head -5
  rm -f /tmp/test_kmers.txt
else
  echo "  (integration binary 'query_kmer_bitmap_test' not found; skipping integration step)"
fi

echo ""
echo "=== Done ==="
rm -f /tmp/test_motif.cpp /tmp/test_motif /tmp/test_kmers.txt
rm -f ./query_kmer_bitmap_test ./query_substring_bitmap_stream_test 2>/dev/null; true
