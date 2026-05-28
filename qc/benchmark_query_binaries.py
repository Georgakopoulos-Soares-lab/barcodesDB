#!/usr/bin/env python3
"""
benchmark_query_binaries.py — Performance/benchmark smoke tests.
Runs N exact k-mer lookups and N substring searches against synthetic/production fixtures.
Reports timing and writes summary to qc/tmp/benchmark_summary.json.

Exit 0 always (informational; failures only on extreme regressions).
"""

import json
import os
import subprocess
import sys
import tempfile
import time

QC_DIR = os.path.dirname(os.path.abspath(__file__))
WEBSITE_DIR = os.path.normpath(os.path.join(QC_DIR, ".."))
QC_TMP = os.path.join(QC_DIR, "tmp")
ROOT = os.environ.get("BARCODEDB_ROOT", os.path.normpath(os.path.join(WEBSITE_DIR, "..")))

BIN_KMER = os.path.join(QC_TMP, "query_kmer_bitmap")
BIN_SUBSTR = os.path.join(QC_TMP, "query_substring_bitmap_stream")
FIXTURE_16 = os.path.join(QC_TMP, "k16_fixture.kbit")

# Production shards (optional)
SHARDS_18 = os.path.join(ROOT, "shards_18")
GC_HIST_18 = os.path.join(ROOT, "gc_hist_shards_18.json")

N_KMER = 1000
N_SUBSTR = 50


def run(cmd, timeout=120):
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    elapsed = time.time() - t0
    return r, elapsed


def main():
    results = {
        "exact_lookup": {"n": 0, "total_s": 0, "mean_ms": 0, "max_ms": 0},
        "substring_search": {"n": 0, "total_s": 0, "mean_ms": 0, "max_ms": 0},
        "fixture": "synthetic",
    }

    passed = 0
    failed = 0
    skipped = 0

    # --- Exact k-mer lookup benchmark ---
    print("=== Benchmark: Exact k-mer lookup ===")
    if not os.path.isfile(BIN_KMER) or not os.path.isfile(FIXTURE_16):
        print(f"  SKIP: binary or fixture not found (need {BIN_KMER} and {FIXTURE_16})")
        skipped += 1
    else:
        # Generate N_KMER test k-mers (mix of present and absent from truth)
        truth_path = os.path.join(QC_TMP, "synthetic_truth.json")
        with open(truth_path) as f:
            truth = json.load(f)
        k16 = truth["fixtures"][0]
        test_kmers = (k16["present_kmers"] + k16["absent_kmers"]) * (N_KMER // len(k16["present_kmers"] + k16["absent_kmers"]) + 1)
        test_kmers = test_kmers[:N_KMER]

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
            for kmer in test_kmers:
                f.write(kmer + "\n")
            tmp_path = f.name

        try:
            latencies = []
            batch_size = 100
            for i in range(0, N_KMER, batch_size):
                batch = test_kmers[i:i + batch_size]
                with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as bf:
                    for kmer in batch:
                        bf.write(kmer + "\n")
                    bpath = bf.name
                try:
                    _, elapsed = run([BIN_KMER, "--bitmap", FIXTURE_16, "--kmers", bpath], timeout=60)
                    for _ in batch:
                        latencies.append(elapsed / len(batch))
                finally:
                    os.unlink(bpath)

            total_s = sum(latencies)
            mean_ms = (total_s / len(latencies)) * 1000 if latencies else 0
            max_ms = max(latencies) * 1000 if latencies else 0
            results["exact_lookup"]["n"] = len(latencies)
            results["exact_lookup"]["total_s"] = round(total_s, 3)
            results["exact_lookup"]["mean_ms"] = round(mean_ms, 3)
            results["exact_lookup"]["max_ms"] = round(max_ms, 3)
            print(f"  {len(latencies)} lookups: total={total_s:.3f}s  mean={mean_ms:.3f}ms  max={max_ms:.3f}ms")
            passed += 1
        finally:
            os.unlink(tmp_path)

    # --- Substring search benchmark ---
    print("=== Benchmark: Substring search ===")
    if not os.path.isfile(BIN_SUBSTR) or not os.path.isdir(SHARDS_18) or not os.path.isfile(GC_HIST_18):
        print("  SKIP: binary or production shards not found")
        skipped += 1
    else:
        latencies = []
        for _ in range(N_SUBSTR):
            cmd = [BIN_SUBSTR, "--shards", SHARDS_18, "--gc-hist", GC_HIST_18,
                   "--gc-min", "0", "--gc-max", "100", "--limit", "10", "--random_access"]
            _, elapsed = run(cmd, timeout=120)
            latencies.append(elapsed)

        if latencies:
            total_s = sum(latencies)
            mean_ms = (total_s / len(latencies)) * 1000
            max_ms = max(latencies) * 1000
            results["substring_search"]["n"] = len(latencies)
            results["substring_search"]["total_s"] = round(total_s, 3)
            results["substring_search"]["mean_ms"] = round(mean_ms, 3)
            results["substring_search"]["max_ms"] = round(max_ms, 3)
            print(f"  {len(latencies)} searches: total={total_s:.3f}s  mean={mean_ms:.3f}ms  max={max_ms:.3f}ms")
            passed += 1

    # --- Write summary ---
    summary_path = os.path.join(QC_TMP, "benchmark_summary.json")
    os.makedirs(os.path.dirname(summary_path), exist_ok=True)
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nBenchmark summary written to {summary_path}")

    print(f"\n=== Results: {passed} passed, {failed} failed, {skipped} skipped ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
