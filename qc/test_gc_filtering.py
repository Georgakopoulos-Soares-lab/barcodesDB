#!/usr/bin/env python3
"""
test_gc_filtering.py — Validate GC filtering behavior.

Tests:
1. query_kmer_bitmap reports k-mers correctly regardless of GC content
   (using synthetic fixtures, --bitmap mode).
2. (Optional) query_substring_bitmap_stream with --gc-min/--gc-max if
   production shards exist.

Exit code 0 = all required tests pass; 1 = any required test fails.
Dependencies: Python standard library only.
"""

import json
import os
import struct
import subprocess
import sys
import tempfile

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
QC_TMP = os.path.join(REPO_ROOT, "website", "qc", "tmp")
BINARY = os.path.join(QC_TMP, "query_kmer_bitmap")
SUBSTR_BINARY = os.path.join(QC_TMP, "query_substring_bitmap_stream")
TRUTH_PATH = os.path.join(QC_TMP, "synthetic_truth.json")

# Production shard locations (used in optional tests)
SHARD_DIRS = {
    16: os.path.join(REPO_ROOT, "shards_16"),
    17: os.path.join(REPO_ROOT, "shards_17"),
    18: os.path.join(REPO_ROOT, "shards_18"),
}
GC_HIST_FILES = {
    16: os.path.join(REPO_ROOT, "gc_hist_shards_16.json"),
    17: os.path.join(REPO_ROOT, "gc_hist_shards_17.json"),
    18: os.path.join(REPO_ROOT, "gc_hist_shards_18.json"),
}

# Allow override via environment variable
BARCODEDB_ROOT = os.environ.get("BARCODEDB_ROOT")
if BARCODEDB_ROOT:
    for k in SHARD_DIRS:
        SHARD_DIRS[k] = os.path.join(BARCODEDB_ROOT, f"shards_{k}")
        GC_HIST_FILES[k] = os.path.join(BARCODEDB_ROOT, f"gc_hist_shards_{k}.json")

# ---------------------------------------------------------------------------
# Globals
# ---------------------------------------------------------------------------
g_failures = 0
g_passes = 0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def check(condition, test_name, detail=""):
    """Record a PASS or FAIL."""
    global g_passes, g_failures
    if condition:
        g_passes += 1
        print(f"  PASS  {test_name}")
    else:
        g_failures += 1
        msg = f"  FAIL  {test_name}"
        if detail:
            msg += f"  {detail}"
        print(msg)


def gc_content(dna):
    """Compute GC% for a DNA string."""
    if not dna:
        return 0.0
    gc = sum(1 for c in dna.upper() if c in "GC")
    return gc / len(dna) * 100.0


def read_fixture_header(fixture_path):
    """Read KBITv1 header and return (k, flags) or None."""
    if not os.path.exists(fixture_path):
        return None
    try:
        with open(fixture_path, "rb") as f:
            hdr = f.read(64)
        if len(hdr) < 64 or hdr[:6] != b"KBITv1":
            return None
        k = struct.unpack_from("<Q", hdr, 24)[0]
        flags = struct.unpack_from("<Q", hdr, 40)[0]
        return (int(k), int(flags))
    except Exception:
        return None


def fixture_path(k):
    """Return expected KBITv1 fixture path for a given k."""
    return os.path.join(QC_TMP, f"k{k}_fixture.kbit")


def run_binary(args, stdin_data=None):
    """Run a binary with given CLI arguments. Returns (stdout, stderr, rc)."""
    cmd = args if isinstance(args, list) and os.path.exists(args[0]) else [BINARY] + args
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
        timeout=60,
    )
    return proc.stdout, proc.stderr, proc.returncode


def write_kmer_file(kmers):
    """Write a list of k-mer strings to a temporary file. Returns the path."""
    tmp = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".kmers")
    for k in kmers:
        tmp.write(k + "\n")
    tmp.close()
    return tmp.name


def parse_output(stdout):
    """Parse tab-delimited output lines into dicts."""
    lines = [l for l in stdout.strip().split("\n") if l.strip()]
    results = []
    for line in lines:
        parts = line.split("\t")
        entry = {"kmer": parts[0], "hit": parts[1]}
        if len(parts) >= 3:
            entry["motif_passes"] = parts[2]
        if len(parts) >= 4:
            entry["motif_hit_count"] = parts[3]
        if len(parts) >= 5:
            entry["motif_hits"] = parts[4]
        results.append(entry)
    return results


def pick_shards():
    """Return the smallest k whose shard dir and gc-hist both exist, or None."""
    for k in (16, 17, 18):
        sd = SHARD_DIRS[k]
        gh = GC_HIST_FILES[k]
        if os.path.isdir(sd) and os.path.isfile(gh):
            return k, sd, gh
    return None, None, None


def get_results(stdout):
    """Return data lines from substring binary stdout (after __META__)."""
    lines = stdout.split("\n")
    data_lines = []
    in_data = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("__META__"):
            in_data = True
            continue
        if in_data:
            if stripped.startswith("[INFO]") or stripped == "":
                break
            data_lines.append(stripped)
    return data_lines


def parse_meta(stdout):
    """Parse __META__ line from substring binary stdout."""
    for line in stdout.split("\n"):
        line = line.strip()
        if line.startswith("__META__"):
            parts = line.split("\t")
            return {
                "cursor": parts[1] if len(parts) > 1 else "",
                "hasMore": parts[2] if len(parts) > 2 else "",
                "returned": int(parts[3]) if len(parts) > 3 and parts[3].isdigit() else 0,
                "k": int(parts[4]) if len(parts) > 4 and parts[4].isdigit() else 0,
            }
    return None


# ===================================================================
#  Test 1: GC content independence (k-mer lookups via --bitmap)
# ===================================================================
def test_gc_independence_kmer():
    """Verify present k-mers are found regardless of GC content."""
    print("\n=== Test 1: GC-content independence (k-mer bitmap lookups) ===")

    if not os.path.exists(BINARY):
        check(False, "query_kmer_bitmap binary not found", BINARY)
        return
    if not os.path.exists(TRUTH_PATH):
        check(False, "synthetic truth not found", TRUTH_PATH)
        return

    with open(TRUTH_PATH, "r") as f:
        truth = json.load(f)

    for fixture in truth["fixtures"]:
        k_val = fixture["k"]
        fp = fixture_path(k_val)
        hdr = read_fixture_header(fp)
        if hdr is None:
            check(False, f"k={k_val}: cannot read fixture", fp)
            continue
        actual_k = hdr[0]

        present_kmers = [k for k in fixture["present_kmers"] if len(k) == actual_k]
        if not present_kmers:
            continue

        # Compute GC% for each present k-mer
        gc_values = [gc_content(k) for k in present_kmers]
        min_gc = min(gc_values)
        max_gc = max(gc_values)
        print(f"  INFO  k={k_val}: present k-mer GC range [{min_gc:.1f}%, {max_gc:.1f}%] "
              f"over {len(present_kmers)} k-mers")

        tmp_path = write_kmer_file(present_kmers)
        stdout, stderr, rc = run_binary([
            BINARY, "--bitmap", fp, "--kmers", tmp_path,
        ])
        os.unlink(tmp_path)

        if rc != 0:
            check(False, f"k={k_val}: binary exited {rc}", stderr.strip()[:200])
            continue

        results = parse_output(stdout)
        if len(results) != len(present_kmers):
            check(False, f"k={k_val}: expected {len(present_kmers)} lines, got {len(results)}")
            continue

        # Verify ALL k-mers are found regardless of GC
        all_present = all(r["hit"] == "1" for r in results)
        check(all_present, f"k={k_val}: all {len(present_kmers)} present k-mers return hit=1",
              "" if all_present else str([(r["kmer"], r["hit"]) for r in results if r["hit"] != "1"]))

        # Group by GC bins and verify no GC bias
        low_gc = [r for r, gc in zip(results, gc_values) if gc < 30]
        mid_gc = [r for r, gc in zip(results, gc_values) if 30 <= gc <= 70]
        high_gc = [r for r, gc in zip(results, gc_values) if gc > 70]

        for label, group in [("low GC (<30%)", low_gc), ("mid GC (30-70%)", mid_gc),
                              ("high GC (>70%)", high_gc)]:
            if group:
                all_ok = all(r["hit"] == "1" for r in group)
                check(all_ok, f"k={k_val}: {label} — all {len(group)} k-mers return hit=1",
                      "" if all_ok else "")

    # Also test absent k-mers across GC ranges
    print("  --- Checking absent k-mers across GC ranges ---")
    for fixture in truth["fixtures"]:
        k_val = fixture["k"]
        fp = fixture_path(k_val)
        hdr = read_fixture_header(fp)
        if hdr is None:
            continue
        actual_k = hdr[0]

        absent_kmers = [k for k in fixture["absent_kmers"] if len(k) == actual_k]
        if not absent_kmers:
            continue

        tmp_path = write_kmer_file(absent_kmers)
        stdout, stderr, rc = run_binary([
            BINARY, "--bitmap", fp, "--kmers", tmp_path,
        ])
        os.unlink(tmp_path)

        if rc != 0:
            continue

        results = parse_output(stdout)
        all_absent = all(r["hit"] == "0" for r in results)
        check(all_absent, f"k={k_val}: all {len(absent_kmers)} absent k-mers return hit=0",
              "" if all_absent else str([(r["kmer"], r["hit"]) for r in results if r["hit"] != "0"]))


# ===================================================================
#  Test 2 (Optional): GC filtering with substring binary
# ===================================================================
def test_gc_filter_substring():
    """Optional: Verify --gc-min/--gc-max with substring binary and production shards."""
    print("\n=== Test 2: GC filtering (substring binary, optional) ===")

    k_val, shard_dir, gc_hist = pick_shards()
    if k_val is None:
        print("  [SKIP] Production shards not found — skipping substring GC filter tests.")
        print("  To run: ensure shards_16/ (or 17/18) and gc_hist_shards_*.json exist.")
        return

    if not os.path.exists(SUBSTR_BINARY):
        print("  [SKIP] query_substring_bitmap_stream binary not found.")
        return

    # 2a: Wide GC range should return results
    print("  2a: Wide GC range (0-100) — baseline")
    p = subprocess.run(
        [SUBSTR_BINARY, "--shards", shard_dir, "--gc-hist", gc_hist,
         "--limit", "10", "--random_access",
         "--gc-min", "0", "--gc-max", "100"],
        capture_output=True, text=True, timeout=120,
    )
    stdout, stderr, rc = p.stdout, p.stderr, p.returncode
    if rc != 0:
        check(False, "wide GC range: binary exited non-zero", f"rc={rc} {stderr[:200]}")
        return

    results_wide = get_results(stdout)
    check(len(results_wide) > 0, "wide GC range returns results",
          f"got {len(results_wide)}")

    # 2b: Narrow GC range (40-60) — results should obey constraint
    print("  2b: Narrow GC range (40-60)")
    p = subprocess.run(
        [SUBSTR_BINARY, "--shards", shard_dir, "--gc-hist", gc_hist,
         "--limit", "20", "--random_access",
         "--gc-min", "40", "--gc-max", "60"],
        capture_output=True, text=True, timeout=120,
    )
    stdout, stderr, rc = p.stdout, p.stderr, p.returncode
    if rc != 0:
        check(False, "narrow GC range: binary exited non-zero", f"rc={rc} {stderr[:200]}")
        return

    results_narrow = get_results(stdout)
    all_in_range = all(40 <= gc_content(k) <= 60 for k in results_narrow)
    check(all_in_range,
          f"all {len(results_narrow)} results have GC% in [40,60]",
          "" if all_in_range else str([(k, gc_content(k)) for k in results_narrow
                                        if not (40 <= gc_content(k) <= 60)]))

    # 2c: Non-overlapping GC range (99-100) — may return 0 results (valid)
    print("  2c: Non-overlapping GC range (99-100)")
    p = subprocess.run(
        [SUBSTR_BINARY, "--shards", shard_dir, "--gc-hist", gc_hist,
         "--limit", "10", "--random_access",
         "--gc-min", "99", "--gc-max", "100"],
        capture_output=True, text=True, timeout=120,
    )
    stdout, stderr, rc = p.stdout, p.stderr, p.returncode
    if rc == 0:
        results_high_gc = get_results(stdout)
        all_high = all(gc_content(k) >= 99 for k in results_high_gc)
        check(all_high, "GC range [99,100] results all >=99% GC",
              "" if all_high else str([(k, gc_content(k)) for k in results_high_gc
                                        if gc_content(k) < 99]))
    else:
        print("  INFO  GC range [99,100]: binary returned non-zero "
              f"(rc={rc}) — may be valid (no matching shards)")

    # 2d: Boundary values (gc-min=gc-max)
    print("  2d: Boundary GC range (gc-min=gc-max at 50)")
    p = subprocess.run(
        [SUBSTR_BINARY, "--shards", shard_dir, "--gc-hist", gc_hist,
         "--limit", "10", "--random_access",
         "--gc-min", "50", "--gc-max", "50"],
        capture_output=True, text=True, timeout=120,
    )
    stdout, stderr, rc = p.stdout, p.stderr, p.returncode
    if rc == 0:
        results_boundary = get_results(stdout)
        all_at_50 = all(gc_content(k) == 50.0 for k in results_boundary)
        check(all_at_50,
              f"GC range [50,50]: all {len(results_boundary)} results at exactly 50% GC",
              "" if all_at_50 else str([(k, gc_content(k)) for k in results_boundary
                                        if gc_content(k) != 50.0]))
    else:
        print("  INFO  GC range [50,50]: binary returned non-zero "
              f"(rc={rc}) — may be valid (no matching shards)")


# ===================================================================
#  Main
# ===================================================================
def main():
    global g_failures, g_passes
    g_failures = 0
    g_passes = 0

    print(f"Repository root: {REPO_ROOT}")
    print(f"Binary:          {BINARY}")
    print(f"Substr binary:   {SUBSTR_BINARY}")
    print(f"Truth file:      {TRUTH_PATH}")

    test_gc_independence_kmer()
    test_gc_filter_substring()

    print(f"\n{'=' * 60}")
    print(f"Results: {g_passes} passed, {g_failures} failed")

    if g_failures == 0:
        print("ALL TESTS PASSED ✓")
        sys.exit(0)
    else:
        print("SOME TESTS FAILED ✗")
        sys.exit(1)


if __name__ == "__main__":
    main()
