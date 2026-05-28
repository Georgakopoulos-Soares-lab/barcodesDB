#!/usr/bin/env python3
"""
test_query_substring_bitmap_stream.py — Validate query_substring_bitmap_stream binary.

This binary requires sharded data (--shards <dir> with index.json + .kbit files)
and a GC histogram JSON (--gc-hist).  The main tests are OPTIONAL — they run only
when production shard directories exist on disk.  One REQUIRED test always runs:
verifying that missing/invalid arguments produce a non-zero exit code.

Run from anywhere (paths are relative to the script location):
    python3 website/qc/test_query_substring_bitmap_stream.py

Exit code 0 = all tests pass (or optional tests are skipped); 1 = any failure.
Dependencies: Python standard library only.
"""

import os
import subprocess
import sys

# ---------------------------------------------------------------------------
# Paths (relative to this script's location → repository root)
# ---------------------------------------------------------------------------
QC_DIR = os.path.abspath(os.path.dirname(__file__))
REPO_ROOT = os.path.abspath(os.path.join(QC_DIR, "../.."))
BINARY = os.path.join(QC_DIR, "tmp", "query_substring_bitmap_stream")

# Production shard locations (relative to repo root)
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
    """Record a PASS or FAIL for the given test."""
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


def pick_shards():
    """Return the smallest k whose shard dir and gc-hist both exist, or None."""
    for k in (16, 17, 18):
        sd = SHARD_DIRS[k]
        gh = GC_HIST_FILES[k]
        if os.path.isdir(sd) and os.path.isfile(gh):
            return k, sd, gh
    return None, None, None


def run_binary(args, timeout=120):
    """
    Run query_substring_bitmap_stream with the given CLI arguments.

    Returns (stdout, stderr, returncode).
    """
    cmd = [BINARY] + args
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    return proc.stdout.strip(), proc.stderr.strip(), proc.returncode


def parse_meta(stdout):
    """Parse the __META__ line from stdout.

    Returns a dict with keys: cursor, count1, count_returned, k, motif_info
    or None if no __META__ line found.
    """
    for line in stdout.split("\n"):
        line = line.strip()
        if line.startswith("__META__"):
            parts = line.split("\t")
            # __META__  <cursor>  <count1>  <returned_count>  <k>  [motif...]
            meta = {
                "cursor": parts[1] if len(parts) > 1 else "",
                "count1": parts[2] if len(parts) > 2 else "",
                "returned": int(parts[3]) if len(parts) > 3 and parts[3].isdigit() else 0,
                "k": int(parts[4]) if len(parts) > 4 and parts[4].isdigit() else 0,
                "motif_info": "\t".join(parts[5:]) if len(parts) > 5 else "",
            }
            return meta
    return None


def get_results(stdout):
    """Return the data lines from stdout (lines before any [INFO] prefix, after __META__)."""
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


def gc_content(dna):
    """Return GC percentage for a DNA string."""
    if not dna:
        return 0.0
    gc = sum(1 for c in dna.upper() if c in "GC")
    return gc / len(dna) * 100.0


# ===================================================================
#  Test 0: REQUIRED — --help / missing args produce non-zero exit
# ===================================================================
def test_required_usage_errors():
    """Verify missing/invalid args produce exit code != 0 and a usage message."""
    print("\n=== Test 0: REQUIRED — usage errors ===")

    # 0a: No arguments at all
    stdout, stderr, rc = run_binary([])
    ok_noargs = rc != 0 and ("Usage:" in stdout or "Usage:" in stderr)
    check(ok_noargs, "0a: no args → non-zero exit + usage message",
          f"rc={rc}" if not ok_noargs else "")

    # 0b: Missing required --shards (but providing --gc-hist)
    stdout, stderr, rc = run_binary(["--gc-hist", "dummy.json"])
    ok_noshards = rc != 0 and ("Usage:" in stdout or "Usage:" in stderr or "--shards" in stdout or "--shards" in stderr)
    check(ok_noshards, "0b: missing --shards → non-zero exit",
          f"rc={rc}" if not ok_noshards else "")

    # 0c: Missing --gc-hist (but providing --shards)
    stdout, stderr, rc = run_binary(["--shards", "dummy_dir"])
    ok_noghist = rc != 0 and ("Usage:" in stdout or "Usage:" in stderr or "--gc-hist" in stdout or "--gc-hist" in stderr)
    check(ok_noghist, "0c: missing --gc-hist → non-zero exit",
          f"rc={rc}" if not ok_noghist else "")

    # 0d: Invalid flag
    stdout, stderr, rc = run_binary(["--nonexistent-flag"])
    ok_badflag = rc != 0 and ("Usage:" in stdout or "Usage:" in stderr or "Unknown arg" in stdout or "Unknown arg" in stderr)
    check(ok_badflag, "0d: invalid flag → non-zero exit",
          f"rc={rc}" if not ok_badflag else "")


# ===================================================================
#  OPTIONAL shard-based tests
# ===================================================================
def test_basic_search(k_val, shard_dir, gc_hist):
    """Basic search with --random_access --limit 10 — verify __META__ and results."""
    print(f"\n=== Test 1: BASIC SEARCH (k={k_val}) ===")
    stdout, stderr, rc = run_binary([
        "--shards", shard_dir,
        "--gc-hist", gc_hist,
        "--limit", "10",
        "--random_access",
    ])
    if rc != 0:
        check(False, "binary exited non-zero", f"rc={rc} {stderr[:200]}")
        return

    meta = parse_meta(stdout)
    check(meta is not None, "__META__ line present")
    if meta:
        check(meta["returned"] > 0, f"returned count > 0 (got {meta['returned']})")
        check(meta["cursor"] != "", "cursor is non-empty")
        check(meta["k"] == k_val, f"k matches ({meta['k']} == {k_val})")

    results = get_results(stdout)
    check(len(results) > 0, f"got {len(results)} result k-mers")
    check(len(results) <= 10, f"at most 10 results (got {len(results)})")


def test_limit(k_val, shard_dir, gc_hist):
    """Verify --limit N returns at most N results."""
    print(f"\n=== Test 2: LIMIT (k={k_val}) ===")
    limit_n = 5
    stdout, stderr, rc = run_binary([
        "--shards", shard_dir,
        "--gc-hist", gc_hist,
        "--limit", str(limit_n),
        "--random_access",
    ])
    if rc != 0:
        check(False, "binary exited non-zero", f"rc={rc} {stderr[:200]}")
        return

    results = get_results(stdout)
    check(len(results) <= limit_n,
          f"at most {limit_n} results (got {len(results)})")


def test_gc_filter(k_val, shard_dir, gc_hist):
    """Verify --gc-min / --gc-max range is respected."""
    print(f"\n=== Test 3: GC FILTER (k={k_val}) ===")
    gc_min, gc_max = 40, 60
    stdout, stderr, rc = run_binary([
        "--shards", shard_dir,
        "--gc-hist", gc_hist,
        "--limit", "10",
        "--random_access",
        "--gc-min", str(gc_min),
        "--gc-max", str(gc_max),
    ])
    if rc != 0:
        check(False, "binary exited non-zero", f"rc={rc} {stderr[:200]}")
        return

    results = get_results(stdout)
    if not results:
        check(False, "no results returned for GC filter 40-60")
        return

    all_in_range = True
    violators = []
    for kmer in results:
        gc = gc_content(kmer)
        if gc < gc_min or gc > gc_max:
            all_in_range = False
            violators.append((kmer, gc))

    check(all_in_range,
          f"all {len(results)} k-mers have GC% in [{gc_min},{gc_max}]",
          f"violations: {violators}" if violators else "")


def test_motif_flag(k_val, shard_dir, gc_hist):
    """Verify --motif-mode flag --filter-homopolymers adds motif info to __META__."""
    print(f"\n=== Test 4: MOTIF FLAG (k={k_val}) ===")
    stdout, stderr, rc = run_binary([
        "--shards", shard_dir,
        "--gc-hist", gc_hist,
        "--limit", "10",
        "--random_access",
        "--motif-mode", "flag",
        "--filter-homopolymers",
        "--include-motif-metadata",
    ])
    if rc != 0:
        check(False, "binary exited non-zero", f"rc={rc} {stderr[:200]}")
        return

    meta = parse_meta(stdout)
    check(meta is not None, "__META__ line present")
    if meta:
        has_motif_info = "motif_filter_applied" in meta["motif_info"]
        check(has_motif_info, "__META__ includes motif info",
              f"got motif_info={meta['motif_info']!r}" if not has_motif_info else "")

    # Also check that [INFO] lines mention motif
    has_motif_info_lines = "Motif mode" in stdout or "motif" in stdout.lower()
    check(has_motif_info_lines, "output contains motif-related info lines")


def test_output_has_kmers(k_val, shard_dir, gc_hist):
    """Verify each data line is a valid DNA string (A/C/G/T only)."""
    print(f"\n=== Test 5: OUTPUT HAS VALID KMERS (k={k_val}) ===")
    stdout, stderr, rc = run_binary([
        "--shards", shard_dir,
        "--gc-hist", gc_hist,
        "--limit", "20",
        "--random_access",
    ])
    if rc != 0:
        check(False, "binary exited non-zero", f"rc={rc} {stderr[:200]}")
        return

    results = get_results(stdout)
    if not results:
        check(False, "no results returned")
        return

    valid_bases = set("ACGT")
    all_valid = True
    bad_kmers = []
    for kmer in results:
        if not all(c in valid_bases for c in kmer.upper()):
            all_valid = False
            bad_kmers.append(kmer)

    check(all_valid,
          f"all {len(results)} k-mers contain only A/C/G/T",
          f"invalid: {bad_kmers}" if bad_kmers else "")


def test_cursor(k_val, shard_dir, gc_hist):
    """Verify cursor-based pagination returns different results."""
    print(f"\n=== Test 6: CURSOR TEST (k={k_val}) ===")
    # First page
    stdout1, stderr1, rc1 = run_binary([
        "--shards", shard_dir,
        "--gc-hist", gc_hist,
        "--limit", "5",
        "--random_access",
    ])
    if rc1 != 0:
        check(False, "first page: binary exited non-zero", f"rc={rc1} {stderr1[:200]}")
        return

    meta1 = parse_meta(stdout1)
    if meta1 is None:
        check(False, "first page: no __META__ line")
        return

    cursor = meta1["cursor"]
    if not cursor:
        check(False, "first page: cursor is empty")
        return

    results1 = get_results(stdout1)
    if not results1:
        check(False, "first page: no results")
        return

    # Second page with cursor
    stdout2, stderr2, rc2 = run_binary([
        "--shards", shard_dir,
        "--gc-hist", gc_hist,
        "--limit", "5",
        "--random_access",
        "--cursor", cursor,
    ])
    if rc2 != 0:
        check(False, "second page: binary exited non-zero", f"rc={rc2} {stderr2[:200]}")
        return

    results2 = get_results(stdout2)

    # Results should be different (no overlap expected)
    overlap = set(results1) & set(results2)
    all_different = len(overlap) == 0

    check(all_different,
          "second page returns different k-mers than first",
          f"overlap={overlap}" if not all_different else "")


# ===================================================================
#  Main
# ===================================================================
def main():
    # Check binary exists
    if not os.path.exists(BINARY):
        print(f"ERROR: binary not found at {BINARY}")
        sys.exit(1)

    # --- REQUIRED TESTS (always run) ---
    test_required_usage_errors()

    # --- OPTIONAL TESTS (need production shards) ---
    k_val, shard_dir, gc_hist = pick_shards()
    if k_val is None:
        print("\n---")
        print("[SKIP] Production shards not found — skipping shard-based tests.")
        print("To run full tests, ensure shards_16/ (or 17/18) and")
        print("corresponding gc_hist_shards_*.json are present at the")
        print("repository root, or set BARCODEDB_ROOT.")
    else:
        print(f"\n--- Found shards: k={k_val} at {shard_dir} ---")
        test_basic_search(k_val, shard_dir, gc_hist)
        test_limit(k_val, shard_dir, gc_hist)
        test_gc_filter(k_val, shard_dir, gc_hist)
        test_motif_flag(k_val, shard_dir, gc_hist)
        test_output_has_kmers(k_val, shard_dir, gc_hist)
        test_cursor(k_val, shard_dir, gc_hist)

    # --- Summary ---
    total = g_passes + g_failures
    print(f"\n{'=' * 50}")
    print(f"Results: {g_passes} passed, {g_failures} failed out of {total} tests")
    if g_failures > 0:
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    main()
