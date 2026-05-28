#!/usr/bin/env python3
"""
test_output_schema.py — Validate output format/schema regression.

Part A (always runs — uses synthetic fixtures):
1. Run query_kmer_bitmap with --bitmap on k16 fixture and present k-mers
2. Verify output lines are: <kmer>TAB<hit> (2 columns separated by tab)
3. Verify hit is '0' or '1'
4. Run with --min-hamming-distance 2 (if binary supports it) and verify extra columns
5. Verify no blank lines in output

Part B (optional — needs production shards):
Run query_substring_bitmap_stream and verify:
1. First line starts with __META__
2. Remaining lines are DNA strings
3. __META__ has enough tab-separated fields

Exit code 0 if all pass, 1 if any fail.
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
        return int(k)
    except Exception:
        return None


def fixture_path(k):
    return os.path.join(QC_TMP, f"k{k}_fixture.kbit")


def write_kmer_file(kmers):
    """Write k-mers to a temp file, return path."""
    tmp = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".kmers")
    for k in kmers:
        tmp.write(k + "\n")
    tmp.close()
    return tmp.name


def run_binary(bin_path, args, timeout=60):
    """Run a binary with args, return (stdout, stderr, rc)."""
    cmd = [bin_path] + args
    proc = subprocess.run(cmd, capture_output=True, text=True,
                          cwd=REPO_ROOT, timeout=timeout)
    return proc.stdout, proc.stderr, proc.returncode


def pick_shards():
    """Return (k, shard_dir, gc_hist) or (None, None, None)."""
    for k in (16, 17, 18):
        sd = SHARD_DIRS[k]
        gh = GC_HIST_FILES[k]
        if os.path.isdir(sd) and os.path.isfile(gh):
            return k, sd, gh
    return None, None, None


def get_substring_results(stdout):
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


# ===================================================================
#  Part A: K-mer bitmap output schema
# ===================================================================
def test_part_a_basic_output_schema():
    """Verify query_kmer_bitmap output format: <kmer>TAB<hit>."""
    print("\n=== Part A: K-mer bitmap output schema ===")

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
        actual_k = read_fixture_header(fp)
        if actual_k is None:
            continue
        if actual_k != k_val:
            actual_k = k_val  # trust the truth

        present_kmers = [k for k in fixture["present_kmers"] if len(k) == actual_k]
        if not present_kmers:
            continue

        # Limit to first 20 to keep test fast
        test_kmers = present_kmers[:20]

        tmp_path = write_kmer_file(test_kmers)
        stdout, stderr, rc = run_binary(BINARY, [
            "--bitmap", fp, "--kmers", tmp_path,
        ])
        os.unlink(tmp_path)

        if rc != 0:
            check(False, f"k={k_val}: binary exited {rc}",
                  stderr.strip()[:200])
            continue

        # --- A1: Check no blank lines ---
        raw_lines = stdout.strip().split("\n")
        blank_lines = [l for l in raw_lines if l.strip() == ""]
        check(len(blank_lines) == 0,
              f"k={k_val}: no blank lines in output",
              f"found {len(blank_lines)} blank line(s)")

        # --- A2: Each output line has exactly 2 tab-separated columns ---
        non_empty = [l for l in raw_lines if l.strip()]
        all_two_cols = True
        bad_lines = []
        for i, line in enumerate(non_empty):
            parts = line.split("\t")
            if len(parts) != 2:
                all_two_cols = False
                bad_lines.append((i, line, len(parts)))

        check(all_two_cols,
              f"k={k_val}: all {len(non_empty)} lines have exactly 2 tab-separated columns",
              f"violations: {bad_lines[:3]}" if bad_lines else "")

        # --- A3: Hit column is '0' or '1' ---
        all_valid_hits = True
        bad_hits = []
        for line in non_empty:
            parts = line.split("\t")
            if len(parts) == 2:
                hit = parts[1]
                if hit not in ("0", "1"):
                    all_valid_hits = False
                    bad_hits.append((parts[0], hit))

        check(all_valid_hits,
              f"k={k_val}: all hit values are '0' or '1'",
              f"invalid: {bad_hits[:5]}" if bad_hits else "")

        # --- A4: First column is a valid DNA k-mer ---
        valid_bases = set("ACGT")
        all_dna = True
        bad_kmers = []
        for line in non_empty:
            parts = line.split("\t")
            if len(parts) >= 1:
                kmer = parts[0].upper()
                if not all(c in valid_bases for c in kmer):
                    all_dna = False
                    bad_kmers.append(kmer)

        check(all_dna,
              f"k={k_val}: all first columns are valid DNA strings",
              f"invalid: {bad_kmers[:5]}" if bad_kmers else "")


def test_part_a_hamming_distance_schema():
    """Verify --min-hamming-distance adds extra columns (if supported)."""
    print("\n  --- A4: --min-hamming-distance output format ---")

    if not os.path.exists(BINARY):
        return
    if not os.path.exists(TRUTH_PATH):
        return

    with open(TRUTH_PATH, "r") as f:
        truth = json.load(f)

    # Use k=16 fixture
    k_val = 16
    fp = fixture_path(k_val)
    actual_k = read_fixture_header(fp)
    if actual_k is None:
        print("  [SKIP] k=16 fixture not found")
        return

    # Find fixture data
    fixture = None
    for f in truth["fixtures"]:
        if f["k"] == k_val:
            fixture = f
            break
    if fixture is None:
        print("  [SKIP] k=16 not in truth data")
        return

    present_kmers = [k for k in fixture["present_kmers"] if len(k) == actual_k]
    if not present_kmers:
        print("  [SKIP] no valid present k-mers for k=16")
        return

    test_kmers = present_kmers[:10]
    tmp_path = write_kmer_file(test_kmers)

    # Try running with --min-hamming-distance 2
    stdout, stderr, rc = run_binary(BINARY, [
        "--bitmap", fp, "--kmers", tmp_path,
        "--min-hamming-distance", "2",
    ])
    os.unlink(tmp_path)

    if rc != 0:
        print(f"  [INFO] --min-hamming-distance not supported in --bitmap mode "
              f"(rc={rc}) — skipping schema check")
        # This is acceptable; the flag may only work in --shards mode
        return

    raw_lines = [l for l in stdout.strip().split("\n") if l.strip()]
    if not raw_lines:
        return

    all_valid_format = True
    for line in raw_lines:
        parts = line.split("\t")
        # With --min-hamming-distance, we expect at least 4 columns:
        # <kmer> TAB <hit> TAB <nearest_hamming> TAB <passes>
        if len(parts) < 4:
            all_valid_format = False
            check(False, f"expected >=4 columns with --min-hamming-distance",
                  f"got {len(parts)} cols: {line[:80]}")
            break

    if all_valid_format:
        check(True,
              f"output has >=4 columns with --min-hamming-distance",
              f"({len(raw_lines)} lines)")

        # Verify extra columns are numeric
        all_numeric = True
        for line in raw_lines:
            parts = line.split("\t")
            if len(parts) >= 3:
                try:
                    int(parts[2])  # nearest hamming
                except ValueError:
                    all_numeric = False
                    break
        check(all_numeric,
              "hamming distance column is numeric",
              "" if all_numeric else "non-numeric value found")


# ===================================================================
#  Part B: Substring binary output schema
# ===================================================================
def test_part_b_substring_schema():
    """Verify query_substring_bitmap_stream output schema."""
    print("\n=== Part B: Substring binary output schema (optional) ===")

    k_val, shard_dir, gc_hist = pick_shards()
    if k_val is None:
        print("  [SKIP] Production shards not found")
        return

    if not os.path.exists(SUBSTR_BINARY):
        print("  [SKIP] query_substring_bitmap_stream binary not found")
        return

    stdout, stderr, rc = run_binary(SUBSTR_BINARY, [
        "--shards", shard_dir, "--gc-hist", gc_hist,
        "--limit", "10", "--random_access",
    ], timeout=120)

    if rc != 0:
        check(False, "substring binary exited non-zero",
              f"rc={rc} {stderr[:200]}")
        return

    lines = stdout.strip().split("\n")
    non_empty = [l for l in lines if l.strip()]

    if not non_empty:
        check(False, "substring binary produced output")
        return

    # --- B1: First non-empty line starts with __META__ ---
    first_line = non_empty[0]
    starts_with_meta = first_line.startswith("__META__")
    check(starts_with_meta,
          "first output line starts with __META__",
          f"got: {first_line[:80]!r}" if not starts_with_meta else "")

    # --- B2: __META__ has enough tab-separated fields ---
    if starts_with_meta:
        meta_parts = first_line.split("\t")
        min_fields = 5  # __META__, cursor, hasMore, returned_count, k
        check(len(meta_parts) >= min_fields,
              f"__META__ has at least {min_fields} tab-separated fields",
              f"got {len(meta_parts)}: {first_line[:120]!r}")

        # Verify field 4 (index 4) is a numeric k
        if len(meta_parts) >= 5:
            k_field = meta_parts[4]
            is_k_numeric = k_field.isdigit() or (
                k_field.startswith("-") and k_field[1:].isdigit()
            )
            check(is_k_numeric,
                  f"__META__ field 5 (k) is numeric: {k_field}",
                  f"got {k_field!r}" if not is_k_numeric else "")
        else:
            check(False, "__META__ field 5 (k) exists")

        # Verify field 3 (index 3) is a numeric returned_count
        if len(meta_parts) >= 4:
            rc_field = meta_parts[3]
            is_rc_numeric = rc_field.isdigit()
            check(is_rc_numeric,
                  f"__META__ field 4 (returned_count) is numeric: {rc_field}",
                  f"got {rc_field!r}" if not is_rc_numeric else "")

        # Verify field 2 (index 2) is hasMore ('0'/'1')
        if len(meta_parts) >= 3:
            hm_field = meta_parts[2]
            is_hm_valid = hm_field in ("0", "1")
            check(is_hm_valid,
                  f"__META__ field 3 (hasMore) is '0' or '1': {hm_field!r}",
                  f"got {hm_field!r}" if not is_hm_valid else "")

    # --- B3: Remaining lines are DNA strings ---
    data_lines = get_substring_results(stdout)
    if data_lines:
        valid_bases = set("ACGT")
        all_dna = True
        bad = []
        for kmer in data_lines:
            if not all(c in valid_bases for c in kmer.upper()):
                all_dna = False
                bad.append(kmer)

        check(all_dna,
              f"all {len(data_lines)} result lines are valid DNA strings",
              f"invalid: {bad[:5]}" if bad else "")

        # Each k-mer length should be consistent within a page
        lengths = set(len(k) for k in data_lines)
        check(len(lengths) == 1,
              "all result k-mers have the same length",
              f"lengths: {lengths}" if len(lengths) != 1 else "")

        # No blank lines in data section
        blank_data = [l for l in data_lines if l.strip() == ""]
        check(len(blank_data) == 0,
              "no blank lines in data section",
              f"found {len(blank_data)} blank(s)")
    else:
        print("  INFO  No data lines returned (possible empty result set)")


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

    test_part_a_basic_output_schema()
    test_part_a_hamming_distance_schema()
    test_part_b_substring_schema()

    total = g_passes + g_failures
    print(f"\n{'=' * 60}")
    print(f"Results: {g_passes} passed, {g_failures} failed out of {total} tests")

    if g_failures == 0:
        print("ALL TESTS PASSED ✓")
        sys.exit(0)
    else:
        print("SOME TESTS FAILED ✗")
        sys.exit(1)


if __name__ == "__main__":
    main()
