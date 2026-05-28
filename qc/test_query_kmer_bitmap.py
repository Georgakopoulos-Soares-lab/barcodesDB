#!/usr/bin/env python3
"""
test_query_kmer_bitmap.py — Validate query_kmer_bitmap binary.

Run from the repository root (barcodes/):
    python3 website/qc/test_query_kmer_bitmap.py

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
# Paths (relative to repository root)
# ---------------------------------------------------------------------------
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
QC_TMP = os.path.join(REPO_ROOT, "website", "qc", "tmp")
BINARY = os.path.join(QC_TMP, "query_kmer_bitmap")
TRUTH_PATH = os.path.join(QC_TMP, "synthetic_truth.json")

# ---------------------------------------------------------------------------
# Globals
# ---------------------------------------------------------------------------
g_truth = None
g_failures = 0
g_passes = 0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_truth():
    """Read synthetic_truth.json and return parsed dict."""
    global g_truth
    if g_truth is not None:
        return g_truth
    with open(TRUTH_PATH, "r") as f:
        g_truth = json.load(f)
    return g_truth


def read_fixture_header(fixture_path):
    """Read KBITv1 header and return (k, flags) or None on failure."""
    if not os.path.exists(fixture_path):
        return None
    try:
        with open(fixture_path, "rb") as f:
            hdr = f.read(64)
        if len(hdr) < 64:
            return None
        if hdr[:6] != b"KBITv1":  # magic is "KBITv1\0\0" (8 bytes)
            return None
        k = struct.unpack_from("<Q", hdr, 24)[0]
        flags = struct.unpack_from("<Q", hdr, 40)[0]
        return (int(k), int(flags))
    except Exception:
        return None


def fixture_path(k):
    """Return the expected KBITv1 fixture path for a given k."""
    return os.path.join(QC_TMP, f"k{k}_fixture.kbit")


def run_binary(args, stdin_data=None):
    """
    Run query_kmer_bitmap with the given list of CLI arguments.

    Args:
        args: List of extra CLI arguments.
        stdin_data: Optional bytes/string to pipe to stdin (not used here).

    Returns:
        (stdout, stderr, returncode)
    """
    cmd = [BINARY] + args
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
        timeout=60,
    )
    return proc.stdout, proc.stderr, proc.returncode


def write_kmer_file(kmers):
    """Write a list of k-mer strings to a temporary file.  Returns the path."""
    tmp = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".kmers")
    for k in kmers:
        tmp.write(k + "\n")
    tmp.close()
    return tmp.name


def parse_output(stdout):
    """
    Parse tab-delimited output lines into a list of dictionaries.

    Without motif flags:  <kmer>TAB<hit:0|1>
    With --motif-mode:    <kmer>TAB<hit>TAB<passes>TAB<hit_count>TAB<hits_detail>

    Returns [{"kmer": ..., "hit": "0"|"1", ...}, ...]
    """
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


# ---------------------------------------------------------------------------
# Test runners
# ---------------------------------------------------------------------------
def check(condition, test_name, detail=""):
    """Record a PASS or FAIL for the given test."""
    global g_passes, g_failures
    if condition:
        g_passes += 1
        print(f"  PASS  {test_name}")
    else:
        g_failures += 1
        print(f"  FAIL  {test_name}  {detail}")


# ===================================================================
#  1. KNOWN-PRESENT
# ===================================================================
def test_known_present():
    """Verify present kmers (matching fixture k) return hit=1."""
    print("\n=== Test 1: KNOWN-PRESENT k-mers ===")
    truth = load_truth()
    for fixture in truth["fixtures"]:
        k_val = fixture["k"]
        fp = fixture_path(k_val)
        hdr = read_fixture_header(fp)
        if hdr is None:
            check(False, f"k={k_val}: cannot read fixture", fp)
            continue
        actual_k = hdr[0]

        present_kmers = fixture["present_kmers"]
        # Filter to only those matching the actual fixture k
        valid_kmers = [k for k in present_kmers if len(k) == actual_k]
        skipped = [k for k in present_kmers if len(k) != actual_k]
        for k in skipped:
            print(f"  INFO  k={k_val}: skipping k-mer len={len(k)} (fixture k={actual_k}) [{k}]")

        if not valid_kmers:
            check(False, f"k={k_val}: no valid present kmers match k={actual_k}")
            continue

        tmp_path = write_kmer_file(valid_kmers)
        stdout, stderr, rc = run_binary([
            "--bitmap", fp,
            "--kmers", tmp_path,
        ])
        os.unlink(tmp_path)

        if rc != 0:
            check(False, f"k={k_val}: binary exited {rc}", stderr.strip()[:200])
            continue

        results = parse_output(stdout)
        if len(results) != len(valid_kmers):
            check(False, f"k={k_val}: expected {len(valid_kmers)} lines, got {len(results)}")
            continue

        all_present = all(r["hit"] == "1" for r in results)
        check(all_present, f"k={k_val}: all {len(valid_kmers)} present k-mers returned 1",
              "" if all_present else str([(r["kmer"], r["hit"]) for r in results if r["hit"] != "1"]))


# ===================================================================
#  2. KNOWN-ABSENT
# ===================================================================
def test_known_absent():
    """Verify absent kmers (matching fixture k) return hit=0."""
    print("\n=== Test 2: KNOWN-ABSENT k-mers ===")
    truth = load_truth()
    for fixture in truth["fixtures"]:
        k_val = fixture["k"]
        fp = fixture_path(k_val)
        hdr = read_fixture_header(fp)
        if hdr is None:
            check(False, f"k={k_val}: cannot read fixture", fp)
            continue
        actual_k = hdr[0]

        absent_kmers = fixture["absent_kmers"]
        valid_kmers = [k for k in absent_kmers if len(k) == actual_k]
        skipped = [k for k in absent_kmers if len(k) != actual_k]
        for k in skipped:
            print(f"  INFO  k={k_val}: skipping k-mer len={len(k)} (fixture k={actual_k}) [{k}]")

        if not valid_kmers:
            check(False, f"k={k_val}: no valid absent kmers match k={actual_k}")
            continue

        tmp_path = write_kmer_file(valid_kmers)
        stdout, stderr, rc = run_binary([
            "--bitmap", fp,
            "--kmers", tmp_path,
        ])
        os.unlink(tmp_path)

        if rc != 0:
            check(False, f"k={k_val}: binary exited {rc}", stderr.strip()[:200])
            continue

        results = parse_output(stdout)
        if len(results) != len(valid_kmers):
            check(False, f"k={k_val}: expected {len(valid_kmers)} lines, got {len(results)}")
            continue

        all_absent = all(r["hit"] == "0" for r in results)
        check(all_absent, f"k={k_val}: all {len(valid_kmers)} absent k-mers returned 0",
              "" if all_absent else str([(r["kmer"], r["hit"]) for r in results if r["hit"] != "0"]))


# ===================================================================
#  3. LOWERCASE
# ===================================================================
def test_lowercase():
    """Test that lowercase inputs are handled."""
    print("\n=== Test 3: LOWERCASE input handling ===")
    truth = load_truth()
    for fixture in truth["fixtures"]:
        k_val = fixture["k"]
        fp = fixture_path(k_val)
        hdr = read_fixture_header(fp)
        if hdr is None:
            continue
        actual_k = hdr[0]

        present_kmers = fixture["present_kmers"]
        valid_kmers = [k for k in present_kmers if len(k) == actual_k]
        if not valid_kmers:
            continue

        lowercase_kmers = [k.lower() for k in valid_kmers]
        tmp_path = write_kmer_file(lowercase_kmers)
        stdout, stderr, rc = run_binary([
            "--bitmap", fp,
            "--kmers", tmp_path,
        ])
        os.unlink(tmp_path)

        if rc == 0:
            results = parse_output(stdout)
            if len(results) == len(lowercase_kmers):
                all_correct = all(r["hit"] == "1" for r in results)
                check(all_correct,
                      f"k={k_val}: lowercase inputs accepted, all hits=1",
                      "" if all_correct else str([(r["kmer"], r["hit"]) for r in results]))
            else:
                check(False, f"k={k_val}: lowercase accepted but wrong line count",
                      f"expected {len(lowercase_kmers)} got {len(results)}")
        else:
            # Binary rejected lowercase — that's also acceptable; just document
            print(f"  INFO  k={k_val}: binary exited {rc} on lowercase input "
                  f"(stderr: {stderr.strip()[:120]})")


# ===================================================================
#  4. INVALID CHARACTERS
# ===================================================================
def test_invalid_characters():
    """Test invalid nucleotide characters and empty lines are handled safely."""
    print("\n=== Test 4: INVALID CHARACTERS ===")
    truth = load_truth()
    test_cases = [
        ("N-containing", "ACGTACGTACGTACGN"),  # 16 chars with trailing N
        ("X-containing", "ACGTACGTACGTACGX"),  # 16 chars with trailing X
    ]
    test_case_empty = ("empty line", "")

    for fixture in truth["fixtures"]:
        k_val = fixture["k"]
        fp = fixture_path(k_val)
        hdr = read_fixture_header(fp)
        if hdr is None:
            continue
        actual_k = hdr[0]

        # --- Test N/X with k-mer sized to match fixture k ---
        for label, template in test_cases:
            # Adjust template length to match fixture k
            if len(template) > actual_k:
                kmer = template[:actual_k]
            elif len(template) < actual_k:
                # Pad with 'A'
                kmer = template + "A" * (actual_k - len(template))
            else:
                kmer = template

            tmp_path = write_kmer_file([kmer])
            stdout, stderr, rc = run_binary([
                "--bitmap", fp,
                "--kmers", tmp_path,
            ])
            os.unlink(tmp_path)

            check(rc != 0,
                  f"k={k_val}, {label}: binary exits nonzero (rc={rc})",
                  stderr.strip()[:200] if stderr else "no stderr")

            if rc != 0:
                has_error_msg = len(stderr.strip()) > 0
                check(has_error_msg,
                      f"k={k_val}, {label}: error message present",
                      "" if has_error_msg else "missing stderr")

        # --- Test empty line (just a newline): binary silently skips it ---
        tmp_path = write_kmer_file([test_case_empty[1]])
        stdout, stderr, rc = run_binary([
            "--bitmap", fp,
            "--kmers", tmp_path,
        ])
        os.unlink(tmp_path)

        # Empty lines are skipped by the binary (FastLineReader skips empty lines)
        check(rc == 0,
              f"k={k_val}, empty line: binary exits 0 (skips empty lines)",
              f"rc={rc}, stderr={stderr.strip()[:200]}")
        # Output should be empty (no lines for empty input)
        check(len(stdout.strip()) == 0,
              f"k={k_val}, empty line: no output",
              f"got {len(stdout.strip())} chars")


# ===================================================================
#  5. WRONG LENGTH
# ===================================================================
def test_wrong_length():
    """Test shorter and longer k-mers are rejected."""
    print("\n=== Test 5: WRONG LENGTH k-mers ===")
    truth = load_truth()
    for fixture in truth["fixtures"]:
        k_val = fixture["k"]
        fp = fixture_path(k_val)
        hdr = read_fixture_header(fp)
        if hdr is None:
            continue
        actual_k = hdr[0]

        shorter = "ACGT"  # always too short
        # longer: use a length that is clearly different
        longer = "ACGT" * 10  # 40 chars

        for label, kmer in [("shorter", shorter), ("longer", longer)]:
            tmp_path = write_kmer_file([kmer])
            stdout, stderr, rc = run_binary([
                "--bitmap", fp,
                "--kmers", tmp_path,
            ])
            os.unlink(tmp_path)

            check(rc != 0,
                  f"k={k_val}, {label}: binary exits nonzero (rc={rc})",
                  stderr.strip()[:200] if stderr else "no stderr")


# ===================================================================
#  6. MISSING FIXTURE
# ===================================================================
def test_missing_fixture():
    """Test that a nonexistent bitmap path produces an error."""
    print("\n=== Test 6: MISSING FIXTURE ===")
    nonexistent = os.path.join(QC_TMP, "nonexistent_fixture.kbit")
    tmp_path = write_kmer_file(["ACGTACGTACGTACGT"])
    stdout, stderr, rc = run_binary([
        "--bitmap", nonexistent,
        "--kmers", tmp_path,
    ])
    os.unlink(tmp_path)

    check(rc != 0,
          f"nonexistent bitmap: binary exits nonzero (rc={rc})",
          stderr.strip()[:200] if stderr else "no stderr")
    if rc != 0:
        has_error_msg = len(stderr.strip()) > 0
        check(has_error_msg,
              "nonexistent bitmap: error message present",
              "" if has_error_msg else "missing stderr")


# ===================================================================
#  7. MOTIF FILTERS (preserved)
# ===================================================================
def test_motif_filters():
    """Verify motif filter flags add extra columns but do not alter hit status."""
    print("\n=== Test 7: MOTIF FILTERS ===")

    # Use k=16 fixture
    k_val = 16
    fp = fixture_path(k_val)
    hdr = read_fixture_header(fp)
    if hdr is None:
        check(False, "k=16 fixture not found", fp)
        return
    actual_k = hdr[0]

    truth = load_truth()
    fixture = None
    for f in truth["fixtures"]:
        if f["k"] == k_val:
            fixture = f
            break
    if fixture is None:
        check(False, f"k={k_val} not in truth data")
        return

    present_kmers = [k for k in fixture["present_kmers"] if len(k) == actual_k]
    if not present_kmers:
        check(False, f"k={k_val}: no valid present kmers")
        return

    # 7a. Run WITHOUT motif flags — verify output format (kmer TAB hit)
    print("  7a. Without motif flags — basic output format")
    tmp_path = write_kmer_file(present_kmers)
    stdout, stderr, rc = run_binary([
        "--bitmap", fp,
        "--kmers", tmp_path,
    ])
    os.unlink(tmp_path)

    check(rc == 0, "binary exits 0", stderr.strip()[:200] if stderr else "")

    if rc == 0:
        lines = [l for l in stdout.strip().split("\n") if l.strip()]
        all_two_cols = all(len(l.split("\t")) == 2 for l in lines)
        check(all_two_cols,
              f"output has exactly 2 tab-separated columns ({len(lines)} lines)",
              "" if all_two_cols else str([l for l in lines if len(l.split("\t")) != 2]))

        results = parse_output(stdout)
        all_hit = all(r["hit"] == "1" for r in results)
        check(all_hit, "all present k-mers still show hit=1 without motif flags",
              "" if all_hit else str([(r["kmer"], r["hit"]) for r in results]))

    # 7b. Run WITH --motif-mode flag --filter-homopolymers on all present kmers
    print("  7b. With --motif-mode flag --filter-homopolymers (all present kmers)")
    tmp_path = write_kmer_file(present_kmers)
    stdout, stderr, rc = run_binary([
        "--bitmap", fp,
        "--kmers", tmp_path,
        "--motif-mode", "flag",
        "--filter-homopolymers",
    ])
    os.unlink(tmp_path)

    check(rc == 0, "binary exits 0", stderr.strip()[:200] if stderr else "")

    if rc == 0:
        results = parse_output(stdout)
        check(len(results) == len(present_kmers),
              f"output has {len(present_kmers)} result lines",
              f"got {len(results)}")

        # Verify all present kmers still have hit=1
        all_hit = all(r["hit"] == "1" for r in results)
        check(all_hit,
              "all present k-mers still show hit=1 with motif flags",
              "" if all_hit else str([(r["kmer"], r["hit"]) for r in results if r["hit"] != "1"]))

        # Verify motif metadata columns present
        has_motif_cols = all("motif_passes" in r for r in results)
        check(has_motif_cols,
              "motif_passes column present for all lines",
              "" if has_motif_cols else str([r["kmer"] for r in results if "motif_passes" not in r]))

        # In flag mode, motif_passes should be '1' for all (just flagging)
        all_pass = all(r.get("motif_passes") == "1" for r in results)
        check(all_pass,
              "motif_passes = 1 for all (flag mode does not exclude)",
              "" if all_pass else str([(r["kmer"], r.get("motif_passes")) for r in results]))

    # 7c. Test with a known homopolymer k-mer (absent from bitmap is fine)
    print("  7c. With --motif-mode flag on a homopolymer k-mer (absent)")
    # Create a 16-mer that exceeds homopolymer threshold: 6 As followed by 10 Ts
    homopolymer_kmer_absent = "AAAAAATTTTTTTTTT"  # 16 chars, run of 6 A's (>4)
    tmp_path = write_kmer_file([homopolymer_kmer_absent, present_kmers[0]])
    stdout, stderr, rc = run_binary([
        "--bitmap", fp,
        "--kmers", tmp_path,
        "--motif-mode", "flag",
        "--filter-homopolymers",
    ])
    os.unlink(tmp_path)

    check(rc == 0, "binary exits 0", stderr.strip()[:200] if stderr else "")

    if rc == 0:
        results = parse_output(stdout)
        check(len(results) == 2, "output has 2 result lines", f"got {len(results)}")
        if len(results) >= 2:
            r_hp = results[0]  # first result: the homopolymer k-mer
            # Hit status is independent of motif
            hit_unchanged = r_hp["hit"] in ("0", "1")
            check(hit_unchanged, "hit status present (unchanged by motif mode)",
                  f"hit={r_hp['hit']}")
            # Should have motif hit count >= 1 (homopolymer detected)
            hp_count = int(r_hp.get("motif_hit_count", "0"))
            check(hp_count >= 1,
                  f"homopolymer detected (motif_hit_count={hp_count})",
                  f"got motif_hit_count={hp_count}")
            motif_hits = r_hp.get("motif_hits", "")
            check("homopolymer" in motif_hits.lower(),
                  f"motif hits detail mentions 'homopolymer'",
                  f"got: {motif_hits}")


# ===================================================================
#  8. OUTPUT STABILITY
# ===================================================================
def test_output_stability():
    """Verify that running the same query twice produces identical output."""
    print("\n=== Test 8: OUTPUT STABILITY ===")
    truth = load_truth()
    for fixture in truth["fixtures"]:
        k_val = fixture["k"]
        fp = fixture_path(k_val)
        hdr = read_fixture_header(fp)
        if hdr is None:
            continue
        actual_k = hdr[0]

        all_kmers = (
            [k for k in fixture["present_kmers"] if len(k) == actual_k]
            + [k for k in fixture["absent_kmers"] if len(k) == actual_k]
        )
        if not all_kmers:
            check(False, f"k={k_val}: no valid kmers for stability test")
            continue

        tmp_path = write_kmer_file(all_kmers)

        stdout1, _, rc1 = run_binary([
            "--bitmap", fp,
            "--kmers", tmp_path,
        ])
        stdout2, _, rc2 = run_binary([
            "--bitmap", fp,
            "--kmers", tmp_path,
        ])
        os.unlink(tmp_path)

        both_ok = (rc1 == 0 and rc2 == 0)
        check(both_ok, f"k={k_val}: both runs exit 0",
              f"rc1={rc1}, rc2={rc2}")

        if both_ok:
            identical = (stdout1 == stdout2)
            check(identical, f"k={k_val}: output identical across two runs",
                  "" if identical else "output differs between runs")
        elif rc1 == 0:
            check(False, f"k={k_val}: first run OK (rc=0), second failed (rc={rc2})")
        elif rc2 == 0:
            check(False, f"k={k_val}: first run failed (rc={rc1}), second OK (rc=0)")
        # else both failed — that's ok


# ===================================================================
# Main
# ===================================================================
def main():
    global g_failures, g_passes
    g_failures = 0
    g_passes = 0

    print(f"Repository root: {REPO_ROOT}")
    print(f"Binary:          {BINARY}")
    print(f"Truth file:      {TRUTH_PATH}")

    if not os.path.exists(BINARY):
        print(f"\nERROR: Binary not found at {BINARY}")
        print("Run 'bash website/qc/build_query_binaries.sh' first.")
        sys.exit(1)

    if not os.path.exists(TRUTH_PATH):
        print(f"\nERROR: Truth file not found at {TRUTH_PATH}")
        print("Ensure build_query_binaries.sh has been run to generate fixtures.")
        sys.exit(1)

    load_truth()
    truth = g_truth
    print(f"Found {len(truth['fixtures'])} fixture(s): "
          + ", ".join(f"k={f['k']}" for f in truth["fixtures"]))

    test_known_present()
    test_known_absent()
    test_lowercase()
    test_invalid_characters()
    test_wrong_length()
    test_missing_fixture()
    test_motif_filters()
    test_output_stability()

    print(f"\n{'='*60}")
    print(f"Results: {g_passes} passed, {g_failures} failed")

    if g_failures == 0:
        print("ALL TESTS PASSED ✓")
        sys.exit(0)
    else:
        print("SOME TESTS FAILED ✗")
        sys.exit(1)


if __name__ == "__main__":
    main()
