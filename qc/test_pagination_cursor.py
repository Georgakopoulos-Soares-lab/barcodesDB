#!/usr/bin/env python3
"""
test_pagination_cursor.py — Validate cursor-based pagination behavior.

Requires production shards (optional — skips gracefully if not present).

Tests:
1. Run with --limit 5, capture cursor and hasMore from __META__
2. Run with --cursor <cursor> --limit 5, verify results differ
3. Collect all pages and verify no duplicate k-mers across pages
4. Verify metadata consistency (returned_count matches actual row count)

Exit code 0 always (optional tests).
Dependencies: Python standard library only.
"""

import os
import subprocess
import sys

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
QC_DIR = os.path.abspath(os.path.dirname(__file__))
BINARY = os.path.join(QC_DIR, "tmp", "query_substring_bitmap_stream")

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


def pick_shards():
    """Return (k, shard_dir, gc_hist) or (None, None, None)."""
    for k in (16, 17, 18):
        sd = SHARD_DIRS[k]
        gh = GC_HIST_FILES[k]
        if os.path.isdir(sd) and os.path.isfile(gh):
            return k, sd, gh
    return None, None, None


def get_results(stdout):
    """Return data lines from stdout (after __META__, before [INFO])."""
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
    """Parse __META__ line. Returns dict or None."""
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


def run_substr(args, timeout=120):
    """Run query_substring_bitmap_stream with given args."""
    cmd = [BINARY] + args
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    return proc.stdout, proc.stderr, proc.returncode


# ===================================================================
#  Tests
# ===================================================================
def test_cursor_pagination(k_val, shard_dir, gc_hist):
    """Run pagination tests using production shards."""
    print(f"\n=== Test 1: CURSOR PAGINATION (k={k_val}) ===")
    page_size = 5

    # --- 1a: First page ---
    print("  1a: First page (--limit 5)")
    stdout1, stderr1, rc1 = run_substr([
        "--shards", shard_dir, "--gc-hist", gc_hist,
        "--limit", str(page_size), "--random_access",
    ])
    if rc1 != 0:
        check(False, "first page: binary exited non-zero",
              f"rc={rc1} {stderr1[:200]}")
        return False

    meta1 = parse_meta(stdout1)
    if meta1 is None:
        check(False, "first page: __META__ line present")
        return False

    cursor1 = meta1["cursor"]
    has_more1 = meta1["hasMore"]
    returned1 = meta1["returned"]
    results1 = get_results(stdout1)

    check(cursor1 != "", "cursor is non-empty",
          f"got empty cursor")
    check(len(results1) > 0, f"first page has {len(results1)} results",
          "got 0 results")
    check(returned1 == len(results1),
          f"returned_count ({returned1}) matches actual row count ({len(results1)})",
          f"meta says {returned1}, got {len(results1)} rows")

    if len(results1) < page_size and has_more1:
        # If fewer results than page_size, hasMore should be falsy
        check(not has_more1 or has_more1 == "0" or has_more1 == "",
              "hasMore is false when results < page_size",
              f"hasMore={has_more1!r} with {len(results1)} results < limit={page_size}")
    elif len(results1) >= page_size and not has_more1:
        print(f"  INFO  hasMore is false even with {len(results1)} >= {page_size} "
              f"(end of data — valid)")

    if not results1:
        print("  [SKIP] No results on first page — cannot test pagination further.")
        return False

    # --- 1b: Second page with cursor ---
    print("  1b: Second page (--cursor from first page)")
    stdout2, stderr2, rc2 = run_substr([
        "--shards", shard_dir, "--gc-hist", gc_hist,
        "--limit", str(page_size), "--random_access",
        "--cursor", cursor1,
    ])
    if rc2 != 0:
        check(False, "second page: binary exited non-zero",
              f"rc={rc2} {stderr2[:200]}")
        return False

    meta2 = parse_meta(stdout2)
    results2 = get_results(stdout2)

    check(meta2 is not None, "second page: __META__ line present")

    if results2:
        # Verify no overlap between pages
        overlap = set(results1) & set(results2)
        check(len(overlap) == 0,
              "second page has no duplicate k-mers from first page",
              f"overlap={overlap}" if overlap else "")

        check(meta2 is None or meta2["returned"] == len(results2),
              "second page: returned_count matches row count",
              f"meta={meta2['returned'] if meta2 else '?'}, got {len(results2)}")
    else:
        print("  INFO  Second page returned no results (end of data reached)")

    # --- 1c: Collect all pages and verify no duplicates ---
    print("  1c: Collect all pages — verify no duplicates across all pages")
    all_kmers = list(results1)
    current_cursor = cursor1
    page_num = 2
    max_pages = 100  # safety limit

    while page_num <= max_pages:
        stdout, stderr, rc = run_substr([
            "--shards", shard_dir, "--gc-hist", gc_hist,
            "--limit", str(page_size), "--random_access",
            "--cursor", current_cursor,
        ])
        if rc != 0:
            check(False, f"page {page_num}: binary exited non-zero",
                  f"rc={rc} {stderr[:200]}")
            break

        meta = parse_meta(stdout)
        results = get_results(stdout)

        if not results:
            break

        # Check for duplicates
        new_kmers = [k for k in results if k not in all_kmers]
        if len(new_kmers) < len(results):
            dups = set(results) & set(all_kmers)
            check(False, f"page {page_num}: found {len(dups)} duplicate(s)",
                  str(dups))
            break

        all_kmers.extend(new_kmers)

        # Check metadata consistency
        if meta and meta["returned"] != len(results):
            check(False, f"page {page_num}: returned_count ({meta['returned']}) != "
                  f"actual rows ({len(results)})")

        # Get next cursor
        if meta and meta["cursor"]:
            current_cursor = meta["cursor"]
            # hasMore check
            has_more = meta.get("hasMore", "")
            if not has_more or has_more == "0":
                break
        else:
            break
        page_num += 1

    check(page_num <= max_pages, "pagination completed within safety limit",
          f"page_num={page_num}, max_pages={max_pages}")

    total_unique = len(set(all_kmers))
    total_collected = len(all_kmers)
    check(total_unique == total_collected,
          f"no duplicates across all {page_num} pages "
          f"({total_unique} unique, {total_collected} total)",
          f"unique={total_unique}, total={total_collected}")

    print(f"  INFO  Collected {total_collected} unique k-mers across "
          f"{page_num} page(s)")

    return True


# ===================================================================
#  Main
# ===================================================================
def main():
    global g_failures, g_passes
    g_failures = 0
    g_passes = 0

    if not os.path.exists(BINARY):
        print(f"ERROR: binary not found at {BINARY}")
        print("Run 'bash website/qc/build_query_binaries.sh' first.")
        sys.exit(1)

    k_val, shard_dir, gc_hist = pick_shards()
    if k_val is None:
        print("[SKIP] Production shards not found — skipping pagination tests.")
        print("To run full tests, ensure shards_16/ (or 17/18) and")
        print("corresponding gc_hist_shards_*.json are present.")
        print(f"Binary exists at: {BINARY}")
    else:
        print(f"Found shards: k={k_val} at {shard_dir}")
        test_cursor_pagination(k_val, shard_dir, gc_hist)

    total = g_passes + g_failures
    print(f"\n{'=' * 50}")
    if total > 0:
        print(f"Results: {g_passes} passed, {g_failures} failed out of {total} tests")
    else:
        print("No tests were run.")

    sys.exit(0)


if __name__ == "__main__":
    main()
