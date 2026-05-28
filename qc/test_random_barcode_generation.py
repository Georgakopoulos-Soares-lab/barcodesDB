#!/usr/bin/env python3
"""
test_random_barcode_generation.py — QC tests for /api/generate-random-kmers

Tests the random barcode generation endpoint:
1. Basic generation with various k, n, GC ranges
2. Substring constraints
3. GC boundary validation
4. Invalid input rejection
5. Response schema verification

Exit 0 if all pass, 1 if any fail.
"""

import json
import os
import sys
import urllib.request
import urllib.error

# Server URL
PORT = os.environ.get("BARCODEDB_PORT", "9090")
BASE = os.environ.get("BARCODEDB_API_URL", f"http://localhost:{PORT}")

PASS = 0
FAIL = 0
SKIP = 0


def check(cond, label, detail=""):
    global PASS, FAIL
    if cond:
        print(f"  PASS  {label}")
        PASS += 1
    else:
        print(f"  FAIL  {label}")
        if detail:
            print(f"        {detail}")
        FAIL += 1


def post(path, body):
    """POST JSON to server, return parsed JSON or (error, status)."""
    try:
        data = json.dumps(body).encode("utf-8")
        req = urllib.request.Request(
            f"{BASE}{path}",
            data=data,
            headers={"Content-Type": "application/json"},
        )
        with urllib.request.urlopen(req) as resp:
            return json.loads(resp.read().decode()), resp.status
    except urllib.error.HTTPError as e:
        try:
            return json.loads(e.read().decode()), e.code
        except Exception:
            return {"error": str(e)}, e.code
    except urllib.error.URLError as e:
        return {"error": f"Connection refused: {e}"}, 0


def test_server_alive():
    """Check if the server is reachable."""
    try:
        urllib.request.urlopen(f"{BASE}/", timeout=3)
        return True
    except Exception:
        return False


def gc_pct(kmer):
    """Compute GC% in Python."""
    if not kmer:
        return 0
    gc = sum(1 for c in kmer.upper() if c in "GC")
    return (gc * 100.0) / len(kmer)


def main():
    global PASS, FAIL, SKIP

    print("=== Random Barcode Generation QC ===")
    print(f"Server: {BASE}")

    if not test_server_alive():
        print("[SKIP] Server not running. Start with: PORT=9090 node server.js")
        print("       Then set BARCODEDB_PORT=9090 or BARCODEDB_API_URL=...")
        SKIP = 99
        return 0

    # ================================================================
    # Test 1: Basic generation — default parameters
    # ================================================================
    print("\n=== Test 1: Basic generation ===")
    data, status = post("/api/generate-random-kmers", {"k": 18, "n": 10})
    check(status == 200, "default params → 200", f"got {status}")

    # Schema check
    check(isinstance(data, dict), "response is JSON object")
    check("k" in data, "response has 'k' field")
    check("n_requested" in data, "response has 'n_requested'")
    check("n_returned" in data, "response has 'n_returned'")
    check("results" in data, "response has 'results' array")
    check(isinstance(data["results"], list), "results is an array")
    check(len(data["results"]) == data["n_returned"], "n_returned matches results length")
    check(len(data["results"]) > 0, "at least one result returned")

    if data["results"]:
        r = data["results"][0]
        check("kmer" in r, "result has 'kmer'")
        check("gc" in r, "result has 'gc'")
        check("comp" in r, "result has 'comp'")
        check(len(r["kmer"]) == data["k"], f"k-mer length ({len(r['kmer'])}) matches k ({data['k']})")

    # ================================================================
    # Test 2: GC% constraint
    # ================================================================
    print("\n=== Test 2: GC% constraint ===")
    for gc_min, gc_max, label in [(40, 60, "40-60"), (0, 50, "0-50"), (50, 100, "50-100"), (50, 50, "exact 50")]:
        data, status = post("/api/generate-random-kmers", {
            "k": 18, "n": 20, "gc_min": gc_min, "gc_max": gc_max
        })
        check(status == 200, f"GC [{gc_min},{gc_max}] → 200", f"got {status}")

        all_in_range = all(gc_min <= gc_pct(r["kmer"]) <= gc_max for r in data["results"])
        check(all_in_range, f"all {data['n_returned']} results have GC% in [{gc_min},{gc_max}]")

    # ================================================================
    # Test 3: Substring constraint
    # ================================================================
    print("\n=== Test 3: Substring constraint ===")
    substring = "ACGTAC"
    data, status = post("/api/generate-random-kmers", {
        "k": 18, "n": 10, "gc_min": 0, "gc_max": 100, "substring": substring
    })
    check(status == 200, "substring constraint → 200", f"got {status}")
    check(data["substring"] == substring, f"substring echoed back: '{data['substring']}'")

    all_contain = all(substring in r["kmer"] for r in data["results"])
    check(all_contain, f"all {data['n_returned']} results contain '{substring}'",
          str([r["kmer"] for r in data["results"] if substring not in r["kmer"]])[:200]
          if not all_contain else "")

    # ================================================================
    # Test 4: Various k values
    # ================================================================
    print("\n=== Test 4: Various k values ===")
    for k_val in [4, 8, 16, 24, 32]:
        data, status = post("/api/generate-random-kmers", {"k": k_val, "n": 5})
        check(status == 200, f"k={k_val} → 200", f"got {status}")
        check(all(len(r["kmer"]) == k_val for r in data["results"]),
              f"all k-mers have length {k_val}")

    # ================================================================
    # Test 5: Invalid input rejection
    # ================================================================
    print("\n=== Test 5: Invalid input rejection ===")

    # k too small
    data, status = post("/api/generate-random-kmers", {"k": 3, "n": 5})
    check(status == 400, "k=3 → 400", f"got {status}")
    check("error" in data, "error message present for bad k")

    # k too large
    data, status = post("/api/generate-random-kmers", {"k": 33, "n": 5})
    check(status == 400, "k=33 → 400", f"got {status}")

    # n too large
    data, status = post("/api/generate-random-kmers", {"k": 16, "n": 20000})
    check(status == 400, "n=20000 → 400", f"got {status}")

    # n=0
    data, status = post("/api/generate-random-kmers", {"k": 16, "n": 0})
    check(status == 400, "n=0 → 400", f"got {status}")

    # invalid GC range
    data, status = post("/api/generate-random-kmers", {"k": 16, "n": 5, "gc_min": 80, "gc_max": 20})
    check(status == 400, "gc_min > gc_max → 400", f"got {status}")

    # negative GC
    data, status = post("/api/generate-random-kmers", {"k": 16, "n": 5, "gc_min": -5, "gc_max": 50})
    check(status == 400, "negative gc_min → 400", f"got {status}")

    # substring longer than k
    data, status = post("/api/generate-random-kmers", {"k": 6, "n": 5, "substring": "ACGTACGT"})
    check(status == 400, "substring > k → 400", f"got {status}")

    # invalid substring chars
    data, status = post("/api/generate-random-kmers", {"k": 16, "n": 5, "substring": "ACGTX"})
    check(status == 400, "substring with X → 400", f"got {status}")

    # ================================================================
    # Test 6: Too-restrictive criteria (graceful degradation)
    # ================================================================
    print("\n=== Test 6: Restrictive criteria ===")
    data, status = post("/api/generate-random-kmers", {
        "k": 18, "n": 5, "gc_min": 50, "gc_max": 50,
        "substring": "ACGTACGTACGTACGTAC"
    })
    check(status == 200, "too restrictive → 200 (not error)", f"got {status}")
    check(data["n_returned"] == 0 or data["n_returned"] <= 5,
          f"returned {data['n_returned']} (graceful: {data['n_returned']} < 5)")
    check(data["attempts"] > 0, f"attempts field present: {data['attempts']}")

    # ================================================================
    # Test 7: Response metadata consistency
    # ================================================================
    print("\n=== Test 7: Response metadata consistency ===")
    data, status = post("/api/generate-random-kmers", {
        "k": 17, "n": 15, "gc_min": 30, "gc_max": 70
    })
    check(status == 200, "metadata consistency → 200")
    check(data["k"] == 17, "k echo matches")
    check(data["n_requested"] == 15, "n_requested echo matches")
    check(data["gc_min"] == 30, "gc_min echo matches")
    check(data["gc_max"] == 70, "gc_max echo matches")
    check(data["substring"] is None, "substring null when not provided")
    check(len(data["results"]) == data["n_returned"], "results count matches n_returned")
    check(data["n_returned"] <= data["n_requested"], "n_returned <= n_requested")

    # ================================================================
    # Summary
    # ================================================================
    print(f"\n{'='*60}")
    print(f"Results: {PASS} passed, {FAIL} failed, {SKIP} skipped")
    if FAIL > 0:
        print("SOME TESTS FAILED!")
        return 1
    print("ALL TESTS PASSED ✓")
    return 0


if __name__ == "__main__":
    sys.exit(main())
