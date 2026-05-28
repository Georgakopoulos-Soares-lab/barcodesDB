#!/usr/bin/env python3
"""
Smoke tests for the barcodesDB server API (server.js).

Usage:
    # With server already running on the default port (8090 or $PORT):
    python3 test_server_api_smoke.py

    # Specify an explicit port:
    PORT=8090 python3 test_server_api_smoke.py

Exit code: 0 if all tests pass, 1 if any fail.
"""

import json
import os
import sys
import time
import urllib.error
import urllib.request

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
PORT = int(os.environ.get("PORT", "8090"))
BASE = f"http://localhost:{PORT}"
TIMEOUT = 30  # seconds for each HTTP request

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
results_log = []
fail_count = 0
pass_count = 0


def request(method, path, body=None):
    """Make an HTTP request and return (status_code, data_or_None, error_msg_or_None)."""
    url = f"{BASE}{path}"
    data = None
    if body is not None:
        data = json.dumps(body).encode("utf-8")

    req = urllib.request.Request(url, data=data, method=method)
    req.add_header("Content-Type", "application/json")

    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
            raw = resp.read()
            status = resp.status
            try:
                decoded = json.loads(raw)
            except (json.JSONDecodeError, ValueError):
                decoded = raw.decode("utf-8", errors="replace")
            return status, decoded, None
    except urllib.error.HTTPError as e:
        status = e.code
        try:
            decoded = json.loads(e.read())
        except (json.JSONDecodeError, ValueError, OSError):
            decoded = None
        return status, decoded, None
    except urllib.error.URLError as e:
        return None, None, f"Connection failed: {e.reason}"
    except OSError as e:
        return None, None, str(e)


def test(name, ok, detail=""):
    """Log a test result."""
    global pass_count, fail_count
    if ok:
        pass_count += 1
        status = "PASS"
    else:
        fail_count += 1
        status = "FAIL"
    line = f"[{status}] {name}"
    if detail:
        line += f" - {detail}"
    results_log.append(line)
    print(line)


def check(condition, msg=""):
    """Return True if condition is truthy, else False + record reason."""
    if condition:
        return True, ""
    return False, msg


# ---------------------------------------------------------------------------
# Server availability check
# ---------------------------------------------------------------------------
def server_is_running():
    """Try connecting to the server; return True if it responds."""
    try:
        req = urllib.request.Request(f"{BASE}/")
        with urllib.request.urlopen(req, timeout=5) as resp:
            return resp.status == 200
    except (urllib.error.URLError, OSError):
        return False


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
def test_01_post_query_kmer_present():
    """Known-present k-mer -> 200, present=true, full response shape."""
    status, data, err = request("POST", "/api/query-kmer",
                                {"kmers": ["ACGTACGTACGTACGT"]})
    if err:
        return test("T1 POST /api/query-kmer (present)", False, err)
    ok, reason = check(
        status == 200
        and data is not None
        and isinstance(data, dict)
        and "results" in data
        and isinstance(data["results"], list)
        and len(data["results"]) == 1
        and data["results"][0].get("present") is True
        and "total" in data
        and "found" in data
        and "foundPct" in data,
        f"status={status}, data={json.dumps(data)[:200]}"
    )
    test("T1 POST /api/query-kmer (present)", ok, reason)


def test_02_post_query_kmer_absent():
    """Known-absent k-mer -> 200, present=false."""
    status, data, err = request("POST", "/api/query-kmer",
                                {"kmers": ["AAAAAAAAAAAAAAAA"]})
    if err:
        return test("T2 POST /api/query-kmer (absent)", False, err)
    ok, reason = check(
        status == 200
        and data is not None
        and isinstance(data, dict)
        and "results" in data
        and data["results"][0].get("present") is False,
        f"status={status}, data={json.dumps(data)[:200]}"
    )
    test("T2 POST /api/query-kmer (absent)", ok, reason)


def test_03_post_query_kmer_invalid():
    """Invalid k-mer (contains N) -> 400."""
    status, data, err = request("POST", "/api/query-kmer",
                                {"kmers": ["ACGTACGTACGTACGN"]})
    if err:
        return test("T3 POST /api/query-kmer (invalid char)", False, err)
    ok, reason = check(
        status == 400
        and data is not None
        and "error" in data
        and "invalid" in str(data.get("error", "")).lower(),
        f"status={status}, data={json.dumps(data)[:200]}"
    )
    test("T3 POST /api/query-kmer (invalid char)", ok, reason)


def test_04_post_query_kmer_mixed_lengths():
    """Mixed k-mer lengths -> 400."""
    status, data, err = request("POST", "/api/query-kmer",
                                {"kmers": ["ACGTACGTACGTACGT", "ACGTACGTACGTAC"]})
    if err:
        return test("T4 POST /api/query-kmer (mixed lengths)", False, err)
    ok, reason = check(
        status == 400
        and data is not None
        and "error" in data,
        f"status={status}, data={json.dumps(data)[:200]}"
    )
    test("T4 POST /api/query-kmer (mixed lengths)", ok, reason)


def test_05_post_query_kmer_min_hamming():
    """min_hamming_distance parameter -> 200, response includes hamming fields."""
    status, data, err = request("POST", "/api/query-kmer", {
        "kmers": ["ACGTACGTACGTACGT", "ACGTACGTACGTACGG"],
        "min_hamming_distance": 2,
    })
    if err:
        return test("T5 POST /api/query-kmer (min_hamming_distance)", False, err)
    ok, reason = check(
        status == 200
        and data is not None
        and isinstance(data, dict)
        and "min_hamming_distance_requested" in data
        and data.get("min_hamming_distance_requested") == 2
        and "hamming_distance_check_applied" in data
        and data.get("hamming_distance_check_applied") is True,
        f"status={status}, data={json.dumps(data)[:300]}"
    )
    test("T5 POST /api/query-kmer (min_hamming_distance)", ok, reason)


def test_06_post_query_kmer_empty_list():
    """Empty kmers list -> 400."""
    status, data, err = request("POST", "/api/query-kmer", {"kmers": []})
    if err:
        return test("T6 POST /api/query-kmer (empty list)", False, err)
    ok, reason = check(
        status == 400,
        f"status={status}, data={json.dumps(data)[:200]}"
    )
    test("T6 POST /api/query-kmer (empty list)", ok, reason)


def test_07_post_query_substring():
    """POST /api/query-substring with basic params -> 200, results <= limit."""
    status, data, err = request("POST", "/api/query-substring",
                                {"gcMin": 0, "gcMax": 100, "limit": 5})
    if err:
        return test("T7 POST /api/query-substring", False, err)
    ok, reason = check(
        status == 200
        and data is not None
        and isinstance(data, dict)
        and "results" in data
        and isinstance(data["results"], list)
        and len(data["results"]) <= 5,
        f"status={status}, results_count={len(data.get('results', [])) if isinstance(data, dict) else 'N/A'}"
    )
    test("T7 POST /api/query-substring", ok, reason)


def test_08_post_query_substring_motif():
    """POST /api/query-substring with motif params -> 200."""
    status, data, err = request("POST", "/api/query-substring", {
        "gcMin": 0, "gcMax": 100, "limit": 5,
        "motif_mode": "flag", "filter_homopolymers": True,
    })
    if err:
        return test("T8 POST /api/query-substring (motif)", False, err)
    ok, reason = check(
        status == 200,
        f"status={status}, data={json.dumps(data)[:300]}"
    )
    test("T8 POST /api/query-substring (motif)", ok, reason)


def test_09_get_kmer_page():
    """GET /kmer -> 200, body contains 'K-mer'."""
    status, data, err = request("GET", "/kmer")
    if err:
        return test("T9 GET /kmer page", False, err)
    body = data if isinstance(data, str) else str(data)
    ok, reason = check(
        status == 200
        and "K-mer" in body,
        f"status={status}, body_contains_K-mer={'K-mer' in body}"
    )
    test("T9 GET /kmer page", ok, reason)


def test_10_get_substring_page():
    """GET /substring -> 200, body contains 'Search barcodes'."""
    status, data, err = request("GET", "/substring")
    if err:
        return test("T10 GET /substring page", False, err)
    body = data if isinstance(data, str) else str(data)
    ok, reason = check(
        status == 200
        and "Search barcodes" in body,
        f"status={status}, body_contains_Search_barcodes={'Search barcodes' in body}"
    )
    test("T10 GET /substring page", ok, reason)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print(f"=== barcodesDB API Smoke Tests (port {PORT}) ===")
    print()

    if not server_is_running():
        print(f"Server not reachable at {BASE}.")
        print(f"Start it with:")
        print(f"    cd website && PORT={PORT} node server.js")
        print()
        print("=== ALL TESTS SKIPPED ===")
        sys.exit(0)

    print(f"Server reachable at {BASE}. Running tests...\n")

    test_01_post_query_kmer_present()
    test_02_post_query_kmer_absent()
    test_03_post_query_kmer_invalid()
    test_04_post_query_kmer_mixed_lengths()
    test_05_post_query_kmer_min_hamming()
    test_06_post_query_kmer_empty_list()
    test_07_post_query_substring()
    test_08_post_query_substring_motif()
    test_09_get_kmer_page()
    test_10_get_substring_page()

    print()
    print(f"=== Results: {pass_count} passed, {fail_count} failed ===")

    if fail_count > 0:
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    main()
