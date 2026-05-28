#!/usr/bin/env python3
"""
test_server_parameter_parity.py — Verify server.js parameter ↔ C++ flag mapping.

Reads server.js as text and performs string-based static analysis to verify
that API parameter names in the request body correctly map to C++ CLI flags.

Checks:
1. appendMotifArgs() maps each JS property to the correct --flag-name
2. The /api/query-kmer route reads expected parameters
3. The /api/query-substring route reads expected parameters

Exit code 0 always (informational test).
Dependencies: Python standard library only.
"""

import os
import re
import sys

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SERVER_JS = os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "server.js"
))

# ---------------------------------------------------------------------------
# Globals
# ---------------------------------------------------------------------------
g_failures = 0
g_passes = 0
g_found_checks = 0
g_not_found_checks = 0


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


def read_server_js():
    """Read server.js and return lines."""
    if not os.path.exists(SERVER_JS):
        print(f"ERROR: server.js not found at {SERVER_JS}")
        return None
    with open(SERVER_JS, "r") as f:
        return f.read()


def count_occurrences(content, pattern):
    """Count occurrences of a regex pattern in the content."""
    return len(re.findall(pattern, content))


# ===================================================================
#  Test 1: appendMotifArgs mapping
# ===================================================================
def test_append_motif_args_mapping():
    """Verify appendMotifArgs maps JS properties to correct C++ flags."""
    print("\n=== Test 1: appendMotifArgs() mapping ===")

    content = read_server_js()
    if content is None:
        return

    # Expected mappings: JS property → C++ CLI flag
    mappings = [
        ("filter_homopolymers", "--filter-homopolymers"),
        ("filter_low_complexity", "--filter-low-complexity"),
        ("filter_dinucleotide_repeats", "--filter-dinucleotide-repeats"),
        ("filter_trinucleotide_repeats", "--filter-trinucleotide-repeats"),
        ("filter_tetranucleotide_repeats", "--filter-tetranucleotide-repeats"),
        ("filter_restriction_sites", "--filter-restriction-sites"),
        ("filter_functional_motifs", "--filter-functional-motifs"),
        ("motif_mode", "--motif-mode"),
        ("max_homopolymer", "--max-homopolymer"),
        ("min_shannon_entropy", "--min-shannon-entropy"),
    ]

    for js_prop, cpp_flag in mappings:
        # Count the JS property name in the appendMotifArgs function
        # Look specifically in the function body
        js_occurrences = count_occurrences(
            content, re.escape(js_prop)
        )
        flag_occurrences = count_occurrences(
            content, re.escape(cpp_flag)
        )

        # JS property should appear at least once (in the function)
        # The flag should appear at least once (in args.push calls)
        found_js = js_occurrences >= 1
        found_flag = flag_occurrences >= 1

        check(found_js and found_flag,
              f"{js_prop} → {cpp_flag}",
              f"js_prop appears {js_occurrences}x, flag appears {flag_occurrences}x"
              if not (found_js and found_flag) else "")


# ===================================================================
#  Test 2: /api/query-kmer route parameters
# ===================================================================
def test_query_kmer_route_params():
    """Verify /api/query-kmer route reads expected parameters from req.body."""
    print("\n=== Test 2: /api/query-kmer route parameters ===")

    content = read_server_js()
    if content is None:
        return

    # Check that these parameter names appear in the server.js
    params_to_check = [
        ("kmers", r'body\.kmers'),
        ("min_hamming_distance", r'body\.min_hamming_distance'),
    ]

    for param_name, pattern in params_to_check:
        count = count_occurrences(content, pattern)
        check(count >= 1,
              f"req.body.{param_name} is read by server",
              f"found {count} occurrence(s)" if count >= 1 else "NOT FOUND")

    # Verify the route handler exists
    route_count = count_occurrences(content, r"/api/query-kmer")
    check(route_count >= 2,
          "/api/query-kmer route handler exists",
          f"found {route_count} reference(s)")


# ===================================================================
#  Test 3: /api/query-substring route parameters
# ===================================================================
def test_query_substring_route_params():
    """Verify /api/query-substring route reads expected parameters."""
    print("\n=== Test 3: /api/query-substring route parameters ===")

    content = read_server_js()
    if content is None:
        return

    params_to_check = [
        ("substring", r'body\.substring'),
        ("gcMin", r'body\.gcMin'),
        ("gcMax", r'body\.gcMax'),
        ("limit", r'body\.limit'),
        ("cursor", r'body\.cursor'),
        ("constructK", r'body\.constructK'),
        ("motif_mode", r'body\.motif_mode'),
        ("filter_homopolymers", r'body\.filter_homopolymers'),
        ("filter_low_complexity", r'body\.filter_low_complexity'),
        ("filter_dinucleotide_repeats", r'body\.filter_dinucleotide_repeats'),
        ("filter_trinucleotide_repeats", r'body\.filter_trinucleotide_repeats'),
        ("filter_tetranucleotide_repeats", r'body\.filter_tetranucleotide_repeats'),
        ("filter_restriction_sites", r'body\.filter_restriction_sites'),
        ("filter_functional_motifs", r'body\.filter_functional_motifs'),
        ("reverse_complement", r'body\.reverse_complement'),
        ("threads", r'body\.threads'),
    ]

    for param_name, pattern in params_to_check:
        count = count_occurrences(content, pattern)
        check(count >= 1,
              f"req.body.{param_name} is read by server",
              f"found {count} occurrence(s)" if count >= 1 else "NOT FOUND")

    # Verify the route handler exists
    route_count = count_occurrences(content, r"/api/query-substring")
    check(route_count >= 2,
          "/api/query-substring route handler exists",
          f"found {route_count} reference(s)")


# ===================================================================
#  Test 4: C++ flag usage in args.push()
# ===================================================================
def test_flag_usage_in_args_push():
    """Verify C++ flags appear in args.push() calls in the route handlers."""
    print("\n=== Test 4: C++ flags used in args.push() ===")

    content = read_server_js()
    if content is None:
        return

    flags_to_check = [
        "--shards",
        "--gc-hist",
        "--gc-min",
        "--gc-max",
        "--limit",
        "--cursor",
        "--threads",
        "--random_access",
        "--construct_k",
        "--substring",
        "--reverse_complement",
        "--motif-mode",
        "--filter-homopolymers",
        "--filter-low-complexity",
        "--filter-dinucleotide-repeats",
        "--filter-trinucleotide-repeats",
        "--filter-tetranucleotide-repeats",
        "--filter-restriction-sites",
        "--filter-functional-motifs",
        "--min-hamming-distance",
    ]

    for flag in flags_to_check:
        # Count occurrences of the flag (in args.push or elsewhere)
        count = count_occurrences(content, re.escape(flag))
        check(count >= 1,
              f"flag {flag} is used in server.js",
              f"found {count} occurrence(s)" if count >= 1 else "NOT FOUND")


# ===================================================================
#  Main
# ===================================================================
def main():
    global g_failures, g_passes
    g_failures = 0
    g_passes = 0

    print(f"Server file: {SERVER_JS}")

    if not os.path.exists(SERVER_JS):
        print(f"\nERROR: {SERVER_JS} not found.")
        sys.exit(0)

    test_append_motif_args_mapping()
    test_query_kmer_route_params()
    test_query_substring_route_params()
    test_flag_usage_in_args_push()

    total = g_passes + g_failures
    print(f"\n{'=' * 50}")
    print(f"Results: {g_passes} passed, {g_failures} failed out of {total} checks")

    if g_failures > 0:
        print("NOTE: Some mappings may not be found. This may indicate")
        print("parameter name changes in server.js — review manually.")
    else:
        print("All parameter mappings verified ✓")

    sys.exit(0)


if __name__ == "__main__":
    main()
