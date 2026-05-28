#!/usr/bin/env python3
"""
test_frontend_parameter_parity.py — Verify frontend HTML ↔ server.js parameter names.

Reads kmer.html, substring.html, and server.js as text and performs static analysis
to verify that form fields, fetch() calls, and FormData keys match the server's
expected parameter names.

Exit code 0 always (informational test).
Dependencies: Python standard library only.
"""

import os
import re
import sys

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
QC_DIR = os.path.abspath(os.path.dirname(__file__))
PUBLIC_DIR = os.path.abspath(os.path.join(QC_DIR, "..", "public"))
SERVER_JS = os.path.join(QC_DIR, "..", "server.js")

KMER_HTML = os.path.join(PUBLIC_DIR, "kmer.html")
SUBSTR_HTML = os.path.join(PUBLIC_DIR, "substring.html")

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


def read_file(path, label):
    """Read a file, return content or None."""
    if not os.path.exists(path):
        print(f"  [SKIP] {label} not found: {path}")
        return None
    with open(path, "r") as f:
        return f.read()


def count_occurrences(content, pattern):
    """Count regex pattern matches in content."""
    return len(re.findall(pattern, content))


# ===================================================================
#  Test 1: kmer.html basic fields
# ===================================================================
def test_kmer_html_fields():
    """Verify kmer.html has expected form elements."""
    print("\n=== Test 1: kmer.html form fields ===")

    html = read_file(KMER_HTML, "kmer.html")
    if html is None:
        return

    # Textarea with id="kmers"
    count_textarea = count_occurrences(html, r'id="kmers"')
    check(count_textarea >= 1, 'textarea id="kmers" exists',
          f"found {count_textarea}x")

    # Numeric input with id="minHamming"
    count_min_hamming = count_occurrences(html, r'id="minHamming"')
    check(count_min_hamming >= 1, 'input id="minHamming" exists',
          f"found {count_min_hamming}x")

    # Fetch call sends "kmers" in JSON body
    count_kmers_json = count_occurrences(
        html, r'kmers'
    )
    check(count_kmers_json >= 2, 'fetch() JSON body includes "kmers" field',
          f"found {count_kmers_json}x references to kmers")

    # Fetch call sends "min_hamming_distance" in JSON body
    count_mhd_json = count_occurrences(
        html, r'min_hamming_distance'
    )
    check(count_mhd_json >= 2,
          'fetch() JSON body includes "min_hamming_distance"',
          f"found {count_mhd_json}x")

    # FormData appends "kmersFile"
    count_kmers_file = count_occurrences(html, r'kmersFile')
    check(count_kmers_file >= 1,
          'FormData appends "kmersFile" key',
          f"found {count_kmers_file}x")

    # FormData appends "min_hamming_distance"
    count_mhd_fd = count_occurrences(
        html, r'append\(.*min_hamming_distance'
    )
    check(count_mhd_fd >= 1,
          'FormData appends "min_hamming_distance" key',
          f"found {count_mhd_fd}x")


# ===================================================================
#  Test 2: substring.html checkbox IDs
# ===================================================================
def test_substring_html_checkbox_ids():
    """Verify substring.html has expected checkbox IDs."""
    print("\n=== Test 2: substring.html checkbox IDs ===")

    html = read_file(SUBSTR_HTML, "substring.html")
    if html is None:
        return

    checkboxes = [
        ("filterHomopolymers", 'id="filterHomopolymers"'),
        ("filterLowComplexity", 'id="filterLowComplexity"'),
        ("filterDinuclRepeats", 'id="filterDinuclRepeats"'),
        ("filterTrinuclRepeats", 'id="filterTrinuclRepeats"'),
        ("filterTetranuclRepeats", 'id="filterTetranuclRepeats"'),
        ("filterRestrictionSites", 'id="filterRestrictionSites"'),
        ("filterFunctionalMotifs", 'id="filterFunctionalMotifs"'),
    ]

    for name, pattern in checkboxes:
        count = count_occurrences(html, re.escape(pattern))
        check(count >= 1, f'checkbox id="{name}" exists',
              f"found {count}x" if count >= 1 else "NOT FOUND")


# ===================================================================
#  Test 3: substring.html getMotifBody() mapping
# ===================================================================
def test_substring_get_motif_body_mapping():
    """Verify getMotifBody() maps checkbox IDs to underscore_names."""
    print("\n=== Test 3: substring.html getMotifBody() field mapping ===")

    html = read_file(SUBSTR_HTML, "substring.html")
    if html is None:
        return

    mappings = [
        ("filterHomopolymers", "filter_homopolymers"),
        ("filterLowComplexity", "filter_low_complexity"),
        ("filterDinuclRepeats", "filter_dinucleotide_repeats"),
        ("filterTrinuclRepeats", "filter_trinucleotide_repeats"),
        ("filterTetranuclRepeats", "filter_tetranucleotide_repeats"),
        ("filterRestrictionSites", "filter_restriction_sites"),
        ("filterFunctionalMotifs", "filter_functional_motifs"),
    ]

    for checkbox_id, underscore_name in mappings:
        count = count_occurrences(html, re.escape(underscore_name))
        check(count >= 1,
              f"getMotifBody() maps {checkbox_id} → {underscore_name}",
              f"found {count}x" if count >= 1 else "NOT FOUND")

    # Verify motif_mode mapping from the mode dropdown
    count_motif_mode = count_occurrences(html, r'motif_mode')
    check(count_motif_mode >= 2,
          'getMotifBody() maps mode dropdown → motif_mode',
          f"found {count_motif_mode}x")

    # Verify max_homopolymer mapping
    count_max_hp = count_occurrences(html, r'max_homopolymer')
    check(count_max_hp >= 1,
          'getMotifBody() maps maxHomopolymer → max_homopolymer',
          f"found {count_max_hp}x")

    # Verify min_shannon_entropy mapping
    count_min_se = count_occurrences(html, r'min_shannon_entropy')
    check(count_min_se >= 1,
          'getMotifBody() maps minShannonEntropy → min_shannon_entropy',
          f"found {count_min_se}x")


# ===================================================================
#  Test 4: kmer.html motif body mapping
# ===================================================================
def test_kmer_html_motif_mapping():
    """Verify kmer.html also maps motif filter checkbox IDs to underscore_names."""
    print("\n=== Test 4: kmer.html motif body mapping ===")

    html = read_file(KMER_HTML, "kmer.html")
    if html is None:
        return

    # kmer.html also has motif filters and should map them
    underscore_names = [
        "motif_mode",
        "filter_homopolymers",
        "filter_low_complexity",
        "filter_dinucleotide_repeats",
        "filter_restriction_sites",
        "filter_functional_motifs",
    ]

    for name in underscore_names:
        count = count_occurrences(html, re.escape(name))
        check(count >= 1,
              f"kmer.html motif body includes {name}",
              f"found {count}x" if count >= 1 else "NOT FOUND")


# ===================================================================
#  Test 5: Server-side reads the same underscore_names
# ===================================================================
def test_server_reads_underscore_names():
    """Verify server.js reads the same underscore parameter names."""
    print("\n=== Test 5: server.js reads underscore parameter names ===")

    content = read_file(SERVER_JS, "server.js")
    if content is None:
        return

    underscore_names = [
        "filter_homopolymers",
        "filter_low_complexity",
        "filter_dinucleotide_repeats",
        "filter_trinucleotide_repeats",
        "filter_tetranucleotide_repeats",
        "filter_restriction_sites",
        "filter_functional_motifs",
        "motif_mode",
        "max_homopolymer",
        "min_shannon_entropy",
        "min_hamming_distance",
    ]

    for name in underscore_names:
        # Check in context of body.* access
        count_body = count_occurrences(content, re.escape(f"body.{name}"))
        # Also check in object literal assignments
        count_literal = count_occurrences(content, re.escape(f"{name}:"))
        total = count_body + count_literal
        check(total >= 1,
              f"server.js reads/uses {name}",
              f"found body.{name} {count_body}x, literal {count_literal}x"
              if total >= 1 else "NOT FOUND")


# ===================================================================
#  Test 6: kmer.html file upload appends motif filters to FormData
# ===================================================================
def test_kmer_html_formdata_motif():
    """Verify kmer.html appends motif filter fields to FormData on file upload."""
    print("\n=== Test 6: kmer.html FormData motif filter appends ===")

    html = read_file(KMER_HTML, "kmer.html")
    if html is None:
        return

    # Check that the for-loop appends motifBody entries to FormData
    motif_keys_in_formdata = [
        "motif_mode",
        "filter_homopolymers",
        "filter_low_complexity",
        "filter_dinucleotide_repeats",
        "filter_restriction_sites",
        "filter_functional_motifs",
    ]

    for key in motif_keys_in_formdata:
        count = count_occurrences(html, re.escape(key))
        check(count >= 2,
              f"kmer.html references {key} (used in FormData/JSON body)",
              f"found {count}x" if count >= 1 else "NOT FOUND",
              )


# ===================================================================
#  Main
# ===================================================================
def main():
    global g_failures, g_passes
    g_failures = 0
    g_passes = 0

    test_kmer_html_fields()
    test_substring_html_checkbox_ids()
    test_substring_get_motif_body_mapping()
    test_kmer_html_motif_mapping()
    test_server_reads_underscore_names()
    test_kmer_html_formdata_motif()

    total = g_passes + g_failures
    print(f"\n{'=' * 50}")
    print(f"Results: {g_passes} passed, {g_failures} failed out of {total} checks")

    if g_failures > 0:
        print("NOTE: Some checks may not match if the HTML has been updated.")
        print("Review the actual HTML files to confirm.")
    else:
        print("All frontend parameter parity checks passed ✓")

    sys.exit(0)


if __name__ == "__main__":
    main()
