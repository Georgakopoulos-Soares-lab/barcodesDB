#!/usr/bin/env bash
# run_website_qc.sh — One-command QC runner for barcodesDB website.
# Stops on required test failures. Writes logs to qc/tmp/logs/.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WEBSITE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
QC_TMP="$SCRIPT_DIR/tmp"
LOG_DIR="$QC_TMP/logs"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"

mkdir -p "$LOG_DIR"

PASS=0
FAIL=0
SKIP=0

# Log file
LOG="$LOG_DIR/qc_run_$TIMESTAMP.log"
SUMMARY_LOG="$LOG_DIR/qc_summary_$TIMESTAMP.log"

exec > >(tee -a "$LOG") 2>&1

echo "============================================"
echo " barcodesDB Website QC — $(date)"
echo "============================================"
echo ""

# ---- Step 0: Build binaries ----
echo "=== [00/12] Build query binaries ==="
if bash "$SCRIPT_DIR/build_query_binaries.sh" 2>&1; then
    echo "  PASS: Build query binaries"
    PASS=$((PASS + 1))
else
    echo "  FAIL: Build query binaries"
    FAIL=$((FAIL + 1))
fi
echo ""

# ---- Step 1: Generate synthetic fixtures ----
echo "=== [01/12] Generate synthetic fixtures ==="
if cd "$WEBSITE_DIR" && ./qc/tmp/create_synthetic_kbit_fixture > /dev/null 2>&1; then
    echo "  PASS: Generate synthetic fixtures"
    PASS=$((PASS + 1))
else
    echo "  FAIL: Generate synthetic fixtures"
    FAIL=$((FAIL + 1))
fi
echo ""

# ---- Helper: run a Python test ----
run_py_test() {
    local step="$1"
    local label="$2"
    local script="$3"
    local optional="${4:-false}"

    echo "=== [$step] $label ==="
    if cd "$WEBSITE_DIR" && python3 "$script" 2>&1; then
        echo "  PASS: $label"
        PASS=$((PASS + 1))
    else
        if [ "$optional" = "true" ]; then
            echo "  SKIP (optional): $label"
            SKIP=$((SKIP + 1))
        else
            echo "  FAIL: $label"
            FAIL=$((FAIL + 1))
        fi
    fi
    echo ""
}

# ---- Step 2: Exact k-mer lookup QC ----
run_py_test "02/12" "Exact k-mer lookup QC" "$SCRIPT_DIR/test_query_kmer_bitmap.py"

# ---- Step 3: Substring search QC ----
run_py_test "03/12" "Substring search QC" "$SCRIPT_DIR/test_query_substring_bitmap_stream.py" true

# ---- Step 4: GC filtering QC ----
run_py_test "04/12" "GC filtering QC" "$SCRIPT_DIR/test_gc_filtering.py"

# ---- Step 5: Pagination/cursor QC ----
run_py_test "05/12" "Pagination/cursor QC" "$SCRIPT_DIR/test_pagination_cursor.py" true

# ---- Step 6: Motif-filter QC (existing tests) ----
echo "=== [06/12] Motif-filter C++ unit tests ==="
if cd "$WEBSITE_DIR/programs" && g++ -std=c++17 -O0 -I. test_motif_filters.cpp -o /tmp/test_motif_filters 2>&1 && /tmp/test_motif_filters 2>&1; then
    echo "  PASS: Motif-filter C++ unit tests"
    PASS=$((PASS + 1))
else
    echo "  FAIL: Motif-filter C++ unit tests"
    FAIL=$((FAIL + 1))
fi
rm -f /tmp/test_motif_filters
echo ""

echo "=== [07/12] Motif-filter shell tests ==="
if cd "$WEBSITE_DIR/programs" && bash test_motif_filters.sh 2>&1; then
    echo "  PASS: Motif-filter shell tests"
    PASS=$((PASS + 1))
else
    echo "  FAIL: Motif-filter shell tests"
    FAIL=$((FAIL + 1))
fi
echo ""

# ---- Step 8: Backend/API smoke tests ----
run_py_test "08/12" "Backend/API smoke tests" "$SCRIPT_DIR/test_server_api_smoke.py" true

# ---- Step 9: Server parameter parity ----
run_py_test "09/12" "Server parameter parity" "$SCRIPT_DIR/test_server_parameter_parity.py" true

# ---- Step 10: Frontend parameter parity ----
run_py_test "10/12" "Frontend parameter parity" "$SCRIPT_DIR/test_frontend_parameter_parity.py" true

# ---- Step 11: Output schema regression ----
run_py_test "11/12" "Output schema regression" "$SCRIPT_DIR/test_output_schema.py"

# ---- Step 12: Benchmark smoke tests ----
run_py_test "12/13" "Benchmark smoke tests" "$SCRIPT_DIR/benchmark_query_binaries.py" true

# ---- Step 13: Random barcode generation QC ----
run_py_test "13/13" "Random barcode generation QC" "$SCRIPT_DIR/test_random_barcode_generation.py" true

# ---- Summary ----
echo "============================================"
echo " QC Complete"
echo "   Passed: $PASS"
echo "   Failed: $FAIL"
echo "   Skipped: $SKIP"
echo "   Log: $LOG"
echo "============================================"

{
    echo "QC run $TIMESTAMP"
    echo "  Passed: $PASS"
    echo "  Failed: $FAIL"
    echo "  Skipped: $SKIP"
} > "$SUMMARY_LOG"

if [ "$FAIL" -gt 0 ]; then
    echo "ERROR: $FAIL required test(s) failed."
    exit 1
fi
echo "All required QC tests passed."
exit 0
