#!/usr/bin/env bash
# build_query_binaries.sh
# Compile website binaries and QC tools into qc/tmp/.
# Does NOT overwrite production binaries in the parent directory.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WEBSITE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
ROOT_DIR="$(cd "$WEBSITE_DIR/.." && pwd)"
QC_TMP="$WEBSITE_DIR/qc/tmp"
PROGRAMS="$WEBSITE_DIR/programs"

# Roaring installation
ROARING_INC="/home1/10899/kimopro/WORK/tools/include"
ROARING_LIB="/home1/10899/kimopro/WORK/tools/lib64/libroaring.a"

if [ ! -f "$ROARING_LIB" ]; then
    echo "ERROR: Roaring library not found at $ROARING_LIB"
    echo "Run ../install_roaring.sh first."
    exit 1
fi

if [ ! -d "$ROARING_INC/roaring" ]; then
    echo "ERROR: Roaring headers not found at $ROARING_INC/roaring"
    exit 1
fi

mkdir -p "$QC_TMP"
CXX="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:--O3 -march=native -flto -std=c++17 -pthread}"
LINKFLAGS="-lm"

echo "=== Building query_kmer_bitmap (QC) ==="
$CXX $CXXFLAGS \
    -I"$ROARING_INC" \
    "$PROGRAMS/query_kmer_bitmap.cpp" \
    "$ROARING_LIB" \
    $LINKFLAGS \
    -o "$QC_TMP/query_kmer_bitmap"
echo "  -> $QC_TMP/query_kmer_bitmap"

echo ""
echo "=== Building query_substring_bitmap_stream (QC) ==="
$CXX $CXXFLAGS \
    -I"$ROARING_INC" \
    "$PROGRAMS/query_substring_bitmap_stream.cpp" \
    "$ROARING_LIB" \
    $LINKFLAGS \
    -o "$QC_TMP/query_substring_bitmap_stream"
echo "  -> $QC_TMP/query_substring_bitmap_stream"

echo ""
echo "=== Building synthetic fixture generator ==="
$CXX -O3 -std=c++17 \
    -I"$ROARING_INC" \
    "$SCRIPT_DIR/create_synthetic_kbit_fixture.cpp" \
    "$ROARING_LIB" \
    -pthread -lm \
    -o "$QC_TMP/create_synthetic_kbit_fixture"
echo "  -> $QC_TMP/create_synthetic_kbit_fixture"

echo ""
echo "=== Build complete ==="
ls -lh "$QC_TMP"/query_kmer_bitmap "$QC_TMP"/query_substring_bitmap_stream "$QC_TMP"/create_synthetic_kbit_fixture 2>/dev/null
