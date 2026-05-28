# barcodesDB

[![barcodesDB](https://img.shields.io/badge/barcodesDB-DNA%20Barcode%20Database-blue)](https://barcodesdb.com)

barcodesDB is a web application designed for researchers and bioinformaticians working with DNA barcodes. It provides high-performance access to k-mer existence bitmaps and barcode search capabilities through optimized C++ tools and a user-friendly web interface.

## Features

### K-mer Existence Lookup
- Verify the existence of k-mers in the barcode database
- Support for batch processing via file uploads
- Detailed statistics including composition and GC content
- Flexible input options: direct text entry or file upload

### Barcode Search
- Filter barcodes by substring and GC content ranges
- Real-time streaming of results with pagination
- Export filtered results for further analysis

### Performance Optimizations
- Sharded bitmap architecture for handling large datasets
- Absent-mode queries for efficient gap analysis
- Multithreaded processing for improved speed
- Memory-efficient Roaring bitmap compression

### Sequence Motif Filters (Optional)
- Flag or exclude barcodes containing short motifs: homopolymer runs, low-complexity sequences, dinucleotide repeats, common restriction enzyme sites, and selected functional-like motifs
- All filters are OFF by default — existing behavior is unchanged
- Operate directly on candidate barcode strings; not organism-specific functional annotation
- Three modes: **Off** (no-op), **Flag only** (report but do not exclude), **Exclude** (filter out sequences that fail)
- Integrated into both C++ binaries and accessible via the web UI and API

## Technical Overview

- **Backend**: C++17 applications utilizing Roaring bitmap libraries
- **Frontend**: Responsive web interface using vanilla JavaScript
- **APIs**: RESTful endpoints for programmatic integration
- **Data Format**: KBITv1 format supporting dense and compressed bitmaps
- **Scalability**: Supports k-mer lengths up to 31

## Core Programs

barcodesDB relies on two primary C++ programs:

- **query_kmer_bitmap**: Performs batch k-mer existence queries against sharded bitmaps, with multithreaded processing and efficient lookups
- **query_substring_bitmap_stream**: Performs substring queries with support for GC filtering, expansion, and multithreaded processing

Both programs include an optional shared motif filter module (`motif_filters.hpp`) that provides sequence-level checks.

## Sequence Motif Filtering

barcodesDB includes optional sequence-level filters for short motifs that may be undesirable in some barcode-design contexts. These include:

| Filter | Description | Default |
|--------|------------|---------|
| **Homopolymer filter** | Detects runs of the same base longer than a threshold | max run = 4 |
| **Low-complexity filter** | Computes Shannon entropy over A/C/G/T; flags sequences below threshold | min entropy = 1.5 |
| **Dinucleotide repeat filter** | Detects repeated dinucleotide patterns (≥3 consecutive repeats = ≥6 bp) | Off |
| **Restriction-site filter** | Flags common restriction enzyme recognition sites (EcoRI, BamHI, HindIII, NotI, XhoI, and 10 more) | Off |
| **Functional-motif filter** | Flags short functional-like motifs (TATA-box, polyA signal, splice-donor/acceptor-like, Kozak-like core, CpG) | Off |

**Important caveats:**
- These filters do **not** prove functional neutrality.
- They are **not** organism-specific transcription-factor or splice-site predictions.
- They are simple sequence-level guards intended as practical safeguards.
- All filters are **off by default** — existing website behavior is unchanged.

### Motif modes
- **`off`** (default): No motif checking is performed.
- **`flag`**: Motif hits are reported but sequences are still returned.
- **`exclude`**: Sequences that fail enabled checks are excluded from output.

### For exact k-mer lookup (`query_kmer_bitmap`)
- Motif filters annotate the queried sequence only.
- They do **not** alter the present/absent existence result.
- Filter metadata is appended as extra columns in the output.

### For barcode search (`query_substring_bitmap_stream`)
- In **flag** mode: candidates pass through normally; aggregate metadata is reported in `__META__`.
- In **exclude** mode: candidates that fail filters are skipped; `__META__` reports the filtered count.
- Per-row motif metadata requires `--include-motif-metadata` flag.

## Usage

Access the web interface at [barcodesdb.com](https://barcodesdb.com) to perform k-mer lookups or barcode searches. Input your data and retrieve results instantly.

## Website/query-binary QC

The `qc/` directory contains a comprehensive end-to-end validation suite for the website-facing implementation.

### Run all QC

```bash
bash qc/run_website_qc.sh
```

### Individual QC commands

```bash
# Build binaries and generate synthetic fixtures
bash qc/build_query_binaries.sh
python3 qc/create_synthetic_kbit_fixture.py        # if a Python generator is used
# (or the C++ generator, already compiled by the build script)

# Exact k-mer lookup validation
python3 qc/test_query_kmer_bitmap.py

# Substring/barcode search validation (requires production shards)
python3 qc/test_query_substring_bitmap_stream.py

# GC filtering validation
python3 qc/test_gc_filtering.py

# Pagination/cursor validation (requires production shards)
python3 qc/test_pagination_cursor.py

# Backend/API smoke tests (requires server running on localhost:8090)
python3 qc/test_server_api_smoke.py

# Parameter parity checks
python3 qc/test_server_parameter_parity.py
python3 qc/test_frontend_parameter_parity.py

# Output schema regression tests
python3 qc/test_output_schema.py

# Benchmark smoke tests
python3 qc/benchmark_query_binaries.py

# Random barcode generation QC
python3 qc/test_random_barcode_generation.py
```

### Test coverage

See `qc/QC_TEST_MATRIX.md` for the full list of covered features, including:
- Exact k-mer present/absent lookup across k=16, 17, 18
- Invalid-input and wrong-length handling
- Substring/candidate search with GC filtering
- Pagination and cursor-based navigation
- Optional short-motif filtering (homopolymers, low-complexity, dinucleotide/trinucleotide/tetranucleotide repeats, restriction sites, functional motifs)
- Backend-to-binary parameter parity
- Frontend-to-API parameter parity
- Output schema regression prevention
- Performance benchmark smoke tests

### Notes

- Tests that require production shards are **optional** and are skipped with `[SKIP]` when shard data is not present.
- All generated test data and logs go under `qc/tmp/`, which is gitignored.
- The QC suite uses **Python standard library only** — no external dependencies.

## Building the Programs

To compile the C++ tools, use the following Makefile. Ensure Roaring is installed and adjust paths as necessary.

```makefile
CXX = g++
CXXFLAGS = -O3 -march=native -flto -std=c++17 -pthread -lm

# Adjust these paths to your Roaring installation
ROARING_INCLUDE = /usr/local/include
ROARING_LIB = /usr/local/lib/libroaring.a

all: query_kmer_bitmap query_substring_bitmap_stream

query_kmer_bitmap: query_kmer_bitmap.cpp
	$(CXX) $(CXXFLAGS) -I$(ROARING_INCLUDE) $< $(ROARING_LIB) -o $@

query_substring_bitmap_stream: query_substring_bitmap_stream.cpp
	$(CXX) $(CXXFLAGS) -I$(ROARING_INCLUDE) $< $(ROARING_LIB) -o $@

clean:
	rm -f query_kmer_bitmap query_substring_bitmap_stream
```

Run `make` to build the executables.

## CLI Examples with Motif Filters

### query_kmer_bitmap — Flag mode
```bash
# Annotate queried k-mers with motif warnings (does NOT alter present/absent result)
./query_kmer_bitmap --shards ../shards_18 --k 18 --kmers query_kmers.txt \
  --motif-mode flag \
  --filter-homopolymers --max-homopolymer 4 \
  --filter-low-complexity --min-shannon-entropy 1.5 \
  --filter-dinucleotide-repeats \
  --filter-restriction-sites \
  --filter-functional-motifs
```

### query_kmer_bitmap — Exclude mode (annotate only)
```bash
# Exact lookup: exclude mode sets motif_passes=false but still returns the exact result
./query_kmer_bitmap --shards ../shards_18 --k 18 --kmers query_kmers.txt \
  --motif-mode exclude \
  --filter-homopolymers --filter-restriction-sites
```

### query_substring_bitmap_stream — Exclude mode
```bash
# Exclude candidates that fail motif filters; __META__ shows filtered count
./query_substring_bitmap_stream --shards ../shards_18 --gc-hist ../gc_hist_shards_18.json \
  --gc-min 40 --gc-max 60 --limit 200 --threads 16 \
  --motif-mode exclude \
  --filter-homopolymers \
  --filter-low-complexity \
  --filter-restriction-sites
```

### query_substring_bitmap_stream — Flag mode with metadata
```bash
# Flag motif hits and include per-row metadata in output
./query_substring_bitmap_stream --shards ../shards_18 --gc-hist ../gc_hist_shards_18.json \
  --gc-min 30 --gc-max 70 --limit 100 \
  --motif-mode flag \
  --filter-homopolymers \
  --filter-functional-motifs \
  --include-motif-metadata
```

## API Examples with Motif Filters

### K-mer Lookup (POST /api/query-kmer)
```json
{
  "kmers": ["AAAAAAAAAAAAAAAA", "ACGTACGTACGTACGT"],
  "motif_mode": "exclude",
  "filter_homopolymers": true,
  "max_homopolymer": 4,
  "filter_low_complexity": true,
  "min_shannon_entropy": 1.5,
  "filter_dinucleotide_repeats": false,
  "filter_restriction_sites": true,
  "filter_functional_motifs": false
}
```

### Barcode Search (POST /api/query-substring)
```json
{
  "gcMin": 40,
  "gcMax": 60,
  "limit": 200,
  "motif_mode": "exclude",
  "filter_homopolymers": true,
  "filter_restriction_sites": true
}
```</content>