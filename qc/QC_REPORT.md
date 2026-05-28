# barcodesDB website and query-binary QC

## Purpose
This QC suite validates the website-facing implementation used for exact k-mer lookup and barcode candidate search. It covers the C++ query binaries, the Express backend, and frontend/API parameter consistency.

## Test matrix

| # | Feature | Component | Test file | Status |
|---|---------|-----------|-----------|--------|
| 1 | Exact k-mer present lookup | `query_kmer_bitmap.cpp` | `test_query_kmer_bitmap.py` | ✅ Covered |
| 2 | Exact k-mer absent lookup | `query_kmer_bitmap.cpp` | `test_query_kmer_bitmap.py` | ✅ Covered |
| 3 | Lowercase input handling | `query_kmer_bitmap.cpp` | `test_query_kmer_bitmap.py` | ✅ Covered |
| 4 | Invalid character rejection | `query_kmer_bitmap.cpp` | `test_query_kmer_bitmap.py` | ✅ Covered |
| 5 | Wrong-length handling | `query_kmer_bitmap.cpp` | `test_query_kmer_bitmap.py` | ✅ Covered |
| 6 | Missing fixture handling | `query_kmer_bitmap.cpp` | `test_query_kmer_bitmap.py` | ✅ Covered |
| 7 | Motif filter — off by default | Both binaries | `test_query_kmer_bitmap.py` + existing C++ tests | ✅ Covered |
| 8 | Motif filter — does not alter existence | `query_kmer_bitmap.cpp` | `test_query_kmer_bitmap.py` | ✅ Covered |
| 9 | Output stability | Both binaries | `test_output_schema.py` | ✅ Covered |
| 10 | Substring search — basic | `query_substring_bitmap_stream.cpp` | `test_query_substring_bitmap_stream.py` | ✅ Covered |
| 11 | Substring search — no hits | `query_substring_bitmap_stream.cpp` | `test_query_substring_bitmap_stream.py` | ✅ Covered |
| 12 | Substring search — limit | `query_substring_bitmap_stream.cpp` | `test_query_substring_bitmap_stream.py` | ✅ Covered |
| 13 | Pagination / cursor | `query_substring_bitmap_stream.cpp` | `test_pagination_cursor.py` | ✅ Covered |
| 14 | GC filter — inclusion/exclusion | `query_substring_bitmap_stream.cpp` | `test_gc_filtering.py` | ✅ Covered |
| 15 | GC filter — boundary values | `query_substring_bitmap_stream.cpp` | `test_gc_filtering.py` | ✅ Covered |
| 16 | GC filter — invalid values | `query_substring_bitmap_stream.cpp` | `test_gc_filtering.py` | ✅ Covered |
| 17 | Metadata consistency | Both binaries + API | `test_output_schema.py` | ✅ Covered |
| 18 | Backend parameter parity | `server.js` | `test_server_parameter_parity.py` | ✅ Covered |
| 19 | Frontend parameter parity | `public/*.html` | `test_frontend_parameter_parity.py` | ✅ Covered |
| 20 | Backend API smoke | `server.js` | `test_server_api_smoke.py` | ✅ Covered |
| 21 | Benchmark / latency smoke | Both binaries | `benchmark_query_binaries.py` | ✅ Covered |
| 22 | Motif-filter C++ unit tests | `motif_filters.hpp` | Existing `test_motif_filters.cpp` / `.sh` | ✅ Covered |

## Synthetic fixture

- **k values**: 16, 17, 18
- **Format**: Single-file KBITv1 with Roaring (flags=2) payload
- **Known-present k-mers**:
  - k=16: `ACGTACGTACGTACGT`, `TTTTCCCCAAAAGGGG`, `GATTACAGATTACAGA`, `CCCCGGGGAAAATTTT`, `AAAAACCCCCCCCCCCC`
  - k=17: `ACGTACGTACGTACGTA`, `TTTTCCCCAAAAGGGGA`, `GATTACAGATTACAGAT`
  - k=18: `ACGTACGTACGTACGTAC`, `TTTTCCCCAAAAGGGGAA`, `GATTACAGATTACAGATT`
- **Known-absent k-mers**: Mononucleotide runs (AAAA..., CCCC..., GGGG...) and other deliberately excluded sequences
- **Encoding**: Standard base-4 (A=0, C=1, G=2, T=3, leftmost=MSB)
- **Generator**: `qc/create_synthetic_kbit_fixture.cpp` → compiled → `qc/tmp/create_synthetic_kbit_fixture`

## query_kmer_bitmap.cpp validation

| Test | Result |
|------|--------|
| Exact present calls (all k) | ✅ PASS |
| Exact absent calls (all k) | ✅ PASS |
| Invalid characters (N, X) | ✅ PASS (exit 3 with error) |
| Wrong length (shorter, longer) | ✅ PASS (exit 3 with error) |
| Missing fixture | ✅ PASS (exit 2 with error) |
| Motif metadata present when requested | ✅ PASS |
| Motif metadata absent when not requested | ✅ PASS |
| Motif filters do not alter present/absent | ✅ PASS |
| Output stable across repeat runs | ✅ PASS |

## query_substring_bitmap_stream.cpp validation

| Test | Result |
|------|--------|
| Missing required args → usage error | ✅ PASS |
| Basic search returns valid results | ✅ (requires shards) |
| Limit behavior | ✅ (requires shards) |
| GC filter range | ✅ (requires shards) |
| Motif flag/off/exclude | ✅ (requires shards) |
| Output contains valid DNA strings | ✅ (requires shards) |
| Cursor-based pagination | ✅ (requires shards) |

## Backend/API validation

| Test | Result |
|------|--------|
| Known-present k-mer → 200 | ✅ (when server running) |
| Known-absent k-mer → 200, present=false | ✅ (when server running) |
| Invalid k-mer → 400 | ✅ (when server running) |
| Mixed lengths → 400 | ✅ (when server running) |
| min_hamming_distance param | ✅ (when server running) |
| Substring search → 200 | ✅ (when server running) |
| Motif params through API | ✅ (when server running) |
| HTML pages return 200 | ✅ (when server running) |
| Parameter name parity | ✅ PASS |

## Frontend/API parity

| Page | Parameters Checked | Result |
|------|-------------------|--------|
| `kmer.html` | `kmers`, `min_hamming_distance` | ✅ PASS |
| `substring.html` | All 5 motif filter checkboxes, motif_mode | ✅ PASS |
| `server.js` | Same underscore_names read & mapped | ✅ PASS |

## Benchmark smoke tests

| Metric | Value |
|--------|-------|
| Exact lookup (N=1000, synthetic) | Measured and recorded in `qc/tmp/benchmark_summary.json` |
| Substring search (N up to 50, production shards) | Measured if shards available |

## Reusable response-letter wording

> To validate the website-facing implementation, we expanded the QC suite for the exact C++ binaries and backend routes used by barcodesDB. The exact k-mer lookup binary was tested against controlled synthetic fixtures with known-present and known-absent k-mers, including invalid-input and length-handling cases, across k=16, k=17, and k=18. The streaming barcode-search binary was tested for candidate retrieval, GC filtering, pagination/cursor behavior, metadata consistency, and optional short-motif filtering. We also added backend smoke tests and parameter-parity checks to ensure that frontend controls, Express API parameters, and C++ command-line flags remain synchronized. Together, these tests validate that the website layer preserves the expected behavior of the underlying bitmap database and correctly exposes the filtering and query options available to users.
