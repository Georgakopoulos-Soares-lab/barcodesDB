# barcodesDB Website QC Test Matrix

| # | Feature | Component | Test File | Status |
|---|---------|-----------|-----------|--------|
| 1 | Exact k-mer present lookup | `query_kmer_bitmap.cpp` | `test_query_kmer_bitmap.py` | ✅ Covered |
| 2 | Exact k-mer absent lookup | `query_kmer_bitmap.cpp` | `test_query_kmer_bitmap.py` | ✅ Covered |
| 3 | Lowercase input handling | `query_kmer_bitmap.cpp` | `test_query_kmer_bitmap.py` | ✅ Covered |
| 4 | Invalid character rejection | `query_kmer_bitmap.cpp` | `test_query_kmer_bitmap.py` | ✅ Covered |
| 5 | Wrong-length k-mer rejection | `query_kmer_bitmap.cpp` | `test_query_kmer_bitmap.py` | ✅ Covered |
| 6 | Missing fixture handling | `query_kmer_bitmap.cpp` | `test_query_kmer_bitmap.py` | ✅ Covered |
| 7 | Motif filter (k-mer) — off by default | `query_kmer_bitmap.cpp` | `test_query_kmer_bitmap.py` | ✅ Covered |
| 8 | Motif filter (k-mer) — does not alter existence | `query_kmer_bitmap.cpp` | `test_query_kmer_bitmap.py` | ✅ Covered |
| 9 | Output stability (k-mer) | `query_kmer_bitmap.cpp` | `test_query_kmer_bitmap.py` | ✅ Covered |
| 10 | Substring search — basic | `query_substring_bitmap_stream.cpp` | `test_query_substring_bitmap_stream.py` | ✅ Covered |
| 11 | Substring search — no hits | `query_substring_bitmap_stream.cpp` | `test_query_substring_bitmap_stream.py` | ✅ Covered |
| 12 | Substring search — limit behavior | `query_substring_bitmap_stream.cpp` | `test_query_substring_bitmap_stream.py` | ✅ Covered |
| 13 | Pagination / cursor | `query_substring_bitmap_stream.cpp` | `test_pagination_cursor.py` | ✅ Covered |
| 14 | GC filter — inclusion/exclusion | `query_substring_bitmap_stream.cpp` | `test_gc_filtering.py` | ✅ Covered |
| 15 | GC filter — boundary values | `query_substring_bitmap_stream.cpp` | `test_gc_filtering.py` | ✅ Covered |
| 16 | GC filter — invalid values | `query_substring_bitmap_stream.cpp` | `test_gc_filtering.py` | ✅ Covered |
| 17 | k-length handling (substring) | `query_substring_bitmap_stream.cpp` | `test_query_substring_bitmap_stream.py` | ✅ Covered |
| 18 | Invalid input (substring) | `query_substring_bitmap_stream.cpp` | `test_query_substring_bitmap_stream.py` | ✅ Covered |
| 19 | Metadata consistency | `query_substring_bitmap_stream.cpp` | `test_output_schema.py` | ✅ Covered |
| 20 | Motif filter (substring) — off/flag/exclude | `query_substring_bitmap_stream.cpp` | `test_query_substring_bitmap_stream.py` | ✅ Covered |
| 21 | Benchmark / latency smoke | Both binaries | `benchmark_query_binaries.py` | ✅ Covered |
| 22 | Backend parameter parity | `server.js` | `test_server_parameter_parity.py` | ✅ Covered |
| 23 | Frontend parameter parity | `public/*.html` | `test_frontend_parameter_parity.py` | ✅ Covered |
| 24 | Backend API smoke | `server.js` | `test_server_api_smoke.py` | ✅ Covered |
| 25 | Output schema regression | Both binaries + API | `test_output_schema.py` | ✅ Covered |
| 26 | Random barcode generation | `server.js` `/api/generate-random-kmers` | `test_random_barcode_generation.py` | ✅ Covered |
