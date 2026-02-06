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

## Technical Overview

- **Backend**: C++17 applications utilizing Roaring bitmap libraries
- **Frontend**: Responsive web interface using vanilla JavaScript
- **APIs**: RESTful endpoints for programmatic integration
- **Data Format**: KBITv1 format supporting dense and compressed bitmaps
- **Scalability**: Supports k-mer lengths up to 31

## Core Programs

barcodesDB relies on two primary C++ programs:

- **gen_kmer_bitmap**: Generates k-mer existence bitmaps with configurable density and reproducible random sampling
- **query_substring_bitmap_stream**: Performs substring queries with support for GC filtering, expansion, and multithreaded processing

## Usage

Access the web interface at [barcodesdb.com](https://barcodesdb.com) to perform k-mer lookups or barcode searches. Input your data and retrieve results instantly.

## Building the Programs

To compile the C++ tools, use the following Makefile. Ensure Roaring is installed and adjust paths as necessary.

```makefile
CXX = g++
CXXFLAGS = -O3 -march=native -flto -std=c++17 -pthread -lm

# Adjust these paths to your Roaring installation
ROARING_INCLUDE = /usr/local/include
ROARING_LIB = /usr/local/lib/libroaring.a

all: gen_kmer_bitmap query_substring_bitmap_stream

gen_kmer_bitmap: gen_kmer_bitmap.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

query_substring_bitmap_stream: query_substring_bitmap_stream.cpp
	$(CXX) $(CXXFLAGS) -I$(ROARING_INCLUDE) $< $(ROARING_LIB) -o $@

clean:
	rm -f gen_kmer_bitmap query_substring_bitmap_stream
```

Run `make` to build the executables.</content>