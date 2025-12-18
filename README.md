
# DSA_Case_Study: Genomic Variant Processing System (Creative Conceptual Theme)

## Project Overview
This project demonstrates how multiple data structures work together to manage genomic variant data.
Primary focus: **Data Structures** (Bloom Filter, Priority Queue, Hash Table, BST) and their role in the pipeline.
UI is intentionally omitted to emphasize DSA concepts.

## Repository Structure
- data_structures.py  -> implementations of BloomFilter, PriorityQueue, VariantMetadata, BST
- genomic_system.py   -> integration and persistence logic (uses JSON files in /data)
- test_cases.py       -> automated test cases that demonstrate behavior
- data/               -> initial JSON files (known_variants.json, metadata.json, processed_cases.json)
- presentation.pptx   -> Creative conceptual PPT (not auto-generated)
- DSA_Case_Study.pdf  -> Detailed report (not auto-generated)

## How to run
1. Ensure you are in the repository root (where these files are located)
2. Run tests:
```
python3 test_cases.py
```
3. Or run interactive (import GenomicSystem in a REPL):
```
python3
>>> from genomic_system import GenomicSystem
>>> s = GenomicSystem()
>>> s.add_case(1, 'P1', 'VAR_1_mut', 1500)
>>> s.save_all()
```

- Emphasis is on correct usage of data structures, their reasoning and complexity
- BloomFilter: O(k) insertion/check (k = number of hashes), memory efficient
- PriorityQueue: O(log n) insert/pop
- Hash Table: O(1) average lookup/insert
- BST (unbalanced): average O(log n) insert/search for random data; worst-case O(n)

## Team contributions suggestion
- Karthikeya: Bloom Filter + explanation slides
- Sairam: Priority Queue + test cases
- Mouliiswar: Hash Table + JSON persistence
- Yaswanth: BST + range query demonstration + README

