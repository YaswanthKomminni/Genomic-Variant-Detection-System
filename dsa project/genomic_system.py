"""
genomic_system.py

Integrates BloomFilter, PriorityQueue, VariantMetadata, and BST.
Designed as a DSA case-study: emphasis on how each data structure stores data and is used.
"""

import json
import os
from data_structures import BloomFilter, PriorityQueue, VariantMetadata, BST

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
KNOWN_FILE = os.path.join(DATA_DIR, "known_variants.json")
META_FILE = os.path.join(DATA_DIR, "metadata.json")
CASE_FILE = os.path.join(DATA_DIR, "processed_cases.json")

class GenomicSystem:
    def __init__(self, bloom_size=10000):
        # Data structure instances
        self.bloom = BloomFilter(size=bloom_size, hash_count=3)
        self.queue = PriorityQueue()
        self.metadata = VariantMetadata()
        self.bst = BST()

        # persistent storage holders
        self.known_variants = self._load_json(KNOWN_FILE, [])
        self.metadata.data = self._load_json(META_FILE, {})
        self.cases_log = self._load_json(CASE_FILE, [])

        # initialize bloom from known variants list
        for v in self.known_variants:
            self.bloom.add(v)

    # ---------- File helpers ----------
    def _load_json(self, path, default):
        if os.path.exists(path):
            try:
                with open(path, "r") as f:
                    content = f.read().strip()
                    if not content:
                        return default
                    return json.loads(content)
            except json.JSONDecodeError:
                print(f"Warning: corrupted JSON file {path} — using default")
                return default
        return default

    def _save_json(self, path, data):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            json.dump(data, f, indent=2)

    def save_all(self):
        self._save_json(KNOWN_FILE, self.known_variants)
        self._save_json(META_FILE, self.metadata.data)
        self._save_json(CASE_FILE, self.cases_log)

    # ---------- Core operations ----------
    def add_known_variant(self, variant, info):
        # store in hash table and bloom and list for persistence
        if variant not in self.known_variants:
            self.known_variants.append(variant)
            self.metadata.add_metadata(variant, info)
            self.bloom.add(variant)
            print(f"Added known variant: {variant}")
        else:
            # update metadata if provided
            if info:
                self.metadata.add_metadata(variant, info)
                print(f"Updated metadata for {variant}")

    def add_case(self, priority, patient_id, variant, position):
        case = {"patient": patient_id, "variant": variant, "position": position, "priority": priority}
        self.queue.push(priority, case)
        # immediately process queue for this simple CLI-based case study
        self.process_cases()

    def process_cases(self):
        while not self.queue.is_empty():
            priority, case = self.queue.pop()
            variant = case["variant"]
            position = case["position"]
            print(f"Processing case: {case['patient']} (priority={priority})")

            # Bloom check (fast membership test)
            possibly_known = self.bloom.contains(variant)
            if possibly_known:
                print(f" - BloomFilter: POSSIBLY known: {variant}")
                meta = self.metadata.get_metadata(variant)
                if meta:
                    print(f"   -> Metadata (HashTable): {meta}")
                else:
                    print("   -> No metadata found in HashTable")
            else:
                print(f" - BloomFilter: NOT known: {variant}")

            # Insert into BST for positional queries
            self.bst.insert(position, variant)
            # log the case into persistent log
            self.cases_log.append(case)
            print(f" - Inserted into BST at position {position} and logged case.")

            # show a small range query example after each insertion
            low, high = position - 500, position + 500
            hits = self.bst.range_query(low, high)
            print(f" - Range query around [{low},{high}] -> {len(hits)} hits (showing up to 5): {hits[-5:]}")

    # helpers for testing & demonstration
    def list_known(self):
        return list(self.known_variants)

    def get_metadata(self, variant):
        return self.metadata.get_metadata(variant)

    def list_cases(self):
        return list(self.cases_log)
