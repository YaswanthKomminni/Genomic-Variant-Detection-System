"""
test_cases.py

Automated test cases for the DSA case study. This file runs several scenarios
and prints outputs to demonstrate how each data structure works in the pipeline.
"""

from genomic_system import GenomicSystem

def run_tests():
    system = GenomicSystem()

    print("=== Initial known variants count:", len(system.list_known()))

    # Test 1: Known variant (should be probable 'known' via Bloom + found in HashTable)
    system.add_case(1, "Patient_A", "VAR_1_mut", 1500)

    # Test 2: Known variant, different position
    system.add_case(2, "Patient_B", "VAR_10_mut", 3200)

    # Test 3: Unknown variant (not in initial JSON) -> Bloom says not known
    system.add_case(1, "Patient_C", "NEW_VARIANT_X", 2100)

    # Test 4: Insert many cases to demonstrate priority ordering and BST growth
    for i in range(5):
        variant = f"VAR_{20 + i}_mut"
        system.add_case(3 + i, f"Patient_bulk_{i}", variant, 2500 + i*50)

    # Test 5: Range query demonstration (manual)
    print("\n-- Range Query 2000-2600 --")
    hits = system.bst.range_query(2000, 2600)
    for pos, var in hits:
        print("  ", pos, var)

    # Save data at the end of tests
    system.save_all()
    print("\nSaved data; tests complete.")

if __name__ == '__main__':
    run_tests()
