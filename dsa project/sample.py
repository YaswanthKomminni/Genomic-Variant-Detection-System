import heapq
import hashlib
import json
import os

# ---------------------------
# 1. Bloom Filter
# ---------------------------
class BloomFilter:
    def _init_(self, size=10000, hash_count=3):
        self.size = size
        self.hash_count = hash_count
        self.bit_array = [0] * size

    def _hashes(self, item):
        results = []
        for i in range(self.hash_count):
            h = int(hashlib.md5((item + str(i)).encode()).hexdigest(), 16)
            results.append(h % self.size)
        return results

    def add(self, item):
        for h in self._hashes(item):
            self.bit_array[h] = 1

    def contains(self, item):
        return all(self.bit_array[h] == 1 for h in self._hashes(item))

# ---------------------------
# 2. Priority Queue
# ---------------------------
class PriorityQueue:
    def _init_(self):
        self.heap = []

    def push(self, priority, case):
        heapq.heappush(self.heap, (priority, case))

    def pop(self):
        return heapq.heappop(self.heap) if self.heap else None

    def is_empty(self):
        return len(self.heap) == 0

# ---------------------------
# 3. Variant Metadata (Hash Table)
# ---------------------------
class VariantMetadata:
    def _init_(self):
        self.data = {}

    def add_metadata(self, variant, info):
        self.data[variant] = info

    def get_metadata(self, variant):
        return self.data.get(variant, None)

    def show_all(self):
        print("\n--- Known Variants & Metadata ---")
        if not self.data:
            print("No known variants.")
            return
        for variant, info in self.data.items():
            print(f"{variant}: {info}")

# ---------------------------
# 4. Binary Search Tree
# ---------------------------
class BSTNode:
    def _init_(self, key, variant):
        self.key = key
        self.variants = [variant]  # list of variants at same position
        self.left = None
        self.right = None

class BST:
    def _init_(self):
        self.root = None

    def insert(self, key, variant):
        self.root = self._insert(self.root, key, variant)

    def _insert(self, node, key, variant):
        if not node:
            return BSTNode(key, variant)
        if key < node.key:
            node.left = self._insert(node.left, key, variant)
        elif key > node.key:
            node.right = self._insert(node.right, key, variant)
        else:
            node.variants.append(variant)
        return node

    def inorder(self):
        result = []
        self._inorder(self.root, result)
        return result

    def _inorder(self, node, result):
        if node:
            self._inorder(node.left, result)
            for var in node.variants:
                result.append((node.key, var))
            self._inorder(node.right, result)

    def range_query(self, low, high):
        result = []
        self._range_query(self.root, low, high, result)
        return result

    def _range_query(self, node, low, high, result):
        if not node:
            return
        if low <= node.key <= high:
            for var in node.variants:
                result.append((node.key, var))
        if low < node.key:
            self._range_query(node.left, low, high, result)
        if node.key < high:
            self._range_query(node.right, low, high, result)

# ---------------------------
# 5. Genomic System
# ---------------------------
class GenomicSystem:
    VARIANT_FILE = "known_variants.json"
    META_FILE = "metadata.json"
    CASE_LOG_FILE = "cases_log.json"

    def _init_(self):
        self.bloom = BloomFilter(size=10000, hash_count=3)
        self.queue = PriorityQueue()
        self.metadata = VariantMetadata()
        self.bst = BST()

        self.known_variants = self.load_json(self.VARIANT_FILE, [])
        self.metadata.data = self.load_json(self.META_FILE, {})
        self.cases_log = self.load_json(self.CASE_LOG_FILE, [])

        print(f"\n✅ Loaded {len(self.known_variants)} known variants.")
        print(f"✅ Loaded {len(self.cases_log)} previous cases.\n")

        for var in self.known_variants:
            self.bloom.add(var)

    # -----------------------
    # File Operations
    # -----------------------
    def load_json(self, filename, default):
        if os.path.exists(filename):
            with open(filename, "r") as f:
                return json.load(f)
        return default

    def save_json(self, filename, data):
        with open(filename, "w") as f:
            json.dump(data, f, indent=2)

    def save_all(self):
        self.save_json(self.VARIANT_FILE, self.known_variants)
        self.save_json(self.META_FILE, self.metadata.data)
        self.save_json(self.CASE_LOG_FILE, self.cases_log)
        print("\n💾 All data saved successfully!")

    # -----------------------
    # Core Logic
    # -----------------------
    def add_known_variant(self, variant, info):
        if variant not in self.known_variants:
            self.known_variants.append(variant)
            self.metadata.add_metadata(variant, info)
            self.bloom.add(variant)
            print(f"✅ Added new known variant: {variant}")
        else:
            self.metadata.add_metadata(variant, info)
            print(f"⚠ Variant {variant} updated.")

    def add_case(self, patient_id, variant, position):
        # Check if variant is known
        if self.bloom.contains(variant):
            info = self.metadata.get_metadata(variant)
            if info and "risk_score" in info:
                risk = info["risk_score"]
                print(f"\n⚡ Known variant found: {variant}")
                print(f"   → Disease: {info['disease']}")
                print(f"   → Drug: {info['drug']}")
                print(f"   → Risk Score: {risk}")
            else:
                print(f"\n⚠ Known variant {variant} has no risk info.")
                risk = float(input("Enter risk score (0.0–1.0): "))
                info["risk_score"] = risk
                self.metadata.add_metadata(variant, info)
        else:
            print(f"\n🧬 Unknown variant: {variant}")
            disease = input("Enter disease name: ")
            drug = input("Enter drug name: ")
            risk = float(input("Enter risk score (0.0–1.0): "))
            self.add_known_variant(variant, {"disease": disease, "drug": drug, "risk_score": risk})

        # Compute priority
        priority = max(1, 5 - int(risk * 4))
        print(f"🩺 Auto-priority assigned based on risk: {priority}")

        # Add to priority queue
        case = {"patient": patient_id, "variant": variant, "position": position, "priority": priority}
        self.queue.push(priority, case)
        self.process_cases()

    def process_cases(self):
        while not self.queue.is_empty():
            priority, case = self.queue.pop()
            variant = case["variant"]
            position = case["position"]

            print(f"\n🔬 Processing {case['patient']} (Priority {priority})")

            info = self.metadata.get_metadata(variant)
            if info:
                print(f"✅ Variant {variant} metadata: {info}")

            # Insert into BST
            self.bst.insert(position, variant)

            # Log case
            self.cases_log.append(case)
            print(f"   → Stored in BST at position {position}")
            print(f"   → Case logged successfully.")

    # -----------------------
    # Display Functions
    # -----------------------
    def show_all_cases(self):
        print("\n📜 --- Processed Case Log ---")
        if not self.cases_log:
            print("No cases logged yet.")
            return
        sorted_cases = sorted(self.cases_log, key=lambda c: c['priority'])
        for case in sorted_cases[-10:]:
            print(f"Priority: {case['priority']} | Patient: {case['patient']} | Variant: {case['variant']} | Pos: {case['position']}")

    def range_query_variants(self):
        try:
            low = int(input("Enter start position: "))
            high = int(input("Enter end position: "))
        except ValueError:
            print("❌ Invalid input. Enter integers only.")
            return

        results = self.bst.range_query(low, high)
        print(f"\n📊 Variants in range {low}-{high}:")
        if not results:
            print("No variants found in this range.")
        for pos, var in results:
            info = self.metadata.get_metadata(var)
            print(f"{var} at {pos} | Metadata: {info}")

# ---------------------------
# 6. Interactive Menu
# ---------------------------
if __name__ == "_main_":
    system = GenomicSystem()

    while True:
        print("\n===============================")
        print("1. Add Patient Case")
        print("2. Show Case Log")
        print("3. Show Known Variants & Metadata")
        print("4. Query Variants in a Range")
        print("5. Save & Exit")
        print("===============================")
        choice = input("Choose an option: ").strip()

        if choice == "1":
            patient = input("Enter patient name: ")
            variant = input("Enter variant name: ")
            position = int(input("Enter position (e.g., 1500): "))
            system.add_case(patient, variant, position)

        elif choice == "2":
            system.show_all_cases()

        elif choice == "3":
            system.metadata.show_all()

        elif choice == "4":
            system.range_query_variants()

        elif choice == "5":
            system.save_all()
            print("👋 Exiting program.")
            break

        else:
            print("Invalid choice. Try again.")