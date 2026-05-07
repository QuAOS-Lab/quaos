import random
from sympleq.graph_utils import find_one_permutation, preserves_per_label_multisets, mapping_key


class TestGraphUtils:
    def test_find_one_permutation(self, num_cases: int = 2000, seed: int | None = None):
        """
        Randomized tests that enforce:
        - No duplicates within a label
        - Lists have unique elements
        Asserts that any permutation found preserves per-label multisets.
        """
        random.seed(seed)
        n_permutations_found = 0
        for _ in range(num_cases):
            # Build a small random instance
            # Choose 2–4 labels from a small pool of complex numbers
            label_pool = [complex(a, b) for a in range(2) for b in range(2)]  # 0+0j, 0+1j, 1+0j, 1+1j
            num_labels = random.randint(1, 3)
            labels_chosen = random.sample(label_pool, num_labels)

            lists = []
            labels = []
            universe = list(range(1, 9))  # element universe

            for lab in labels_chosen:
                # 1–3 lists per label
                m = random.randint(1, 3)
                seen = set()  # avoid duplicates within this label
                for _ in range(m):
                    k = random.randint(2, 4)  # list size
                    lst = tuple(sorted(random.sample(universe, k)))
                    # ensure no duplicate lists within this label
                    while lst in seen:
                        lst = tuple(sorted(random.sample(universe, k)))
                    seen.add(lst)
                    lists.append(list(lst))
                    labels.append(lab)

            # Run search
            found = find_one_permutation(lists, labels, previously_found=set(), max_cycle_size=4)
            if found is not None:
                # Verify correctness strictly per label
                assert preserves_per_label_multisets(lists, labels, found), \
                    f"Invalid permutation returned: {found}\nlists={lists}\nlabels={labels}"

                # Verify exclude works
                exclude = {mapping_key(found, domain=sorted({x for lst in lists for x in lst}))}
                found2 = find_one_permutation(lists, labels, previously_found=exclude, max_cycle_size=3)
                if found2 is not None:
                    # If another is found, it must be different under the canonical key
                    assert mapping_key(found2) not in exclude, "Exclude list not respected."
                n_permutations_found += 1
            # else: no permutation found is perfectly fine
