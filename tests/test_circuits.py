import sys
sys.path.append("./")
from quaos.paulis import random_pauli_string
from quaos.circuits.known_circuits import to_x, to_ix
# TODO: Make these proper tests


def test_to_x():
    target_x = 0
    list_of_failures = []
    for _ in range(20000):
        ps = random_pauli_string([3, 3, 3, 3])
        if ps.n_identities() == 4:
            continue
        c = to_x(ps, 0)
        if c.act(ps).x_exp[target_x] == 0 or c.act(ps).z_exp[target_x] != 0:
            print(f"Failed: {ps} -> {c.act(ps)}")
            list_of_failures.append(ps)

    return list_of_failures


def test_to_ix():
    target_x = 0
    list_of_failures = []
    for _ in range(2000):
        ps = random_pauli_string([3, 3, 3, 3])
        if ps.n_identities() == 4:
            continue
        c = to_ix(ps, 0)
        if c is None:
            print(f"Failed: {ps} -> {c}")
            list_of_failures.append(ps)
            continue
        failed = False
        for i in range(ps.n_qudits()):
            if i == target_x and failed is False:
                if c.act(ps).x_exp[target_x] == 0 or c.act(ps).z_exp[target_x] != 0:
                    print(f"Failed target x: {ps} -> {c.act(ps)}")
                    list_of_failures.append(ps)
                    failed = True
            elif failed is False:
                if c.act(ps).x_exp[i] != 0 or c.act(ps).z_exp[i] != 0:
                    print(f"Failed identity: {ps} -> {c.act(ps)}")
                    list_of_failures.append(ps)
                    failed = True

    return list_of_failures


if __name__ == "__main__":
    print(test_to_x())
    print(test_to_ix())
