from gates import Circuit, SUM as CX, PHASE as S, Hadamard as H
from symplectic import Pauli, PauliString, PauliSum, symplectic_product
import numpy as np
from sympy import Matrix, mod_inverse
import random

# def check_mappable_via_clifford(pauli_sum: PauliSum, target_pauli_sum: PauliSum):


def check_mappable_via_clifford(pauli_sum: PauliSum, target_pauli_sum: PauliSum) -> bool:
    if pauli_sum.n_qudits() != target_pauli_sum.n_qudits():
        raise ValueError("Pauli sums must have the same number of qudits")
    if np.any(pauli_sum.dimensions != target_pauli_sum.dimensions):
        raise ValueError("Pauli sums must have the same dimensions")
    if pauli_sum.n_paulis() != target_pauli_sum.n_paulis():
        raise ValueError("Pauli sums must have the same number of paulis")
    return np.all(pauli_sum.symplectic_product_matrix() == target_pauli_sum.symplectic_product_matrix())


def find_circuit(pauli_sum: PauliSum, target_pauli_sum: PauliSum, iterations: int,
                 compare_phases: bool = True, stop_if_found: bool = True) -> list[Circuit]:
    
    mappable = check_mappable_via_clifford(pauli_sum, target_pauli_sum)
    if not mappable:
        return []

    n_qudits = pauli_sum.n_qudits()
    SUMs = [CX(i, j, 2) for i in range(n_qudits) for j in range(n_qudits) if i != j]
    Ss = [S(i, 2) for i in range(n_qudits)]
    Hs = [H(i, 2) for i in range(n_qudits)]
    all_gates = SUMs + Ss + Hs

    goal_circuits = []
    circuits = [Circuit(dimensions=pauli_sum.dimensions)]
    intermediate_paulis = [pauli_sum.copy()]
    for i in range(iterations):
        print(i, len(intermediate_paulis), len(goal_circuits))
        intermediate_paulis_old = intermediate_paulis.copy()
        for i, p in enumerate(intermediate_paulis_old):
            for g in all_gates:
                P_temp = g.act(p)
                if not compare_phases:
                    P_temp.phases = np.zeros(P_temp.n_paulis(), dtype=int)
                C_temp = Circuit(dimensions=[2 for i in range(n_qudits)])
                for g2 in circuits[i].gates:
                    C_temp.add_gate(g2)
                C_temp.add_gate(g)
                if P_temp not in intermediate_paulis:
                    intermediate_paulis.append(P_temp)
                    circuits.append(C_temp)
                if P_temp == target_pauli_sum:
                    goal_circuits.append(C_temp)
                    if stop_if_found:
                        return [C_temp]
    return goal_circuits


def find_agnostic_circuit(pauli_sum: PauliSum, target_paulis: list[Pauli], target_indexes: list[tuple[int, int]],
                          iterations: int, compare_phases: bool = True) -> Circuit:
    """
    Find a circuit that maps pauli_sum to target_paulis with target_indexes. Gives the first circuit found that does
    this, agnostic to any Pauli not in the target list.
    """
    n_trials = 100

    current_target = pauli_sum.copy()  # the target starts as a copy of the initial pauli_sum with target paulis changed only
    for i, indexes in enumerate(target_indexes):
        current_target[indexes[0], indexes[1]] = target_paulis[i]  
    
    for max_iter in range(iterations):

        mappable = check_mappable_via_clifford(pauli_sum, current_target)
        if not mappable:
            # change a pauli in the target - loop through the target pauli_strings and target qudits for change locations
            print('Changing target state')
            pass
        else:
            print('finding circuit with naive state' )
            C = find_circuit(pauli_sum, current_target, max_iter, compare_phases, True)
            if C != []:
                return C[0]
        
        paulis = [Pauli(0, 0, 2), Pauli(0, 1, 2), Pauli(1, 0, 2), Pauli(1, 1, 2)]  # only works for qubits for now
        agnostic_pauli_locations = [(i, j) for i in range(pauli_sum.n_paulis()) for j in range(pauli_sum.n_qudits()) if (i, j) not in target_indexes]
        # we make random changes to the target until we find one that works
        for trial_number in range(n_trials):
            ps, q = random.choice(agnostic_pauli_locations)  # pick a random pauli string and qudit to change
            p = random.choice(paulis)  # pick a random pauli
            current_target[ps, q] = p
            mappable = check_mappable_via_clifford(pauli_sum, current_target)
            if mappable:
                C = find_circuit(pauli_sum, current_target, max_iter, compare_phases, True)
                if C != []:
                    print('Found on trial number ', trial_number + 1)
                    return C[0]
                
    raise Exception('No circuit found')
        
    
if __name__ == "__main__":
    initial_pauli = ['x0z0 x0z1', 'x1z1 x1z0']
    goal_pauli = ['x1z1 x1z1', 'x1z1 x0z1']
    initial_pauli = PauliSum(initial_pauli, dimensions=[2, 2])
    goal_pauli = PauliSum(goal_pauli, dimensions=[2, 2])
    print(initial_pauli)

    # C = find_circuit(initial_pauli, goal_pauli, 8)
    # print(C)

    goal_single_pauli = Pauli(1, 0, 2)
    goal_location = (0, 1)

    C = find_agnostic_circuit(initial_pauli, [goal_single_pauli], [goal_location], 8)
    print(C)
    print(initial_pauli)
    print(C.act(initial_pauli))
