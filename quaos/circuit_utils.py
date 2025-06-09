from gates import Circuit, SUM as CX, PHASE as S, Hadamard as H, GateOperation
from paulis import Pauli, PauliString, PauliSum, symplectic_product
import numpy as np
import sympy as sym
from pauli_utils import are_subsets_equal
from sympy.physics.quantum import TensorProduct, Operator
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


# def find_agnostic_circuit(pauli_sum: PauliSum, target_paulis: list[Pauli], target_indexes: list[tuple[int, int]],
#                           iterations: int, compare_phases: bool = True) -> Circuit | None:
#     """
#     Find a circuit that maps pauli_sum to target_paulis with target_indexes. Gives the first circuit found that does
#     this, agnostic to any Pauli not in the target list.

#     Horribly wasteful approach - would be better to ask 'are subsets equal' rather than check for loads of possible
#     pauli sums. Still, fast enough for now.
#     """
#     n_trials = 100

#     current_target = pauli_sum.copy()  # the target starts as a copy of the initial pauli_sum with target paulis changed only
#     for i, indexes in enumerate(target_indexes):
#         current_target[indexes[0], indexes[1]] = target_paulis[i]
    
#     for max_iter in range(1, iterations):

#         mappable = check_mappable_via_clifford(pauli_sum, current_target)
#         if not mappable:
#             # change a pauli in the target - loop through the target pauli_strings and target qudits for change locations
#             print('Changing target state')
#             pass
#         else:
#             print('finding circuit with naive state')
#             C = find_circuit(pauli_sum, current_target, max_iter, compare_phases, True)
#             if C != []:
#                 return C[0]
        
#         paulis = [Pauli(0, 0, 2), Pauli(0, 1, 2), Pauli(1, 0, 2), Pauli(1, 1, 2)]  # only works for qubits for now
#         agnostic_pauli_locations = [(i, j) for i in range(pauli_sum.n_paulis()) for j in range(pauli_sum.n_qudits()) if (i, j) not in target_indexes]
#         # we make random changes to the target until we find one that works
#         for trial_number in range(n_trials):
#             ps, q = random.choice(agnostic_pauli_locations)  # pick a random pauli string and qudit to change
#             p = random.choice(paulis)  # pick a random pauli
#             current_target[ps, q] = p
#             mappable = check_mappable_via_clifford(pauli_sum, current_target)
#             if mappable:
#                 C = find_circuit(pauli_sum, current_target, max_iter, compare_phases, True)
#                 if C != []:
#                     print('Found on trial number ', trial_number + 1, ', depth = ', max_iter)
#                     return C[0]
                
#     return None


def find_agnostic_circuit(pauli_sum: PauliSum, target_paulis: list[Pauli], target_indexes: list[tuple[int, int]],
                          iterations: int, compare_phases: bool = True, stop_if_found: bool = True) -> list[Circuit]:

    n_qudits = pauli_sum.n_qudits()
    SUMs = [CX(i, j, 2) for i in range(n_qudits) for j in range(n_qudits) if i != j]
    Ss = [S(i, 2) for i in range(n_qudits)]
    Hs = [H(i, 2) for i in range(n_qudits)]
    all_gates = SUMs + Ss + Hs

    goal_circuits = []
    circuits = [Circuit(dimensions=pauli_sum.dimensions)]
    intermediate_paulis = [pauli_sum.copy()]
    for i in range(iterations):
        print('iteration ', i)
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
                if np.all([P_temp[target_indexes[i]] == target_paulis[i] for i in range(len(target_paulis))]):
                    goal_circuits.append(C_temp)
                    if stop_if_found:
                        return [C_temp]
    return goal_circuits


def symplectic_effect(circuit):
    n_qudits = len(circuit.dimensions)
    r_now = list(sym.symbols([f'r{i}' for i in range(1, n_qudits+1)]))
    omega = sym.symbols('omega')
    s_now = list(sym.symbols([f's{i}' for i in range(1, n_qudits+1)]))
    r_next = [r_now[i] for i in range(n_qudits)]
    s_next = [s_now[i] for i in range(n_qudits)]

    X = Operator('X')
    Z = Operator('Z')

    phase = 0
    gates = circuit.gates
    qubits = circuit.indexes
    for i, g in enumerate(gates):
        if g.name == 'SUM':
            r_next[qubits[i][1]] = r_now[qubits[i][0]] + r_now[qubits[i][1]]
            s_next[qubits[i][0]] = s_now[qubits[i][0]] + s_now[qubits[i][1]]
        elif g.name == 'H' or g.name == 'HADAMARD':
            r_next[qubits[i][0]] = s_now[qubits[i][0]]
            s_next[qubits[i][0]] = r_now[qubits[i][0]]
            phase += s_now[qubits[i][0]] * r_now[qubits[i][0]]
        elif g.name == 'S' or g.name == 'PHASE':
            s_next[qubits[i][0]] = s_now[qubits[i][0]] + r_now[qubits[i][0]]
            phase += r_now[qubits[i][0]] * (r_now[qubits[i][0]] - 1) / 2
        r_now = [r_next[i] for i in range(n_qudits)]
        s_now = [s_next[i] for i in range(n_qudits)]
    final = TensorProduct(X**(modulo_2(r_now[0])) * Z**(modulo_2(s_now[0])), X**(modulo_2(r_now[1])) * Z**(modulo_2(s_now[1])))
    for i in range(2, n_qudits):
        final = TensorProduct(final, X**(modulo_2(r_now[i])) * Z**(modulo_2(s_now[i])))
    
    display(omega**modulo_2(reduce_exponents(modulo_2(sym.simplify(phase)))) * final)


def modulo_2(expr):
    """
    Takes a SymPy expression and reduces its coefficients modulo 2.
    """
    # Expand the expression to handle all terms
    expr = expr.expand()
    
    # Iterate through the terms and apply modulo 2 to coefficients
    terms = expr.as_ordered_terms()
    mod_expr = sum(sym.Mod(term.as_coeff_Mul()[0], 2) * term.as_coeff_Mul()[1] for term in terms)
    
    return mod_expr


def reduce_exponents(expr):
    """
    Reduces all exponents in a SymPy expression to zero, assuming symbols are binary (0 or 1).
    
    Args:
        expr (sympy.Expr): The input SymPy expression.
    
    Returns:
        sympy.Expr: The modified expression with all exponents set to zero.
    """
    expr = sym.expand(expr)  # Expand the expression to handle all terms
    return expr.replace(lambda x: x.is_Pow, lambda x: x.base)


def random_gate(dimensions: list[int] | np.ndarray) -> GateOperation:
    gate = random.choice(['SUM', 'H', 'S'])
    qudits = np.arange(len(dimensions))
    
    if gate == 'SUM':
        control = random.choice(qudits)
        target = random.choice(np.delete(qudits, control))
        if dimensions[control] != dimensions[target]:  # reselect at random
            return random_gate(dimensions)
        return CX(control, target, dimensions[control])
    elif gate == 'H':
        qudit = random.choice(qudits)
        return H(qudit, dimensions[qudit])
    elif gate == 'S':
        qudit = random.choice(qudits)
        return S(qudit, dimensions[qudit])


def random_clifford(depth: int, dimensions: list[int] | np.ndarray) -> Circuit:
    circuit = Circuit(dimensions)
    for _ in range(depth):
        gate = random_gate(dimensions)
        circuit.add_gate(gate)
    return circuit


if __name__ == "__main__":
    from hamiltonian import random_pauli_hamiltonian
    initial_pauli = ['x0z0 x0z1', 'x1z1 x1z0']
    goal_pauli = ['x1z1 x1z1', 'x1z1 x0z1']
    initial_pauli = PauliSum(initial_pauli, dimensions=[2, 2])
    goal_pauli = PauliSum(goal_pauli, dimensions=[2, 2])
    print(initial_pauli)

    # C = find_circuit(initial_pauli, goal_pauli, 8)
    # print(C)
    initial_pauli = random_pauli_hamiltonian(6, [2, 2, 2, 2, 2, 2], mode='random')
    goal_single_pauli = Pauli(1, 0, 2)
    goal_location = (0, 1)

    C = find_agnostic_circuit(initial_pauli.copy(), [goal_single_pauli], [goal_location], 8)
    print(C[0])
    print(initial_pauli)
    print(C[0].act(initial_pauli))
