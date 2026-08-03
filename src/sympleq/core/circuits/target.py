"""Codes for finding target Paulis and gates which map a given Pauli to a target Pauli."""
# TODO: move all functions that are required for the class method
#       `input_to_target` to a unique file, that works with mixed dimension.
from __future__ import annotations
from sympleq.core.paulis import PauliSum
from sympleq.core.paulis._typing import TableauType, PhasesType
import numpy as np
from sympleq.core.circuits.helpers_solve_from_target import get_phase_vector, map_tableau_to_target_tableau

# TODO: Remove this file across repo and replace with solve_from_target
def find_map_to_target_pauli_sum(input_pauli: PauliSum, target_pauli: PauliSum) -> tuple[TableauType, PhasesType,
                                                                                         list[int], int]:
    """
    TODO: For efficiency improvement act only on target qudits

    Find a gate that maps Pauli P to target Pauli.

    Args:
        P (Pauli): The Pauli to be mapped.
        target (Pauli): The target Pauli.
        dimension (int): The dimension of the qudit.

    Returns:
        images (list[np.ndarray]): The images of the gate.
        h (np.ndarray): The phase vector of the gate.
        qudit_indices (list[int]): The indices of the qudits acted upon by the gate.
        gate_dimension (int): The dimension of the gate.
        """
    if np.all(input_pauli.dimensions != target_pauli.dimensions):
        raise ValueError("PauliSum and gate must have the same dimension.")

    n_qudits = input_pauli.n_qudits()
    if n_qudits != target_pauli.n_qudits():
        raise ValueError("PauliSum and target must have the same number of qudits.")

    # get list of qudits where input and target differ
    qudit_indices = list(range(n_qudits))
    gate_dimension = input_pauli.dimensions[qudit_indices[0]]

    if not np.all(input_pauli.dimensions[qudit_indices] == gate_dimension):
        raise ValueError("PauliSum must have the same dimension for all qudits acted upon by the gate.")

    if np.all(input_pauli.symplectic_product_matrix() != target_pauli.symplectic_product_matrix()):
        raise ValueError("Input and target PauliSum must be symplectically equivalent.")

    input_tableau = input_pauli.tableau  # [:, qudit_indices]
    target_tableau = target_pauli.tableau  # [:, qudit_indices]

    F = map_tableau_to_target_tableau(input_tableau, target_tableau, p=int(gate_dimension))
    if F is None:
        raise ValueError("No symplectic map found for the target PauliSum.")

    h = get_phase_vector(F, gate_dimension)

    return F, h, qudit_indices, gate_dimension
