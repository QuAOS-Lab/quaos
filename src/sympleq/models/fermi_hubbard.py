# pip install openfermion
import numpy as np

from openfermion.transforms.opconversions.jordan_wigner import jordan_wigner
from openfermion.transforms.opconversions.bravyi_kitaev import bravyi_kitaev
from openfermion.hamiltonians.hubbard import fermi_hubbard

from sympleq.core.paulis import PauliSum


def qubit_pauli_tableau(qubit_op, n_qubits=None):
    """
    Encode each Pauli monomial as a row [x | z] with Y -> (1,1), plus its coeff.
    """
    if n_qubits is None:
        max_index = 0
        for term, _ in qubit_op.terms.items():
            if term:
                max_index = max(max_index, max(q for q, _ in term))
        n_qubits = max_index + 1

    rows, coeffs, labels = [], [], []
    for term, coeff in qubit_op.terms.items():
        x = np.zeros(n_qubits, dtype=int)
        z = np.zeros(n_qubits, dtype=int)
        parts = []
        for q, p in term:
            if p == 'X':
                x[q] = 1
            elif p == 'Y':
                x[q] = 1
                z[q] = 1
            elif p == 'Z':
                z[q] = 1
            else:
                raise ValueError(f"Unknown Pauli {p}")
            parts.append(f"{p}{q}")
        rows.append(np.concatenate([x, z]))
        coeffs.append(coeff)
        labels.append(" ".join(parts) if parts else "I")

    T = np.vstack(rows) if rows else np.zeros((0, 2 * n_qubits), dtype=int)
    return T, np.asarray(coeffs, dtype=complex), labels


def fermi_hubbard_tableau(Lx: int,
                          Ly: int = 1,
                          t: float = 1.0,
                          U: float = 4.0,
                          mu: float = 0.0,
                          periodic: bool = False,
                          spinless: bool = False,
                          mapping: str = "jw"):
    """
    Build the Fermi–Hubbard model on an Lx x Ly lattice, map to qubits,
    and return the tableau [x|z] with coefficients.

    Parameters
    ----------
    Lx, Ly : int
        Lattice dimensions (Ly=1 gives a chain).
    t : float
        Hopping amplitude.
    U : float
        On-site interaction strength.
    mu : float
        Chemical potential (number term).
    periodic : bool
        Periodic boundary conditions along both axes if True.
    spinless : bool
        Use spinless fermions if True (1 spin-orbital per site).
        Default is spinful (2 spin-orbitals per site).
    mapping : {"jw","bk"}
        Jordan-Wigner or Bravyi-Kitaev.
    save_csv : str | None
        If provided, write the tableau DataFrame to this path.

    Returns
    -------
    T : (n_terms, 2*n_qubits) int ndarray
    coeffs : (n_terms,) complex ndarray
    labels : list[str]
    n_qubits : int
    """
    # 1) Fermionic Hamiltonian
    fermion_ham = fermi_hubbard(
        x_dimension=Lx,
        y_dimension=Ly,
        tunneling=t,
        coulomb=U,
        chemical_potential=mu,
        periodic=periodic,
        spinless=spinless
    )

    # 2) Map to qubits
    mapper = jordan_wigner if mapping.lower() == "jw" else bravyi_kitaev
    qubit_op = mapper(fermion_ham)

    # 3) Qubit count (OpenFermion uses 2*sites for spinful; 1*sites for spinless)
    n_sites = Lx * Ly
    n_qubits = n_sites * (1 if spinless else 2)

    # 4) Tableau
    T, coeffs, labels = qubit_pauli_tableau(qubit_op, n_qubits=n_qubits)

    return T, coeffs, labels, n_qubits


def fermi_hubbard_model(x_dimension: int,
                        y_dimension: int = 1,
                        tunneling: float = 1.0,
                        coulomb: float = 4.0,
                        chemical_potential: float = 0.0,
                        periodic: bool = False,
                        spinless: bool = False):

    tableau, coeffs, _, n_qubits = fermi_hubbard_tableau(Lx=x_dimension,
                                                         Ly=y_dimension,
                                                         t=tunneling,
                                                         U=coulomb,
                                                         mu=chemical_potential,
                                                         periodic=periodic,
                                                         spinless=spinless)

    P = PauliSum.from_tableau(tableau, weights=coeffs, dimensions=[2] * n_qubits)
    return P


if __name__ == "__main__":
    fh_model = fermi_hubbard_model(x_dimension=4, y_dimension=2)
    print(fh_model)
    print(fh_model.n_qudits())
    print(fh_model.n_paulis())
