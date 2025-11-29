# pip install openfermion openfermionpyscf pyscf pandas
import numpy as np
import pandas as pd

# from openfermion.chem.molecular_data import MolecularData
# from openfermion.transforms.opconversions.jordan_wigner import jordan_wigner
# from openfermion.transforms.opconversions.conversions import get_fermion_operator
# from openfermionpyscf import run_pyscf


def qubit_pauli_tableau(qubit_op, n_qubits=None):
    """
    Convert an OpenFermion QubitOperator into a tableau [X|Z] and coefficients.

    Parameters
    ----------
    qubit_op : openfermion.ops.QubitOperator
        Sum of Pauli terms.
    n_qubits : int | None
        Total number of qubits; if None, inferred from operator.

    Returns
    -------
    tableau : np.ndarray, shape (n_terms, 2*n_qubits), dtype=int
        Binary tableau rows per term.
    coeffs : np.ndarray, shape (n_terms,), dtype=complex
        Coefficient per row.
    terms  : list[str]
        Human-readable Pauli strings (e.g., 'X0 Y2 Z5').
    """
    # Infer n_qubits
    if n_qubits is None:
        max_index = 0
        for term, _ in qubit_op.terms.items():
            if term:
                max_index = max(max_index, max(q for q, _ in term))
        n_qubits = max_index + 1

    rows = []
    coeffs = []

    for term, coeff in qubit_op.terms.items():
        x = np.zeros(n_qubits, dtype=int)
        z = np.zeros(n_qubits, dtype=int)

        # term is a tuple like ((q, 'X'), (r, 'Y'), ...)
        # Empty term () corresponds to identity
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

        rows.append(np.concatenate([x, z]))
        coeffs.append(coeff)

    return np.vstack(rows), np.array(coeffs, dtype=complex)


def tableau_dataframe(tableau, coeffs):
    n_qubits = tableau.shape[1] // 2
    cols = [f"X{j}" for j in range(n_qubits)] + [f"Z{j}" for j in range(n_qubits)] + ["coeff_re", "coeff_im"]
    df = pd.DataFrame(
        np.hstack([tableau, np.column_stack([coeffs.real, coeffs.imag])]),
        columns=cols
    )
    df[cols] = df[cols].astype(object)  # keep ints where possible
    return df


# def chemistry_hamiltonian(geometry, basis="sto-3g", multiplicity=1, charge=0,
#                           run_scf=True, run_mp2=False, run_cisd=False, run_ccsd=False, run_fci=False):
#     """
#     Generate a molecular Hamiltonian as a PauliSum using OpenFermion and PySCF.

#     Parameters
#     ----------
#     geometry : list of tuples
#         List of (element_symbol, (x,y,z)) tuples defining the molecule.
#     basis : str
#         Basis set to use (default "sto-3g").
#     multiplicity : int
#         Spin multiplicity (default 1 for singlet).
#     charge : int
#         Molecular charge (default 0).
#     run_scf, run_mp2, run_cisd, run_ccsd, run_fci : bool
#         Whether to run these calculations in PySCF.

#     Returns
#     -------
#     df : pd.DataFrame
#         DataFrame with columns ['Pauli_term', 'X0', ..., 'Z0', ..., 'coeff_re', 'coeff_im'].
#     n_qubits : int
#         Number of qubits in the Hamiltonian.
#     n_terms : int
#         Number of Pauli terms in the Hamiltonian.
#     tableau : np.ndarray
#         Binary tableau representation of the Hamiltonian.
#     """
#     molecule = MolecularData(geometry, basis, multiplicity, charge)

#     # Run PySCF to get molecular integrals and build the Fermionic Hamiltonian
#     molecule = run_pyscf(molecule, run_scf=run_scf, run_mp2=run_mp2,
#                          run_cisd=run_cisd, run_ccsd=run_ccsd, run_fci=run_fci)
#     fermion_ham = molecule.get_molecular_hamiltonian()
#     fermion_op = get_fermion_operator(fermion_ham)

#     # Map to qubits (choose jordan_wigner or bravyi_kitaev from openfermion)
#     qubit_op = jordan_wigner(fermion_op)

#     # Convert to tableau
#     T, coeffs = qubit_pauli_tableau(qubit_op)

#     # Present as a dataframe (helpful for inspection/CS
#     n_qubits = T.shape[1] // 2

#     h = PauliSum.from_tableau(T, [2] * n_qubits, coeffs)

#     return h


# def water_molecule(basis='sto-3g'):
#     """
#     Water (angstrom), angle 104.5°, OH bond length 0.9584.
#     Basis: STO-3G.
#     1086 Paulis.
#     14 qubits.
#     Basis: 6-31G
#     12732 Paulis.
#     26 qubits
#     """
#     geometry = [("O", (0.0, 0.0, 0.0)),
#                 ("H", (0.9584, 0.0, 0.0)),
#                 ("H", (0.9584 * np.cos(np.radians(104.5)), 0.9584 * np.sin(np.radians(104.5)), 0.0))]
#     multiplicity = 1
#     charge = 0

#     h = chemistry_hamiltonian(geometry, basis, multiplicity, charge,
#                               run_scf=True, run_mp2=False, run_cisd=False,
#                               run_ccsd=False, run_fci=False)

#     return h


# def beh2_molecule(basis='cc-pVDZ'):
#     """
#     Linear BeH2 (angstrom), Be at origin, H along ±z.
#     Basis: cc-pVDZ.
#     102697 paulis
#     48 qubits.
#     """
#     b_be_h = 1.316  # angstrom
#     geometry = [
#         ("Be", (0.0, 0.0, 0.0)),
#         ("H", (0.0, 0.0, +b_be_h)),
#         ("H", (0.0, 0.0, -b_be_h)),
#     ]
#     multiplicity = 1
#     charge = 0

#     h = chemistry_hamiltonian(
#         geometry, basis, multiplicity, charge,
#         run_scf=True, run_mp2=False, run_cisd=False,
#         run_ccsd=False, run_fci=False
#     )
#     return h


# def ch4_molecule(basis='6-31G'):
#     """
#     Methane (angstrom), ideal tetrahedral CH bond length.
#     Basis: 6-31G.
#     103652 paulis
#     34 qubits.
#     """
#     b_ch = 1.09  # angstrom
#     a = b_ch / np.sqrt(3.0)
#     geometry = [
#         ("C", (0.0, 0.0, 0.0)),
#         ("H", (+a, +a, +a)),
#         ("H", (+a, -a, -a)),
#         ("H", (-a, +a, -a)),
#         ("H", (-a, -a, +a)),
#     ]
#     multiplicity = 1
#     charge = 0

#     h = chemistry_hamiltonian(
#         geometry, basis, multiplicity, charge,
#         run_scf=True, run_mp2=False, run_cisd=False,
#         run_ccsd=False, run_fci=False
#     )
#     return h


# def c2h4_molecule(basis='6-31G'):
#     """
#     Ethylene (angstrom), planar C=C along x, H's ~120° in-plane.
#     Basis: 6-31G.
#     104263 paulis
#     52 qubits.
#     """
#     d_cc = 1.339  # angstrom (C=C)
#     d_ch = 1.09   # angstrom
#     cx = d_cc / 2.0

#     # angles relative to +x and -x directions
#     ang = np.deg2rad(60.0)
#     # Left carbon at (-cx, 0, 0)
#     vL1 = (np.cos(np.pi - ang), np.sin(np.pi - ang), 0.0)
#     vL2 = (np.cos(np.pi + ang), np.sin(np.pi + ang), 0.0)
#     # Right carbon at (+cx, 0, 0)
#     vR1 = (np.cos(+ang), np.sin(+ang), 0.0)
#     vR2 = (np.cos(-ang), np.sin(-ang), 0.0)

#     geometry = [
#         ("C", (-cx, 0.0, 0.0)),
#         ("C", (+cx, 0.0, 0.0)),
#         ("H", (-cx + d_ch * vL1[0], 0.0 + d_ch * vL1[1], 0.0)),
#         ("H", (-cx + d_ch * vL2[0], 0.0 + d_ch * vL2[1], 0.0)),
#         ("H", (+cx + d_ch * vR1[0], 0.0 + d_ch * vR1[1], 0.0)),
#         ("H", (+cx + d_ch * vR2[0], 0.0 + d_ch * vR2[1], 0.0)),
#     ]
#     multiplicity = 1
#     charge = 0

#     h = chemistry_hamiltonian(
#         geometry, basis, multiplicity, charge,
#         run_scf=True, run_mp2=False, run_cisd=False,
#         run_ccsd=False, run_fci=False
#     )
#     return h


# def c6h6_molecule(basis='sto-3g'):
#     """
#     Benzene (angstrom), flat hexagon in xy-plane; H's radial.
#     Basis: STO-3G.
#     500909 paulis
#     72 qubits.
#     """
#     d_cc = 1.397  # angstrom (approx. ring C–C)
#     d_ch = 1.09   # angstrom (C–H)
#     geometry = []

#     for k in range(6):
#         theta = 2.0 * np.pi * k / 6.0
#         cx = d_cc * np.cos(theta)
#         cy = d_cc * np.sin(theta)
#         cz = 0.0
#         # unit radial vector from center
#         r_norm = np.sqrt(cx**2 + cy**2)
#         rx, ry = cx / r_norm, cy / r_norm
#         hx = cx + d_ch * rx
#         hy = cy + d_ch * ry
#         hz = 0.0
#         geometry.append(("C", (cx, cy, cz)))
#         geometry.append(("H", (hx, hy, hz)))

#     multiplicity = 1
#     charge = 0

#     h = chemistry_hamiltonian(
#         geometry, basis, multiplicity, charge,
#         run_scf=True, run_mp2=False, run_cisd=False,
#         run_ccsd=False, run_fci=False
#     )
#     return h


# if __name__ == "__main__":
#     print('Water')
#     h = water_molecule(basis='6-31G')
#     # print(h)
#     print(h.n_paulis())
#     print(h.n_qudits())

#     # print('BeH2')
#     # h = beh2_molecule()
#     # print(h.n_paulis())
#     # print(h.n_qudits())

#     # print('CH4')
#     # h = ch4_molecule()
#     # print(h.n_paulis())
#     # print(h.n_qudits())

#     # print('C2H4')
#     # h = c2h4_molecule()
#     # print(h.n_paulis())
#     # print(h.n_qudits())

#     # print('C6H6')
#     # h = c6h6_molecule()
#     # print(h.n_paulis())
#     # print(h.n_qudits())
