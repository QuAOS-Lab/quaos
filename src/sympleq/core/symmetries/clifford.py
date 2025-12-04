import numpy as np
from sympleq.core.circuits import Gate, Circuit
from sympleq.core.paulis import PauliSum
from sympleq.core.graphs.graph_isomorphism import _build_base_partition, _full_dfs_complete, _labels_union
from sympleq.core.finite_field_solvers import get_linear_dependencies
from scripts.experiments.symmetries.src.block_decomposition import block_decompose, ordered_block_sizes
from typing import Optional


def min_qudit_clifford_symmetry(pauli_sum: PauliSum,
                                check_symmetry: bool = True,
                                ) -> tuple[Gate, Gate, Gate]:

    G = find_clifford_symmetries(pauli_sum, num_symmetries=1,
                                 dynamic_refine_every=0)
    g = G[0]
    if check_symmetry:
        lhs = g.act(pauli_sum).to_standard_form()
        rhs = pauli_sum.to_standard_form()
        assert lhs == rhs, f'Symmetry finder failed\n{lhs.__str__()}\n{rhs.__str__()}'

    S, T = block_decompose(g.symplectic, int(pauli_sum.lcm), min_block_size=4)
    h_S, h_T = clifford_phase_decomposition(g.symplectic, g.phase_vector, S, T, int(pauli_sum.lcm))
    S_gate = Gate('S', g.qudit_indices, S, g.dimensions, h_S)
    T_gate = Gate('T', g.qudit_indices, T, g.dimensions, h_T)

    if check_symmetry:
        lhs = T_gate.act(S_gate.act(T_gate.inv().act(pauli_sum))).to_standard_form()
        rhs = pauli_sum.to_standard_form()
        ts = T_gate.symplectic
        tis = T_gate.inv().symplectic
        ss = S_gate.symplectic
        fs = g.symplectic
        C = Circuit(pauli_sum.dimensions, [T_gate.inv(), S_gate, T_gate])
        F = C.composite_gate()
        assert F.act(pauli_sum).to_standard_form() == pauli_sum.to_standard_form()

        assert np.all(F.symplectic == g.symplectic), f'symplectic failed\n{F.symplectic}\n{g.symplectic}'
        assert np.all(F.phase_vector == g.phase_vector), f'phase vector failed\n{F.phase_vector}\n{g.phase_vector}'

        assert np.array_equal(fs, (ts @ ss @ tis) % 2), f'symplectic failed\n{fs}\n{ts @ ss @ ts.T}'
        assert lhs == rhs, f'Localisation failed\n{lhs.__str__()}\n{rhs.__str__()}'
    return g, S_gate, T_gate


def multiple_clifford_symmetries(pauli_sum: PauliSum,
                                 n_symmetries: int = 1,
                                 check_symmetry: bool = True,
                                 ) -> tuple[list[Gate], list[Gate], list[Gate]]:

    G = find_clifford_symmetries(pauli_sum, num_symmetries=n_symmetries,
                                 dynamic_refine_every=0)

    if check_symmetry:
        for g in G:
            assert g.act(pauli_sum).to_standard_form() == pauli_sum.to_standard_form()

    Ss = []
    Ts = []
    for i, g in enumerate(G):
        S, T = block_decompose(g.symplectic, pauli_sum.lcm)
        h_S, h_T = clifford_phase_decomposition(g.symplectic, g.phase_vector, S, T, int(pauli_sum.lcm))
        S_gate = Gate(f'S{i}', g.qudit_indices, S, g.dimensions, h_S)
        T_gate = Gate(f'T{i}', g.qudit_indices, T, g.dimensions, h_T)
        Ss.append(S_gate)
        Ts.append(T_gate)

    return G, Ss, Ts


def block_structure(gate: Gate):
    symp = gate.symplectic
    sizes = np.asarray(ordered_block_sizes(symp, int(gate.lcm)), dtype=int) / 2
    return sizes


def qudit_cost(gate: Gate):
    return int(max(block_structure(gate)))


def clifford_phase_decomposition(F: np.ndarray, h_F: np.ndarray,
                                 S: np.ndarray, T: np.ndarray, d: int,
                                 l_T: np.ndarray | None = None):
    """
    Inputs:
      F,h_F : target Clifford (symplectic F, phase vector h_F) with phases mod 2d
      S,T   : symplectics satisfying F = T S T^{-1}
      d     : qudit dimension
      l_T   : optional gauge vector added to h_T (default 0)

    Outputs:
      h_S, h_T : phase vectors of S and T (mod 2d)

    Conventions:
      - Pauli exponent rows update as a' = a @ F.T.
    """
    mod = 2 * d
    n2 = F.shape[0]
    n = n2 // 2
    dims = [d] * n

    def diag_U(C: np.ndarray) -> np.ndarray:
        """Return diag(C^T U C) in the same convention used by Gate.act."""
        U = np.zeros((n2, n2), dtype=int)
        U[n:, :n] = np.eye(n, dtype=int)
        return np.diag(C.T @ U @ C) % mod

    # Gauge choice for T: default ℓ_T = 0  ⇒  h_T = diag(T^T U T) + ℓ_T
    if l_T is None:
        l_T = np.zeros(n2, dtype=int)
    h_T = (diag_U(T) + l_T) % mod

    # Helper to build the composite phase for given h_S, with h_T fixed above
    def composite_phase(h_S: np.ndarray) -> np.ndarray:
        T_gate = Gate('T', list(range(n)), T, dims, h_T)
        S_gate = Gate('S', list(range(n)), S, dims, h_S)
        # Order matters: [T^{-1}, S, T] gives composite symplectic T S T^{-1}
        circuit = Circuit(dims, [T_gate.inv(), S_gate, T_gate])
        return circuit.composite_gate().phase_vector % mod

    base = composite_phase(np.zeros(n2, dtype=int))

    # Build the linear map from h_S to the resulting phase vector
    cols = []
    eye = np.eye(n2, dtype=int)
    for i in range(n2):
        cols.append((composite_phase(eye[i]) - base) % mod)
    A = np.stack(cols, axis=1) % mod
    b = (h_F % mod - base) % mod

    def solve_mod(mat: np.ndarray, vec: np.ndarray, modulus: int) -> np.ndarray:
        """Gaussian elimination over Z_mod, requiring unit pivots (gcd=1)."""
        mat = mat.copy().astype(int)
        vec = vec.copy().astype(int)
        m, ncols = mat.shape
        aug = np.concatenate([mat, vec.reshape(-1, 1)], axis=1) % modulus
        row = 0
        for col in range(ncols):
            if row >= m:
                break
            pivot = None
            for r in range(row, m):
                if np.gcd(int(aug[r, col]), modulus) == 1:
                    pivot = r
                    break
            if pivot is None:
                continue
            if pivot != row:
                aug[[row, pivot]] = aug[[pivot, row]]
            inv = pow(int(aug[row, col]) % modulus, -1, modulus)
            aug[row] = (aug[row] * inv) % modulus
            for r in range(m):
                if r == row:
                    continue
                factor = aug[r, col]
                aug[r] = (aug[r] - factor * aug[row]) % modulus
            row += 1
        if row < ncols:
            # Under-determined or singular; fall back to least filled solution if possible
            pass
        return aug[:ncols, -1] % modulus

    h_S = solve_mod(A, b, mod)
    return h_S.astype(int), h_T.astype(int)


def find_clifford_symmetries(
    pauli_sum: PauliSum,
    num_symmetries: int = 1,
    # Strategy
    dynamic_refine_every: int = 0,
    extra_column_invariants: str = "none",
    p2_bitset: str = "auto",
    F_known_debug: Optional[np.ndarray] = None,
    color_mode: str = "wl",   # "wl" | "coeffs_only" | "none"
    max_wl_rounds: int = 10,
) -> list[Gate]:
    """
    Return up to k automorphisms preserving S and the vector set. See flags above.
    """
    independent, dependencies = get_linear_dependencies(pauli_sum.tableau, 2)
    S = pauli_sum.symplectic_product_matrix()
    G, basis_order = pauli_sum.matroid()
    coeffs = pauli_sum.weights

    if not np.all([pauli_sum.dimensions[i] == pauli_sum.dimensions[0] for i in range(1, len(pauli_sum.dimensions))]):
        raise ValueError("All qubits must have same dimension for now. The key things to fix are: "
                         "_gf_solve_one_solution, and the symplectic_solver for F.")
    p = int(pauli_sum.lcm)

    pres_labels = _labels_union(independent, dependencies)
    n = len(pres_labels)

    col_invariants = None
    if extra_column_invariants != "none":
        G_for_inv = G.copy()
        if extra_column_invariants == "hist":
            inv = np.zeros((n, min(p, 16)), dtype=np.int64)
            for j in range(n):
                col = np.array([int(x) for x in G_for_inv[:, j]])
                cnt = np.bincount(col, minlength=p)
                inv[j, :min(p, 16)] = cnt[:min(p, 16)]
            col_invariants = inv
        else:
            raise ValueError("extra_column_invariants must be 'none' or 'hist'.")

    # ---------- choose the base partition via color_mode ----------
    base_colors, base_classes = _build_base_partition(
        S, p,
        coeffs=coeffs,
        col_invariants=col_invariants if color_mode == "wl" else None,
        max_rounds=max_wl_rounds,
        color_mode=color_mode,
    )

    use_bitset = (p == 2 and (p2_bitset is True or (p2_bitset == "auto" and n <= 256)))

    return _full_dfs_complete(
        pauli_sum,
        independent,
        S,
        coeffs=coeffs,
        base_colors=base_colors,
        base_classes=base_classes,
        G=G, basis_order=basis_order, labels=pres_labels,
        k_wanted=num_symmetries,
        p2_bitset=use_bitset,
        dynamic_refine_every=int(dynamic_refine_every),
        F_known_debug=F_known_debug,
    )
