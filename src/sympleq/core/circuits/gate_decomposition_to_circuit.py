import numpy as np
from sympleq.core.circuits import Circuit
from sympleq.core.circuits.circuits import GateSpec
from sympleq.core.circuits.utils import symplectic_form, is_symplectic
from sympleq.core.circuits.gates import GATES, Gate, PauliGate
from sympleq.core.paulis import PauliString
from collections import deque

# Convenience aliases for singleton gates
H = GATES.H
H_inv = GATES.H_inv
S = GATES.S
S_inv = GATES.S_inv
CX = GATES.CX
SWAP = GATES.SWAP

# ---------- GF(p) helpers - should remove most and put in better place ----------


def mod_p(A, p):
    return np.asarray(A, dtype=int) % p


def identity(n):
    return np.eye(n, dtype=int)


def zeros(shape):
    return np.zeros(shape, dtype=int)


def mm_p(A, B, p):
    return mod_p(A @ B, p)


def inv_gfp(A: np.ndarray, p: int) -> np.ndarray:
    """Inverse over GF(p) by Gauss–Jordan."""
    A = mod_p(A.copy(), p)
    n = A.shape[0]
    if A.shape[0] != A.shape[1]:
        raise ValueError("Matrix must be square.")
    Id = identity(n)
    aug = np.hstack((A, Id))
    row = 0
    for col in range(n):
        piv = None
        for r in range(row, n):
            if aug[r, col] % p != 0:
                piv = r
                break
        if piv is None:
            raise ValueError("Singular over GF(p).")
        if piv != row:
            aug[[row, piv]] = aug[[piv, row]]
        inv_pivot = pow(int(aug[row, col] % p), -1, p)
        aug[row, :] = (aug[row, :] * inv_pivot) % p
        for r in range(n):
            if r == row:
                continue
            f = aug[r, col] % p
            if f != 0:
                aug[r, :] = (aug[r, :] - f * aug[row, :]) % p
        row += 1
    return aug[:, n:] % p


def blocks(F: np.ndarray, p: int):
    F = mod_p(F, p)
    n = F.shape[0] // 2
    return F[:n, :n], F[:n, n:], F[n:, :n], F[n:, n:]


def _is_invertible_mod(A: np.ndarray, p: int) -> bool:
    try:
        _ = inv_gfp(A, p)
        return True
    except ValueError:
        return False


# ---------- preconditioner (BFS over Gate objects) ----------


def ensure_invertible_A_circuit(F: np.ndarray, p: int, max_depth: int | None = None) -> Circuit:
    """
    Find a Circuit C_pre such that the A-block of (C_pre.full_symplectic(n) @ F) is invertible mod p.
    Returns the Circuit (does NOT modify F).
    """
    F = mod_p(F, p)
    assert is_symplectic(F, p), "F must be symplectic over GF(p)."
    n = F.shape[0] // 2

    A, _, _, _ = blocks(F, p)
    if _is_invertible_mod(A, p):
        return Circuit.empty([p] * n)

    if max_depth is None:
        max_depth = max(1, 3 * n)

    # Candidate generators as (gate, qudits...) tuples
    candidates: list[GateSpec] = []
    for i in range(n):
        candidates.append((H, i))
        candidates.append((S, i))
        candidates.append((S_inv, i))
    for i in range(n):
        for j in range(i + 1, n):
            candidates.append((SWAP, i, j))
            # Both directions of CX are useful during search
            candidates.append((CX, i, j))
            candidates.append((CX, j, i))

    def apply_left_once(F_mat: np.ndarray, op: GateSpec) -> np.ndarray:
        gate, *qudits = op
        return (gate.full_symplectic(tuple(qudits), n, p) @ F_mat) % p

    # BFS with an index (no deque) for clarity and determinism
    seen = {tuple(F.flatten())}
    queue: list[tuple[np.ndarray, list[GateSpec]]] = [(F, [])]
    head = 0

    while head < len(queue):
        F_cur, ops = queue[head]
        head += 1
        if len(ops) >= max_depth:
            continue
        for g in candidates:
            F_new = apply_left_once(F_cur, g)
            key = tuple(F_new.flatten())
            if key in seen:
                continue
            new_ops = ops + [g]

            Anew, _, _, _ = blocks(F_new, p)
            if _is_invertible_mod(Anew, p):
                # Build the concrete Circuit
                C_pre = Circuit.empty([p] * n)
                for gate, *qudits in new_ops:
                    C_pre.add_gate(gate, *qudits)
                return C_pre

            seen.add(key)
            queue.append((F_new, new_ops))

    raise RuntimeError("Unable to make A invertible within the depth limit.")


#  ############################################################################  #
#                                  M - Block                                     #
#  ############################################################################  #

# ---- synthesize local diagonal D(u) = [[u,0],[0,u^{-1}]] on wire `index` ----
# Cache stores gate names ('S', 'H', 'H_inv') rather than gate instances
_DU_CACHE: dict[int, dict[int, list[str]]] = {}   # per p: u -> word (list of gate names on wire 0)


def _compose_1q_word_as_matrix(p: int, word0: list[str]) -> np.ndarray:
    """Compose a 1-qudit word (on qudit index 0) into its 2x2 symplectic over GF(p)."""
    gate_map = {'S': S, 'H': H, 'H_inv': H_inv}
    gates = [gate_map[name] for name in word0]
    qudits = [(0,)] * len(gates)
    return Circuit.from_gates_and_qudits([p], gates, qudits).composite_gate().symplectic % p


def _all_Du_targets_1q(p: int) -> dict[int, np.ndarray]:
    targets = {}
    for u in range(2, p):
        inv_u = pow(u, -1, p)
        targets[u] = np.array([[u, 0], [0, inv_u]], dtype=int) % p
    return targets


def _build_Du_cache_for_p(p: int, bfs_depth: int = 12) -> None:
    if p in _DU_CACHE:
        return

    # Generators as gate name lists
    gens = [['S'], ['H'], ['H_inv']]

    I2 = np.eye(2, dtype=int) % p
    targets = _all_Du_targets_1q(p)  # now excludes u==1
    remaining = set(targets.keys())

    solutions: dict[int, list[str]] = {}
    # Pre-seed identity: D(1) = I has empty word
    solutions[1] = []
    _DU_CACHE[p] = solutions  # write early so recursion-safe
    # BFS only for the nontrivial diagonals
    from collections import deque
    Q = deque([([], I2)])
    seen = {tuple(I2.flatten())}

    while Q and remaining:
        word, M = Q.popleft()
        if len(word) > bfs_depth:
            continue
        for g in gens:
            w2 = word + g
            M2 = _compose_1q_word_as_matrix(p, w2)
            key = tuple(M2.flatten())
            if key in seen:
                continue
            seen.add(key)
            # check if M2 matches any remaining D(u)
            hit = [u for u in list(remaining) if np.array_equal(M2, targets[u])]
            for u in hit:
                solutions[u] = w2
                remaining.remove(u)
            Q.append((w2, M2))

    if remaining:
        raise RuntimeError(f"1q BFS couldn't reach D(u) for u={sorted(remaining)} over GF({p}) "
                           f"with depth={bfs_depth}; increase depth or add generators.")

    _DU_CACHE[p] = solutions  # finalize


def _emit_local_ops_for_D(index: int, u: int, p: int) -> list[GateSpec]:
    """
    Emit a minimal word for D(u) on the given wire `index`, using a cached BFS word on wire 0.
    """
    u %= p
    if u == 1 % p:
        return []
    _build_Du_cache_for_p(p, bfs_depth=12)  # build once per p
    word0 = _DU_CACHE[p][u]                 # gate names on wire 0
    # remap the cached wire-0 gates to the requested wire `index`
    gate_map = {'S': S, 'H': H, 'H_inv': H_inv}
    remapped: list[GateSpec] = []
    for name in word0:
        remapped.append((gate_map[name], index))
    return remapped


# ---- MAIN: synthesize M = diag(A, (A^T)^{-1}) as a gate list (right-multiplication) ----
def synth_linear_A_to_gates(n: int, A: np.ndarray, p: int) -> list[GateSpec]:
    """
    RIGHT-multiplication consistent synthesis of [[A,0],[0,(A^T)^{-1}]].

    We eliminate on A^T with LEFT row ops, mirroring each row-op L as a RIGHT column-op E on A:
      • Row swap on A^T (r <-> c)         -> SWAP(r, c)
      • Scale row c by λ on A^T           -> scale column c by λ^{-1} on the RIGHT
                                            BUT since we normalized A^T by inv(λ), we must apply D(λ) on A.  # see below
      • Row add: row i ← row i − f*row c  -> col i ← col i − f*col c (RIGHT)
                                            -> use CX(src=c, dst=i) repeated (-f) mod p times.
    """
    A = (A % p).copy()
    At = A.T.copy()
    ops: list[GateSpec] = []

    for c in range(n):
        # --- pivot selection (prefer 1 to avoid scaling) ---
        pivot_row = None
        for r in range(c, n):
            if At[r, c] % p == 1 % p:
                pivot_row = r
                break
        if pivot_row is None:
            for r in range(c, n):
                if At[r, c] % p != 0:
                    pivot_row = r
                    break
        if pivot_row is None:
            for r in range(0, c):
                if At[r, c] % p != 0:
                    pivot_row = r
                    break
        if pivot_row is None:
            raise RuntimeError("A invertible, but no pivot found (unexpected).")

        # Row swap on A^T -> column swap on A
        if pivot_row != c:
            At[[c, pivot_row], :] = At[[pivot_row, c], :]
            ops.append((SWAP, pivot_row, c))

        # Scale row c of A^T to 1: multiply row c by inv(pv)
        # RIGHT effect: scale column c of A by pv (not inv_pv)  # <<< FIX
        pv = int(At[c, c] % p)
        if pv != 1 % p:
            inv_pv = pow(pv, -1, p)
            At[c, :] = (At[c, :] * inv_pv) % p
            ops += _emit_local_ops_for_D(c, pv, p)             # <<< FIX: use pv, not inv_pv

        # Eliminate other rows in column c on A^T:
        # row i <- row i - f * row c  =>  col i <- col i - f * col c (RIGHT)
        for i in range(n):
            if i == c:
                continue
            f = int(At[i, c] % p)
            if f != 0:
                At[i, :] = (At[i, :] - f * At[c, :]) % p
                reps = (f) % p
                for _ in range(reps):
                    ops.append((CX, i, c))

    return ops

#  ############################################################################  #
#                        Symmetric L and R block synthesis                       #
#  ############################################################################  #


def _as_int_mod(A, p: int) -> np.ndarray:
    """Coerce A to a plain integer ndarray and reduce mod p.
    Accepts: ndarray-like, Circuit/Gate with .symplectic, lists of lists, etc.
    """
    # unwrap symplectic matrices if a Circuit/Gate-ish object was passed
    if hasattr(A, "symplectic"):
        A = A.symplectic
    # standardize to ndarray first
    A = np.asarray(A)
    # aggressively coerce dtype to integer
    if A.dtype.kind in ("O", "U", "S"):  # object/str/unicode
        try:
            A = A.astype(np.int64, copy=False)
        except Exception:
            # last resort: element-wise int()
            A = np.vectorize(lambda x: int(x))(A)
    else:
        A = A.astype(np.int64, copy=False)
    return A % int(p)


def _full_from_lower(n, C, p):
    Id = np.eye(n, dtype=int)
    Z = np.zeros((n, n), dtype=int)
    F = np.block([[Id, Z], [_as_int_mod(C, p), Id]])
    return _as_int_mod(F, p)


def _compose_symp(n: int, p: int, ops: list[GateSpec]) -> np.ndarray:
    """Compose a list of (gate, qudits...) tuples into full symplectic matrix."""
    gates = [op[0] for op in ops]
    qudits = [op[1:] for op in ops]
    F = Circuit.from_gates_and_qudits([p] * n, gates, qudits).composite_gate().symplectic
    return _as_int_mod(F, p)


def _gens_2q(i: int, j: int) -> list[list[GateSpec]]:
    """Small generating set on two wires i,j (right-multiplication consistent)."""
    return [
        [(S, i)], [(S, j)],
        [(H, i)], [(H_inv, i)],
        [(H, j)], [(H_inv, j)],
        [(CX, i, j)], [(CX, j, i)],
        [(SWAP, i, j)],
    ]


def _bfs_2q_to_target(n: int, p: int, target_full: np.ndarray, i: int, j: int, max_depth: int = 7) -> list[GateSpec]:
    """Breadth-first synth on the two wires i,j; returns a word (list of (gate, qudits) tuples)."""
    Id = np.eye(2 * n, dtype=int) % p
    if np.array_equal(target_full % p, Id):
        return []
    gens = _gens_2q(i, j)
    Q = deque([([], Id)])
    seen = {tuple(Id.flatten())}
    while Q:
        word, M = Q.popleft()
        if len(word) >= max_depth:
            continue
        for g in gens:
            new_word = word + g
            M2 = _compose_symp(n, p, new_word)
            key = tuple(M2.flatten())
            if key in seen:
                continue
            if np.array_equal(M2, target_full % p):
                return new_word
            seen.add(key)
            Q.append((new_word, M2))
    raise RuntimeError("2q BFS failed (increase depth or check generators).")


def synth_lower_from_symmetric(n: int, C_sym: np.ndarray, p: int) -> list[GateSpec]:
    """
    Build L(C) = [[I,0],[C,I]] with right-multiplication convention.
    - Diagonals: C_ii * PHASE(i)
    - Off-diagonals: for each i<j with c != 0, BFS a minimal 2-qudit word that
      realizes [[I,0],[c(E_ij + E_ji), I]].
    """
    C_sym = mod_p(C_sym, p)
    assert np.array_equal(C_sym, C_sym.T % p), "C_sym must be symmetric mod p."

    ops: list[GateSpec] = []

    # Diagonal entries via local PHASE
    for i in range(n):
        t = int(C_sym[i, i] % p)
        if t:
            # PHASE acts like adding to the lower-left block on wire i
            ops += [(S, i)] * t

    # Off-diagonals via tiny BFS on the 2-wire subspace
    for i in range(n):
        for j in range(i + 1, n):
            c = int(C_sym[i, j] % p)
            if c == 0:
                continue
            # Build the target lower-symmetric elementary
            E = np.zeros((n, n), dtype=int)
            E[i, j] = (E[i, j] + c) % p
            E[j, i] = (E[j, i] + c) % p
            T = _full_from_lower(n, E, p)
            ops += _bfs_2q_to_target(n, p, T, i, j, max_depth=7)

    return ops


def synth_upper_from_symmetric_via_H(n: int, S_sym: np.ndarray, p: int) -> list[GateSpec]:
    """
    Build R(S) = [[I,S],[0,I]] using the H-sandwich:
      R(S) = (H_all) · L(-S) · (H_all)^(-1)   (right-multiplication)
    """
    S_sym = mod_p(S_sym, p)
    assert np.array_equal(S_sym, S_sym.T % p), "S_sym must be symmetric mod p."
    H_all: list[GateSpec] = [(H, q) for q in range(n)]
    H_all_inv: list[GateSpec] = [(H_inv, q) for q in range(n)]
    ops_L = synth_lower_from_symmetric(n, (-S_sym) % p, p)
    return H_all + ops_L + H_all_inv


# =========================
#  FINAL DECOMPOSITION
# =========================

def _add_ops_to_circuit(circuit: Circuit, ops: list[GateSpec]) -> None:
    """Helper to add a list of (gate, qudits...) tuples to a circuit."""
    for gate, *qudits in ops:
        circuit.add_gate(gate, *qudits)


def decompose_symplectic_to_circuit(F: np.ndarray, p: int, *, check: bool = True) -> Circuit:
    """
    Return a Circuit whose composite symplectic equals F (mod p),
    consistent with your RIGHT-multiplication convention for Pauli action.

      1) Find C_pre so that A' of F' := C_pre · F is invertible.
      2) Split F' into blocks A,B,C,D and form the symmetric Schur pieces:
           B' := A^{-1} B,   C' := C A^{-1}
         (both must be symmetric over GF(p) for symplectic F')
      3) Synthesize
           R ~ [[I, B'], [0, I]]
           M ~ [[A, 0], [0, (A^T)^{-1}]]
           L ~ [[I, 0], [C', I]]
      4) Since F' = L M R, we have
           F = C_pre^{-1} · L · M · R
         With your Circuit accumulation rule (left-multiply accumulator),
         the gate list must be [R, M, L, C_pre^{-1}].

    Returns
    -------
    Circuit  (gates ordered so that Circuit.composite_gate().symplectic == F mod p)
    """
    F = np.asarray(F, dtype=int) % p
    assert is_symplectic(F, p), "F must be symplectic over GF(p)."
    n = F.shape[0] // 2

    # 1) Preconditioner
    C_pre = ensure_invertible_A_circuit(F, p)
    F_prime = (C_pre.full_symplectic() @ F) % p
    A, B, C, D = blocks(F_prime, p)

    # 2) Symmetric Schur complements
    A_inv = inv_gfp(A, p)
    Bp = (A_inv @ B) % p
    Cp = (C @ A_inv) % p
    if check:
        assert np.array_equal(Bp % p, Bp.T % p), "A^{-1} B must be symmetric (bug)."
        assert np.array_equal(Cp % p, Cp.T % p), "C A^{-1} must be symmetric (bug)."

    # 3) Synthesize each block
    ops_R = synth_upper_from_symmetric_via_H(n, Bp, p)      # [[I, B'], [0, I]]
    ops_M = synth_linear_A_to_gates(n, A, p)                # [[A, 0], [0, (A^T)^{-1}]]
    ops_L = synth_lower_from_symmetric(n, Cp, p)            # [[I, 0], [C', I]]

    # 4) Assemble with your circuit accumulation rule (left-multiply stack),
    #    but aiming for F = C_pre^{-1} · L · M · R  => append [R, M, L, C_pre^{-1}]
    C_tot = Circuit.empty([p] * n)
    _add_ops_to_circuit(C_tot, ops_R)
    _add_ops_to_circuit(C_tot, ops_M)
    _add_ops_to_circuit(C_tot, ops_L)

    # Add inverse of preconditioner
    C_pre_inv = C_pre.inverse()
    for gate, qudits in zip(C_pre_inv.gates, C_pre_inv.qudit_indices):
        C_tot.add_gate(gate, *qudits)

    if check:
        F_rec = C_tot.composite_gate().symplectic % p
        if not np.array_equal(F_rec, F):
            # Minimal, targeted debug
            A_rec, B_rec, C_rec, D_rec = blocks(F_rec, p)
            A_tgt, B_tgt, C_tgt, D_tgt = blocks(F, p)
            raise AssertionError(
                "Decomposition round-trip failed.\n"
                f"Got:\n{F_rec}\nExp:\n{F}\n"
                f"[A got vs exp]\n{A_rec}\n{A_tgt}\n"
                f"[B got vs exp]\n{B_rec}\n{B_tgt}\n"
                f"[C got vs exp]\n{C_rec}\n{C_tgt}\n"
                f"[D got vs exp]\n{D_rec}\n{D_tgt}\n"
            )

    return C_tot


def solve_g_for_phase_delta(delta_h: np.ndarray, p: int) -> np.ndarray | None:
    """
    Solve 2 Ω g = delta_h (mod p) for prime p (odd only).
    Returns g or None if unsolvable (e.g., p=2 and delta_h not in image).
    """
    delta_h = np.asarray(delta_h, dtype=int) % p
    n2 = delta_h.size
    n = n2 // 2
    Ω = symplectic_form(n, p)
    # For odd p, 2 is invertible; for p=2 there is no solution unless delta_h ≡ 0
    if p % 2 == 0:
        return delta_h if np.all(delta_h % 2 == 0) else None  # informational, not actually g
    inv2 = pow(2, -1, p)
    # Ω is invertible over GF(p); Ω^{-1} = -Ω (since Ω^2 = -I) for odd p
    # You can compute g = (1/2) * Ω^{-1} delta_h; using generic inverse for safety:
    Ω_inv = inv_gfp(Ω, p)
    g = (inv2 * (Ω_inv @ delta_h)) % p
    return g


def pauli_gate_for_phase_fix(h_raw: np.ndarray,
                             h_target: np.ndarray,
                             dimensions: list[int] | np.ndarray) -> Gate | None:
    """
    Build a PauliGate that corrects (part of) delta_h = h_target - h_raw by Pauli conjugation.
    Uses per-wire dims via PauliString with your global 2*lcm phase ring.
    Returns None if there is no nontrivial solvable component (common for p=2).
    """
    dims = list(map(int, dimensions))
    n = len(dims)
    lcm = int(np.lcm.reduce(dims))
    # We’ll fix each prime-power component via CRT in practice; for now, do per-prime p where possible.
    # The simplest usable path in your pipeline (since you decompose per single p) is:
    p = dims[0]  # your code currently uses uniform p per run
    delta_h = (np.asarray(h_target, dtype=int) - np.asarray(h_raw, dtype=int)) % (2 * lcm)

    if p % 2 == 0:
        # Pauli conjugation cannot move phases mod 2; skip if nothing to do
        if np.all(delta_h % 2 == 0):
            return None
        # If you want: raise or log that a Clifford tweak (not Pauli) is needed
        return None

    # Reduce the desired shift into GF(p) and solve for g
    g = solve_g_for_phase_delta(delta_h % p, p)
    if g is None:
        return None

    x = (g[:n] % p).tolist()
    z = (g[n:] % p).tolist()
    pauli = PauliString.from_exponents(x_exp=x, z_exp=z, dimensions=dims)
    return PauliGate(pauli)


def pauli_correction_gate(F: np.ndarray,
                          h_body: np.ndarray,
                          h_target: np.ndarray,
                          dimensions: list[int]) -> "PauliGate":
    """
    Build a PauliGate that fixes the phase vector under right-multiplication:
        h_out = h_body + F^T * (2 Ω v)  (mod 2*lcm)
    Choose v in (Z_p)^{2n} so that h_out == h_target.

    For odd p: exact correction.
    For p=2: corrects the even part (the best Pauli can do).
    """
    import numpy as np

    dims = list(map(int, dimensions))
    assert len(set(dims)) == 1, ("Mixed dimension qudits not supported. As there are no entangling mixed qudit gates,"
                                 " split the system into qudits of equal dimension and run separately.")
    p = dims[0]

    n2 = F.shape[0]
    assert n2 % 2 == 0
    n = n2 // 2

    # modulus for the phase vector
    lcm = int(np.lcm.reduce(dims))
    MOD = 2 * lcm

    Δ = (h_target.astype(int) - h_body.astype(int)) % MOD  # even if realizable by Pauli
    if p % 2 == 1:
        # Must be even component-wise; then we can halve in Z and reduce mod p
        if np.any(Δ % 2 != 0):
            raise AssertionError("For odd p, Δ must be even for Pauli correction.")
        Δ_half = ((Δ // 2) % p).astype(int)  # well-defined because Δ is even

        # Ω over GF(p): [[0, I], [-I, 0]] (note -I ≡ p-1 mod p)
        Id = np.eye(n, dtype=int)
        Z = np.zeros((n, n), dtype=int)
        Om = np.block([[Z, Id], [(-Id) % p, Z]]) % p

        # Solve (F^T Ω) v = Δ_half  over GF(p)
        M = (F.T % p) @ Om % p
        v = (inv_gfp(M, p) @ Δ_half.reshape(-1, 1)) % p
        v = v.flatten()

        x = v[:n]
        z = v[n:]
        pauli = PauliString.from_exponents(x, z, dimensions=dims)
        return PauliGate(pauli)

    else:
        # p == 2: Ω has -I == +I; only even Δ are reachable (gives Δ/2 in {0,1})
        if np.any(Δ % 2 != 0):
            # Can't fix odd residuals with a Pauli; return identity Pauli (no-op)
            x = np.zeros(n, dtype=int)
            z = np.zeros(n, dtype=int)
            return PauliGate(PauliString.from_exponents(x, z, dimensions=dims))

        Δ_half = ((Δ // 2) % 2).astype(int)

        # Ω over GF(2): [[0, I], [ I, 0]]
        Id = np.eye(n, dtype=int)
        Z = np.zeros((n, n), dtype=int)
        Om2 = np.block([[Z, Id], [Id, Z]]) % 2

        M = (F.T % 2) @ Om2 % 2
        v = (inv_gfp(M, 2) @ Δ_half.reshape(-1, 1)) % 2
        v = v.flatten()

        x = v[:n]
        z = v[n:]
        pauli = PauliString.from_exponents(x, z, dimensions=dims)
        return PauliGate(pauli)


def gate_to_circuit(big_gate: "Gate", dimensions: list[int] | np.ndarray) -> "Circuit":
    """
    Decompose a single Gate (F, h) into a Circuit of Clifford generators
    (Hadamard/PHASE/CX/SWAP/…) followed by a final PauliGate that fixes the
    phase vector under the right-multiplication convention p' = p @ F^T.

    For odd p, the resulting circuit reproduces both F and h exactly.
    For p=2, the Pauli correction matches the even part of h; this is the
    known limitation of pure Pauli phase correction in qubit systems.

    Parameters
    ----------
    big_gate : Gate
        The gate to decompose. Can be a _GenericGate from composite_gate() or any other Gate.
    dimensions : list[int] | np.ndarray
        The dimensions of each qudit. Must be uniform (all equal).
    """
    # Pull target data
    F_target = big_gate.symplectic
    dims = list(map(int, dimensions))
    assert len(set(dims)) == 1, "gate_to_circuit assumes a uniform local dimension."
    p = dims[0]
    h_target = big_gate.phase_vector(p)

    n2 = F_target.shape[0]
    assert n2 % 2 == 0
    n = n2 // 2
    assert len(dims) == n, f"dimensions length {len(dims)} must match gate's n_qudits {n}"

    # 1) Build the Clifford "body" that realizes F_target (no Pauli yet)
    decomposition_no_corr = decompose_symplectic_to_circuit(F_target % p, p)

    # 2) Compute the *current* (F_body, h_body) after the body only
    body_gate = decomposition_no_corr.composite_gate()
    F_body = body_gate.symplectic % p
    h_body = body_gate.phase_vector(p)  # lives mod 2*lcm(dims)

    # Sanity: symplectic must already match
    if not np.array_equal(F_body % p, F_target % p):
        raise AssertionError("Body symplectic mismatch—decomposition error.")

    # 3) Append a Pauli correction so that h_body + F_body^T h_Pauli = h_target  (mod 2ℓ)
    P_corr = pauli_correction_gate(F_body, h_body, h_target, dims)

    # 4) Return full circuit: body then Pauli correction (right-multiplication order)
    # PauliGate acts on all qudits by default
    qudits = tuple(range(n))
    decomposition_no_corr.add_gate(P_corr, *qudits)
    return decomposition_no_corr
