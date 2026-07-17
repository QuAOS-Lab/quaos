# flake8: noqa
"""
Build the transvection sequence needed to map one PauliSum in GF(p) to another.
This currently assumes prime qudit dimension.

Future extensions: input symplectic matrix -> output symplectic matrix.
"""

import numpy as np
import galois
from itertools import islice
from sympleq.core.circuits.utils import transvection_matrix, symplectic_product_arrays
from sympleq.core.circuits.find_symplectic import check_mappable_via_clifford


def modinv(a, p):
    """Multiplicative inverse of a mod p (works for prime p)."""
    return pow(int(a), -1, p)


def solve_gfp(A: np.ndarray, b: np.ndarray, p: int):
    """
    TODO: Make this file use the solver in finite_field_solving.py

    Solve Ax = b over GF(p) using Gaussian elimination.
    Returns one valid solution vector x (length n), or None if no solution exists.

    Args:
        A : coefficient matrix (m x n)
        b : right-hand side vector (m,)
        p : prime

    Returns:
        x : one solution in GF(p)^n as np.ndarray, or None
    """
    A = A.astype(int) % p
    b = b.astype(int) % p
    m, n = A.shape

    # Augmented matrix
    Ab = np.hstack([A, b.reshape(-1, 1)]) % p

    pivot_row = 0
    pivot_cols = []

    for col in range(n):
        # Find pivot
        pivot = None
        for row in range(pivot_row, m):
            if Ab[row, col] != 0:
                pivot = row
                break
        if pivot is None:
            continue

        # Swap into place
        if pivot != pivot_row:
            Ab[[pivot_row, pivot]] = Ab[[pivot, pivot_row]]

        # Normalize pivot row
        inv = modinv(Ab[pivot_row, col], p)
        Ab[pivot_row] = (Ab[pivot_row] * inv) % p

        # Eliminate other rows
        for row in range(m):
            if row != pivot_row and Ab[row, col] != 0:
                factor = Ab[row, col]
                Ab[row] = (Ab[row] - factor * Ab[pivot_row]) % p

        pivot_cols.append(col)
        pivot_row += 1
        if pivot_row == m:
            break

    # Check inconsistency
    for row in range(pivot_row, m):
        if np.all(Ab[row, :-1] == 0) and Ab[row, -1] != 0:
            return None  # no solution

    # Build one solution: set free vars = 0
    x = np.zeros(n, dtype=int)
    for i, col in enumerate(pivot_cols):
        row = i
        val = Ab[row, -1]
        for j in range(col + 1, n):
            val = (val - Ab[row, j] * x[j]) % p
        x[col] = val % p

    return x


def pair(vec, i, n):
    """Return 2x1 pair at index i"""
    return np.array([vec[i], vec[i+n]])


def solve_non_proportional(u_i, v_i, p):
    """Solve B x = y for non-proportional 2x2 vectors"""

    GF = galois.GF(p)
    B = np.array([[u_i[0], -u_i[1]], [-v_i[0], v_i[1]]]) % p
    y = GF([1, 1])
    B = GF(B)
    x = np.linalg.solve(B, y).astype(np.int32).view(np.ndarray) % p
    return x


def proportionality_constant_mod(u_i, v_i, p):
    """
    Compute the proportionality constant k mod p such that v_i = k * u_i (mod p),
    assuming u_i and v_i are proportional.
    """
    a, b = u_i
    c, d = v_i

    if not galois.is_prime(p):
        raise ValueError(f"{p} is not prime!")

    # Convert to Python ints to avoid NumPy overflow
    a, b, c, d = int(a), int(b), int(c), int(d)

    # Case 1: all entries nonzero
    if a != 0 and b != 0 and c != 0 and d != 0:
        k1 = (a * pow(c, p-2, p)) % p  # k = a / c mod p
        k2 = (b * pow(d, p-2, p)) % p  # k = b / d mod p
        if k1 != k2:
            raise ValueError("Vectors are not proportional")
        return k1

    # Case 2: only second component nonzero
    elif b != 0 and d != 0:
        k = (b * pow(d, p-2, p)) % p
        return k

    # Case 3: only first component nonzero
    elif a != 0 and c != 0:
        k = (a * pow(c, p-2, p)) % p
        return k

    else:
        raise ValueError("Cannot determine proportionality constant: one of the components is zero")


def solve_proportional(u_i, v_i, p, proportionality_const):
    """Solve B x = y for proportional pair with free parameter"""

    GF = galois.GF(p)
    B = np.array([[u_i[0], -u_i[1]], [-v_i[0], v_i[1]]]) % p
    y = GF([1, proportionality_const])
    B = GF(B)
    row1 = B[0]
    t_free = proportionality_const
    solutions = []

    if row1[0] != 0 and row1[1] != 0:
        x1 = (y[0] - row1[1]*t_free) * GF(int(row1[0]))**(-1)
        solutions.append(GF([x1, t_free]))
    elif row1[0] != 0 and row1[1] == 0:
        x2 = (y[0] - row1[0]*t_free) * GF(int(row1[0]))**(-1)
        solutions.append(GF([t_free, x2]))
    elif row1[0] == 0 and row1[1] != 0:
        x2 = (y[0] - row1[0]*t_free) * GF(int(row1[1]))**(-1)
        solutions.append(GF([x2, t_free]))

    return solutions[0].astype(np.int32).view(np.ndarray)


def build_symplectic_for_transvection(u, v, p):
    """
    Build the intermediate vector w for two-step transvections if [u, v]=0 mod p 
    Assumes qudit dimension p is prime and symplectic_product_arrays(u, v, p) == 0
    """
    # --- Input validation ---
    if np.all(u == 0) or np.all(v == 0):
        raise ValueError("Input vectors cannot be zero")
    if not galois.is_prime(p):
        raise NotImplementedError(f"Prime dimension expected, got p={p}")
    if symplectic_product_arrays(u, v, p) != 0:
        raise ValueError("Symplectic product must be zero")

    n = len(u) // 2
    w = np.zeros(2*n, dtype=int)
    # GF = galois.GF(p)

    # --- First pass: find non-proportional pair ---
    for i, (ui, vi) in enumerate(islice(zip(u, v), n)):
        u_i = pair(u, i, n)
        v_i = pair(v, i, n)
        nz_u, nz_v = np.count_nonzero(u_i), np.count_nonzero(v_i)
        sp_uv = symplectic_product_arrays(u_i, v_i, p)

        if nz_u and nz_v and sp_uv != 0:
            x = solve_non_proportional(u_i, v_i, p)
            w[i+n], w[i] = x[0], x[1]
            break

    # --- Second pass: proportional pair ---
    if np.all(w == 0):
        for i, (ui, vi) in enumerate(islice(zip(u, v), n)):
            u_i = pair(u, i, n)
            v_i = pair(v, i, n)
            nz_u, nz_v = np.count_nonzero(u_i), np.count_nonzero(v_i)
            sp_uv = symplectic_product_arrays(u_i, v_i, p)

            if nz_u and nz_v and sp_uv == 0:
                k = proportionality_constant_mod(u_i, v_i, p)
                x = solve_proportional(u_i, v_i, p,  k)
                w[i+n], w[i] = x[0], x[1]
                break

    # --- Third pass: handle single nonzero entries manually ---
    if np.all(w == 0):
        for i in range(n):
            u_i, v_i = pair(u, i, n), pair(v, i, n)

            # Case u_i nonzero, v_i zero
            if np.count_nonzero(u_i) != 0 and np.count_nonzero(v_i) == 0:
                w_u = np.zeros(2, dtype=int)
                if u_i[0] != 0:
                    w_u[1] = modinv(u_i[0],p)
                else:
                    w_u[0] = modinv(-u_i[1],p)
                assert symplectic_product_arrays(u_i, w_u, p) == 1
                w[i], w[i+n] = w_u[0], w_u[1]
                break

        for i in range(n):
            u_i, v_i = pair(u, i, n), pair(v, i, n)

            # Case u_i zero, v_i nonzero
            if np.count_nonzero(u_i) == 0 and np.count_nonzero(v_i) != 0:
                w_v = np.zeros(2, dtype=int)
                if v_i[0] != 0:
                    w_v[1] = (-modinv(v_i[0],p)) % p
                else:
                    w_v[0] = modinv(v_i[1],p)
                assert symplectic_product_arrays(w_v, v_i, p) == 1
                w[i], w[i+n] = w_v[0], w_v[1]
                break

    # --- Final validation ---
    if symplectic_product_arrays(u, w, p) == 0 or symplectic_product_arrays(w, v, p) == 0:
        raise ValueError("Failed to construct valid symplectic vector w")

    return w


def find_transvection_map(input_ps, output_ps, p):
    """
    Provides the map to transfer one paulisum in GF(p) to another paulisum.
    Currently works only for prime dimension. 
    Returns combined transvection (final map)
    """

    # --- Input validation ---
    if np.all(input_ps == 0) or np.all(output_ps == 0):
        raise ValueError("Vectors cannot be zero")
    if not galois.is_prime(p):
        raise NotImplementedError(f"Prime dimension expected, got p={p}")

    n = int(len(input_ps)//2)
    a = symplectic_product_arrays(input_ps, output_ps,p)

    if a != 0:
        ainv = modinv(a, p)
        h=(-input_ps + output_ps) % p
        F_h= transvection_matrix(h, p, multiplier=ainv)
        if not np.array_equal((input_ps @ F_h) % p, output_ps):
            raise ValueError("Failed to construct valid transvection for nonzero symplectic product")

    elif a == 0:
        w = build_symplectic_for_transvection(input_ps, output_ps, p)
        a_w = symplectic_product_arrays(input_ps, w,p)
        a_w_inv = modinv(a_w, p)
        h = (-input_ps + w) % p
        F_h_1 = transvection_matrix(h, p, multiplier=a_w_inv)
        if not np.array_equal((input_ps @ F_h_1) % p, w):
            raise ValueError("Failed to construct valid transvection for u->w")
        b_w = symplectic_product_arrays(w, output_ps, p)
        b_w_inv = modinv(b_w, p)
        h = (-w + output_ps) % p
        F_h_2 = transvection_matrix(h, p, multiplier=b_w_inv)
        if not np.array_equal((w @ F_h_2) % p, output_ps):
            raise ValueError("Failed to construct valid transvection for w->v")
        F_h = F_h_1 @ F_h_2 %p
        if not np.array_equal((input_ps @ F_h) % p, output_ps):
            raise ValueError("Failed to construct valid transvection map: for u->v")

    return F_h


def intermediate_transvection_solve(u, v, p):
    # solving in GF(p) instead of building for symplectic product = 0 case
    # --- Input validation ---
    if np.all(u == 0) or np.all(v == 0):
        raise ValueError("Input vectors cannot be zero")
    if not galois.is_prime(p):
        raise NotImplementedError(f"Prime dimension expected, got p={p}")
    if symplectic_product_arrays(u, v, p) != 0:
        raise ValueError("Symplectic product must be zero")

    n= int(len(u)//2)

    A = np.zeros((2, 2 * n), dtype=int)           # The structure of A is such that [ux|-uz
                                                                                    #  -vx|vz],
                                                    #i.e, A @ w1 becomes [uxw1z- w1xuz
                                                    #                      w1xvz- w1zvx]

    w= np.zeros(2 * n, dtype=int)

    A[0, :n] = u[:n]  # X part of u multiplies Z part of w
    A[0, n:] = -u[n:]  # Z part of u multiplies X part of w

    A[1, :n] = -v[:n]  # X part of v multiplies Z part of w
    A[1, n:] = v[n:]  # Z part of v multiplies X part of w

    for p1, p2 in zip(range(1, p), range(1, p)):
        b = np.array([p1, p2])                 #The structure of b is : [symplectic product(u, w, p)
                                                    #                         symplectic product(w, v, p)]
        w1=solve_gfp(A, b, p)

        if w1 is not None:                       # The solution will be [w1z
            #                        w1x], so, we will change it to w]
            w[:n]=w1[n:]
            w[n:]= w1[:n]
            break

    #assert symplectic_product_arrays(u, w, p)!=0 and symplectic_product_arrays(w, v,p) !=0, f'{symplectic_product_arrays(u, w, p), symplectic_product_arrays(w, v, p)}'

    return w


def find_transvection_map_solve(input_ps, output_ps, p):
    """
    Provides the map to transfer one PauliString in GF(p) to another PauliString.
    Currently works only for prime dimension. 
    Returns combined transvection (final map)
    """

    # --- Input validation ---
    if np.all(input_ps == 0) or np.all(output_ps == 0):
        raise ValueError("Vectors cannot be zero")
    if not galois.is_prime(p):
        raise NotImplementedError(f"Prime dimension expected, got p={p}")

    a=symplectic_product_arrays(input_ps, output_ps,p)
    if a != 0:
        ainv = modinv(a, p)
        h=(-input_ps + output_ps) % p
        F_h= transvection_matrix(h, p, multiplier=ainv)
        if not np.array_equal((input_ps @ F_h) % p, output_ps):
            raise ValueError("Failed to construct valid transvection for nonzero symplectic product")

    elif a == 0:
        w=intermediate_transvection_solve(input_ps, output_ps, p)
        a_w= symplectic_product_arrays(input_ps, w,p)
        a_w_inv = modinv(a_w, p)
        h=(-input_ps + w) % p
        F_h_1= transvection_matrix(h, p, multiplier=a_w_inv)
        if not np.array_equal((input_ps @ F_h_1) % p, w):
            raise ValueError("Failed to construct valid transvection for u->w")
        b_w= symplectic_product_arrays(w, output_ps, p)
        b_w_inv = modinv(b_w, p)
        h=(-w + output_ps) % p
        F_h_2= transvection_matrix(h, p, multiplier=b_w_inv)
        if not np.array_equal((w @ F_h_2) % p, output_ps):
            raise ValueError("Failed to construct valid transvection for w->v")
        F_h= F_h_1 @ F_h_2 %p
        if not np.array_equal((input_ps @ F_h) % p, output_ps):
            raise ValueError("Failed to construct valid transvection map: for u->v")

    return F_h


def intermediate_transvection_solve_extended(u, v, constraints, sps,  p):

    # solving in GF(p) instead of building for symplectic product = 0 case
    # --- Input validation ---
    if np.all(u == 0) or np.all(v == 0):
        raise ValueError("Input vectors cannot be zero")
    if np.all(u==v):
        raise ValueError("Input PauliSum == Output PauliSum")

    if not galois.is_prime(p):
        raise NotImplementedError(f"Prime dimension expected, got p={p}")

    if symplectic_product_arrays(u, v, p) != 0:
        raise ValueError("Symplectic product must be zero")

    else:
        t= len(constraints)
        n= int(len(u)//2)
        # Set up the linear system A @ w1 = b

        A = np.zeros((2+t, 2 * n), dtype=int)           # The structure of A is such that [ux|-uz
                                                                                        #  -vx|vz],
                                                        #i.e, A @ w1 becomes [uxw1z- w1xuz
                                                        #                      w1xvz- w1zvx]

        w= np.zeros(2 * n, dtype=int)

        A[0, :n] = u[:n]  # X part of u multiplies Z part of w
        A[0, n:] = -u[n:]  # Z part of u multiplies X part of w

        A[1, :n] = -v[:n]  # X part of v multiplies Z part of w
        A[1, n:] = v[n:]  # Z part of v multiplies X part of w

        for el in range(t):
            A[2+el, :n] = constraints[el][:n]
            A[2+el, n:]= -constraints[el][n:]

        b= np.zeros(2+t, dtype=int)
        for el in range(t):
            b[2+el]= sps[el]

        for p1 in range(1, p):
            for p2 in range(1, p):
                b[0] = p1              #The structure of b is : [symplectic product(u, w, p)
                b[1] = p2               #                         symplectic product(w, v, p)]
                w1=solve_gfp(A, b, p)
                if w1 is not None:
                    break
            break
        if w1 is None:
            raise Exception("Could not find transvection vector.")

        w[:n]=w1[n:]              # The solution will be [w1z
                                  #                        w1x], so, we will change it to w]
        w[n:]= w1[:n]

        if np.all(w==0):
            raise ValueError('No solution found')
        else:
            assert symplectic_product_arrays(u, w, p)!=0 and symplectic_product_arrays(
                w, v,p) !=0, f'{symplectic_product_arrays(u, w, p), symplectic_product_arrays(w, v, p)}'

    return w

def find_transvection_map_solve_extended(input_ps, output_ps, constraints=[], sps=[], p=2):
    """
    Provides the map to transfer one PauliString in GF(p) to another PauliString, subject to some constraints.
    Currently works only for prime dimension. 
    Returns combined transvection (final map)
    """

    # --- Input validation ---
    if np.all(input_ps == 0) or np.all(output_ps == 0):
        raise ValueError("Vectors cannot be zero")
    if not galois.is_prime(p):
        raise NotImplementedError(f"Prime dimension expected, got p={p}")

    a= symplectic_product_arrays(input_ps, output_ps,p)
    if a != 0:
        ainv = modinv(a, p)
        h=(-input_ps + output_ps) % p
        F_h= transvection_matrix(h, p, multiplier=ainv)
        if not np.array_equal((input_ps @ F_h) % p, output_ps):
            raise ValueError("Failed to construct valid transvection for nonzero symplectic product")

    elif a == 0:
        if np.all(input_ps==output_ps):
            n=int(len(input_ps)//2)
            F_h= np.eye(2 * n, dtype=int)
        else:
            w=intermediate_transvection_solve_extended(input_ps, output_ps, constraints=constraints, sps=sps, p=p)
            a_w= symplectic_product_arrays(input_ps, w,p)
            a_w_inv = modinv(a_w, p)
            h=(-input_ps + w) % p
            F_h_1= transvection_matrix(h, p, multiplier=a_w_inv)
            if not np.array_equal((input_ps @ F_h_1) % p, w):
                raise ValueError("Failed to construct valid transvection for u->w")
            b_w= symplectic_product_arrays(w, output_ps, p)
            b_w_inv = modinv(b_w, p)
            h=(-w + output_ps) % p
            F_h_2= transvection_matrix(h, p, multiplier=b_w_inv)
            if not np.array_equal((w @ F_h_2) % p, output_ps):
                raise ValueError("Failed to construct valid transvection for w->v")
            F_h= F_h_1 @ F_h_2 %p
            if not np.array_equal((input_ps @ F_h) % p, output_ps):
                raise ValueError("Failed to construct valid transvection map: for u->v")

    return F_h


def map_paulisum_to_target_paulisum(input_tab, output_tab, p):

    if not check_mappable_via_clifford(input_tab, output_tab, p):
        raise Exception(f'Cannot map these tabs')

    # Map the first Pauli
    F_total = find_transvection_map_solve_extended(
        input_tab[0], output_tab[0], constraints=[], sps=[], p=p
    )
    assert (input_tab[0] @ F_total % p == output_tab[0]).all()

    # Iteratively map the rest of the Paulis
    for k in range(1, len(input_tab)):
        u_k = input_tab[k] @ F_total % p

        # Constraints: previously mapped Paulis must stay fixed
        constraints = [output_tab[j] for j in range(k)]
        sps = [symplectic_product_arrays(input_tab[j], input_tab[k], p) for j in range(k)]

        F_k = find_transvection_map_solve_extended(
            u_k, output_tab[k], constraints=constraints, sps=sps, p=p
        )

        # Update cumulative Clifford map
        F_total = (F_total @ F_k) % p

        # Verify this step
        assert (u_k @ F_k % p == output_tab[k]).all()
        for j in range(k):
            assert (output_tab[j] @ F_k % p == output_tab[j]).all()

    return F_total
