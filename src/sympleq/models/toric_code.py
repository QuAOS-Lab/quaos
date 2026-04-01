import numpy as np
from sympleq.core.paulis import PauliSum


class ToricCode:

    def __init__(
        self,
        Nx: int,
        Ny: int,
        c_x: float,
        c_z: float,
        c_g: float | list[float] | np.ndarray,
        periodic: bool = True,
        d: int = 2,
        hermitian: bool = True,
        tol: float = 0.0,
    ):
        self.Nx = int(Nx)
        self.Ny = int(Ny)
        self.periodic = bool(periodic)
        self.c_x = float(c_x)
        self.c_z = float(c_z)

        self.d = int(d)
        if self.d < 2:
            raise ValueError("d must be >= 2.")
        self.hermitian = bool(hermitian)
        self.tol = float(tol)

        self.n_qubits = (
            2 * self.Nx * self.Ny
            if self.periodic
            else self.Nx * (self.Ny - 1) + self.Ny * (self.Nx - 1)
        )
        self.c_g = self._normalize_gauge_coeffs(c_g)

    def _normalize_gauge_coeffs(self, c_g: float | list[float] | np.ndarray) -> float | np.ndarray:
        if np.isscalar(c_g):
            return float(c_g)

        c_g_arr = np.asarray(c_g, dtype=float).reshape(-1)
        if c_g_arr.shape[0] != self.n_qubits:
            raise ValueError(
                f"Vector c_g must have length {self.n_qubits}, got {c_g_arr.shape[0]}."
            )
        return c_g_arr.copy()

    def _gauge_coeffs(self) -> np.ndarray:
        if np.isscalar(self.c_g):
            return np.full(self.n_qubits, float(self.c_g), dtype=float)
        return np.asarray(self.c_g, dtype=float).reshape(-1)

    def list_qubits(self) -> list[tuple[str, int, int]]:
        """
        List all edge qudits.
        ('h', i, j): horizontal edge (i,j)->(i+1,j)
        ('v', i, j): vertical edge   (i,j)->(i,j+1)
        """
        qubits = []
        # Horizontal edges
        for j in range(self.Ny):
            max_i = self.Nx if self.periodic else self.Nx - 1
            for i in range(max_i):
                qubits.append(('h', i, j))
        # Vertical edges
        for i in range(self.Nx):
            max_j = self.Ny if self.periodic else self.Ny - 1
            for j in range(max_j):
                qubits.append(('v', i, j))

        # Check
        if self.periodic:
            assert len(qubits) == 2 * self.Nx * self.Ny, "Periodic qubit count mismatch"
        else:
            assert len(qubits) == self.Nx * (self.Ny - 1) + self.Ny * (self.Nx - 1), "Non-periodic qubit count mismatch"
        return qubits

    def build_star_ops(self) -> list[list[int]]:
        """
        Qubit-style star supports (unsigned edge lists). For d>2, use build_star_ops_signed().
        """
        qubits = self.list_qubits()
        index_of = {q: idx for idx, q in enumerate(qubits)}
        star_ops = []
        for i in range(self.Nx):
            for j in range(self.Ny):
                edges = []
                # Horizontal edge to the right of (i,j)
                if (i < self.Nx - 1) or self.periodic:
                    qi = ('h', i, j)
                    if qi in index_of:
                        edges.append(index_of[qi])
                # Horizontal edge to the left of (i,j)
                i_left = (i - 1) % self.Nx if self.periodic else i - 1
                if (i > 0) or self.periodic:
                    qi = ('h', i_left, j)
                    if qi in index_of:
                        edges.append(index_of[qi])
                # Vertical edge above (i,j)
                if (j < self.Ny - 1) or self.periodic:
                    qi = ('v', i, j)
                    if qi in index_of:
                        edges.append(index_of[qi])
                # Vertical edge below (i,j)
                j_below = (j - 1) % self.Ny if self.periodic else j - 1
                if (j > 0) or self.periodic:
                    qi = ('v', i, j_below)
                    if qi in index_of:
                        edges.append(index_of[qi])
                if edges:
                    star_ops.append(sorted(set(edges)))
        assert len(star_ops) == self.Nx * self.Ny, "Star operator count mismatch"
        return star_ops

    def build_plaquette_ops(self) -> list[list[int]]:
        """
        Qubit-style plaquette supports (unsigned edge lists). For d>2, use build_plaquette_ops_signed().
        """
        qubits = self.list_qubits()
        index_of = {q: idx for idx, q in enumerate(qubits)}
        plaquette_ops = []
        max_i = self.Nx if self.periodic else self.Nx - 1
        max_j = self.Ny if self.periodic else self.Ny - 1
        for i in range(max_i):
            for j in range(max_j):
                edges = []
                # Bottom horizontal edge of plaquette at (i,j)
                qi = ('h', i, j)
                if qi in index_of:
                    edges.append(index_of[qi])
                elif self.periodic and i == self.Nx - 1:  # wrap around right
                    edges.append(index_of[('h', i, j)])
                # Top horizontal edge
                top_j = (j + 1) % self.Ny if self.periodic else j + 1
                if top_j < self.Ny:
                    qi = ('h', i, top_j)
                    if qi in index_of:
                        edges.append(index_of[qi])

                # Left vertical edge
                qi = ('v', i, j)
                if qi in index_of:
                    edges.append(index_of[qi])

                # Right vertical edge
                right_i = (i + 1) % self.Nx if self.periodic else i + 1
                if right_i < self.Nx:
                    qi = ('v', right_i, j)
                    if qi in index_of:
                        edges.append(index_of[qi])

                # Skip incomplete plaquettes on open boundary
                if not self.periodic and (i == self.Nx - 1 or j == self.Ny - 1):
                    continue
                if edges:
                    plaquette_ops.append(sorted(set(edges)))

        if self.periodic:
            assert len(plaquette_ops) == self.Nx * self.Ny, "Periodic plaquette count mismatch"
        else:
            assert len(plaquette_ops) == (self.Nx - 1) * (self.Ny - 1), "Non-periodic plaquette count mismatch"
        return plaquette_ops

    def build_gauge_ops(self) -> list[list[int]]:
        """
        Qubit-style single-site Z supports.
        """
        return [[idx] for idx, _ in enumerate(self.list_qubits())]

    def _index_of(self) -> dict[tuple[str, int, int], int]:
        qubits = self.list_qubits()
        return {q: idx for idx, q in enumerate(qubits)}

    @staticmethod
    def _combine_same_edge_exponents(edge_exp_list: list[tuple[int, int]], d: int) -> list[tuple[int, int]]:
        """
        Combine duplicate edges modulo d (important for tiny periodic lattices like Nx=1/Ny=1).
        """
        acc: dict[int, int] = {}
        for idx, e in edge_exp_list:
            acc[idx] = (acc.get(idx, 0) + int(e)) % d
        # drop zeros after mod reduction
        return [(idx, e) for idx, e in sorted(acc.items()) if (e % d) != 0]

    def build_star_ops_signed(self) -> list[list[tuple[int, int]]]:
        """
        Oriented star operators as Z-type exponent lists [(edge_index, exp_mod_d), ...].

        Convention:
          A_{i,j} = Z(h,i,j) * Z(v,i,j) * Z(h,i-1,j)^(-1) * Z(v,i,j-1)^(-1)
        """
        index_of = self._index_of()
        stars: list[list[tuple[int, int]]] = []
        d = self.d
        minus = (d - 1) % d

        for i in range(self.Nx):
            for j in range(self.Ny):
                terms: list[tuple[int, int]] = []

                # Outgoing horizontal edge to the right: +1
                if (i < self.Nx - 1) or self.periodic:
                    q = ('h', i, j)
                    if q in index_of:
                        terms.append((index_of[q], 1))

                # Incoming horizontal edge from left: -1
                i_left = (i - 1) % self.Nx if self.periodic else i - 1
                if (i > 0) or self.periodic:
                    q = ('h', i_left, j)
                    if q in index_of:
                        terms.append((index_of[q], minus))

                # Outgoing vertical edge upward: +1
                if (j < self.Ny - 1) or self.periodic:
                    q = ('v', i, j)
                    if q in index_of:
                        terms.append((index_of[q], 1))

                # Incoming vertical edge from below: -1
                j_below = (j - 1) % self.Ny if self.periodic else j - 1
                if (j > 0) or self.periodic:
                    q = ('v', i, j_below)
                    if q in index_of:
                        terms.append((index_of[q], minus))

                if terms:
                    stars.append(self._combine_same_edge_exponents(terms, d))

        assert len(stars) == self.Nx * self.Ny, "Star operator count mismatch"
        return stars

    def build_plaquette_ops_signed(self) -> list[list[tuple[int, int]]]:
        """
        Oriented plaquette operators as X-type exponent lists [(edge_index, exp_mod_d), ...].

        Counterclockwise convention:
          B_{i,j} = X(h,i,j) * X(v,i+1,j) * X(h,i,j+1)^(-1) * X(v,i,j)^(-1)
        """
        index_of = self._index_of()
        plaquettes: list[list[tuple[int, int]]] = []
        d = self.d
        minus = (d - 1) % d

        max_i = self.Nx if self.periodic else self.Nx - 1
        max_j = self.Ny if self.periodic else self.Ny - 1

        for i in range(max_i):
            for j in range(max_j):
                if not self.periodic and (i == self.Nx - 1 or j == self.Ny - 1):
                    continue

                terms: list[tuple[int, int]] = []

                # bottom edge h(i,j): +1
                q = ('h', i, j)
                if q not in index_of:
                    continue
                terms.append((index_of[q], 1))

                # right edge v(i+1,j): +1
                i_right = (i + 1) % self.Nx if self.periodic else i + 1
                if i_right >= self.Nx:
                    continue
                q = ('v', i_right, j)
                if q not in index_of:
                    continue
                terms.append((index_of[q], 1))

                # top edge h(i,j+1): -1
                j_top = (j + 1) % self.Ny if self.periodic else j + 1
                if j_top >= self.Ny:
                    continue
                q = ('h', i, j_top)
                if q not in index_of:
                    continue
                terms.append((index_of[q], minus))

                # left edge v(i,j): -1
                q = ('v', i, j)
                if q not in index_of:
                    continue
                terms.append((index_of[q], minus))

                plaquettes.append(self._combine_same_edge_exponents(terms, d))

        if self.periodic:
            assert len(plaquettes) == self.Nx * self.Ny, "Periodic plaquette count mismatch"
        else:
            assert len(plaquettes) == (self.Nx - 1) * (self.Ny - 1), "Non-periodic plaquette count mismatch"

        return plaquettes

    def build_gauge_ops_signed(self) -> list[list[tuple[int, int]]]:
        """Single-edge Z terms as signed exponent supports."""
        return [[(idx, 1)] for idx, _ in enumerate(self.list_qubits())]

    def _add_x_term(
        self,
        terms: dict[tuple[tuple[int, ...], tuple[int, ...]], complex],
        edge_exp: list[tuple[int, int]],
        coeff: complex
    ) -> None:
        n = self.n_qubits
        d = self.d
        x = np.zeros(n, dtype=int)
        z = np.zeros(n, dtype=int)
        for idx, e in edge_exp:
            x[idx] = (x[idx] + e) % d
        key = (tuple(x.tolist()), tuple(z.tolist()))
        terms[key] = terms.get(key, 0.0) + coeff

    def _add_z_term(
        self,
        terms: dict[tuple[tuple[int, ...], tuple[int, ...]], complex],
        edge_exp: list[tuple[int, int]],
        coeff: complex
    ) -> None:
        n = self.n_qubits
        d = self.d
        x = np.zeros(n, dtype=int)
        z = np.zeros(n, dtype=int)
        for idx, e in edge_exp:
            z[idx] = (z[idx] + e) % d
        key = (tuple(x.tolist()), tuple(z.tolist()))
        terms[key] = terms.get(key, 0.0) + coeff

    def _dagger_edge_exp(self, edge_exp: list[tuple[int, int]]) -> list[tuple[int, int]]:
        """Inverse exponents mod d (safe for pure-X / pure-Z strings)."""
        d = self.d
        return [(idx, (-e) % d) for idx, e in edge_exp]

    def _build_tableau_hamiltonian(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Build qudit toric code Hamiltonian as tableau + weights.
        Works for d=2 as well.
        """
        N = self.n_qubits
        d = self.d
        tol = self.tol

        stars = self.build_star_ops_signed()
        plaquettes = self.build_plaquette_ops_signed()
        gauges = self.build_gauge_ops_signed()

        # key: (x_tuple, z_tuple) -> coefficient
        terms: dict[tuple[tuple[int, ...], tuple[int, ...]], complex] = {}

        def add_term_and_hc(kind: str, edge_exp: list[tuple[int, int]], coeff: float):
            if abs(coeff) <= 1e-12:
                return

            if kind == "x":
                self._add_x_term(terms, edge_exp, complex(coeff))
                if self.hermitian:
                    edge_exp_dag = self._dagger_edge_exp(edge_exp)
                    if edge_exp_dag != edge_exp:  # avoids double count at d=2
                        self._add_x_term(terms, edge_exp_dag, complex(coeff))
            elif kind == "z":
                self._add_z_term(terms, edge_exp, complex(coeff))
                if self.hermitian:
                    edge_exp_dag = self._dagger_edge_exp(edge_exp)
                    if edge_exp_dag != edge_exp:
                        self._add_z_term(terms, edge_exp_dag, complex(coeff))
            else:
                raise ValueError("kind must be 'x' or 'z'.")

        # Star terms (Z-type)
        if abs(self.c_z) > 1e-12:
            for edge_exp in stars:
                add_term_and_hc("z", edge_exp, self.c_z)

        # Plaquette terms (X-type)
        if abs(self.c_x) > 1e-12:
            for edge_exp in plaquettes:
                add_term_and_hc("x", edge_exp, self.c_x)

        # Gauge terms (single-edge Z)
        gauge_coeffs = self._gauge_coeffs()
        if np.any(np.abs(gauge_coeffs) > 1e-12):
            for edge_exp, coeff in zip(gauges, gauge_coeffs):
                add_term_and_hc("z", edge_exp, float(coeff))

        items = []
        for (x_t, z_t), c in terms.items():
            if tol > 0.0 and abs(c) <= tol:
                continue
            if c != 0:
                items.append(((x_t, z_t), c))

        items.sort(key=lambda kv: (kv[0][0], kv[0][1]))

        M = len(items)
        tableau = np.zeros((M, 2 * N), dtype=int)
        weights = np.zeros(M, dtype=complex)

        for row, ((x_t, z_t), c) in enumerate(items):
            tableau[row, :N] = np.array(x_t, dtype=int) % d
            tableau[row, N:] = np.array(z_t, dtype=int) % d
            weights[row] = c

        return tableau, weights

    def build_toric_code_hamiltonian(self) -> tuple[list[str], list[float]]:
        """
        Legacy qubit-only string builder.
        For d>2, use hamiltonian() (tableau-based qudit path).
        """
        if self.d != 2:
            raise NotImplementedError(
                "build_toric_code_hamiltonian() is the legacy qubit string builder. "
                "For d>2 use hamiltonian() which builds a tableau-based qudit PauliSum."
            )

        qubits = self.list_qubits()
        Nq = len(qubits)
        stars = self.build_star_ops()
        plaquettes = self.build_plaquette_ops()
        gauges = self.build_gauge_ops()
        terms = []
        coeffs = []

        # Star terms (Z on each edge in the star)
        if abs(self.c_z) > 1e-12:
            for edge_list in stars:
                word = []
                edge_set = set(edge_list)
                for q in range(Nq):
                    word.append('x0z1' if q in edge_set else 'x0z0')
                terms.append(' '.join(word))
                coeffs.append(self.c_z)

        # Plaquette terms (X on each edge in the plaquette)
        if abs(self.c_x) > 1e-12:
            for edge_list in plaquettes:
                word = []
                edge_set = set(edge_list)
                for q in range(Nq):
                    word.append('x1z0' if q in edge_set else 'x0z0')
                terms.append(' '.join(word))
                coeffs.append(self.c_x)

        # Gauge terms (Z on each edge)
        gauge_coeffs = self._gauge_coeffs()
        if np.any(np.abs(gauge_coeffs) > 1e-12):
            for edge_list, coeff in zip(gauges, gauge_coeffs):
                word = []
                edge_set = set(edge_list)
                for q in range(Nq):
                    word.append('x0z1' if q in edge_set else 'x0z0')
                terms.append(' '.join(word))
                coeffs.append(float(coeff))

        return terms, coeffs

    def hamiltonian(self) -> PauliSum:
        """
        Returns a PauliSum.
        - d=2: uses your legacy qubit string path (backward compatible)
        - d>2: uses tableau-based qudit path with oriented star/plaquette operators
        """
        if self.d == 2:
            ps, weights = self.build_toric_code_hamiltonian()
            return PauliSum.from_string(ps, weights=weights, dimensions=[2] * self.n_qubits)

        tableau, weights = self._build_tableau_hamiltonian()
        return PauliSum.from_tableau(
            tableau,
            dimensions=[self.d] * self.n_qubits,
            weights=weights,
        )
