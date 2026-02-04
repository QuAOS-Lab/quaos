from sympleq.core.paulis import PauliSum


class ToricCode:

    def __init__(self, Nx: int, Ny: int, c_x: float, c_z: float, c_g: float, periodic: bool = True):
        self.Nx = Nx
        self.Ny = Ny
        self.periodic = periodic
        self.c_x = c_x
        self.c_z = c_z
        self.c_g = c_g
        self.n_qubits = 2 * Nx * Ny if periodic else Nx * (Ny - 1) + Ny * (Nx - 1)

    def list_qubits(self) -> list[tuple[str, int, int]]:
        """
        List all qubits (edges) on an Nx-by-Ny lattice.
        Each qubit is denoted ('h', i, j) for a horizontal edge (i->i+1 at row j)
        or ('v', i, j) for a vertical edge (i, j->j+1).  Periodic=True wraps around.
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
        elif not self.periodic:
            assert len(qubits) == self.Nx * (self.Ny - 1) + self.Ny * (self.Nx - 1), "Non-periodic qubit count mismatch"
        else:
            raise ValueError("Invalid boundary condition")
        return qubits

    def build_star_ops(self) -> list[list[int]]:
        """
        Build the list of star operators.  Each star is a list of edge-indices
        (from list_qubits) on which Z acts.
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
        # Check
        assert len(star_ops) == self.Nx * self.Ny, "Star operator count mismatch"
        return star_ops

    def build_plaquette_ops(self) -> list[list[int]]:
        """
        Build the list of plaquette operators. Each plaquette is a list of
        edge-indices (from list_qubits) on which X acts.
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
                if True:
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
                    elif self.periodic and i == self.Nx - 1:
                        edges.append(index_of[('h', i, top_j)])
                # Left vertical edge
                qi = ('v', i, j)
                if qi in index_of:
                    edges.append(index_of[qi])
                elif self.periodic and j == self.Ny - 1:
                    edges.append(index_of[('v', i, j)])
                # Right vertical edge
                right_i = (i + 1) % self.Nx if self.periodic else i + 1
                if right_i < self.Nx:
                    qi = ('v', right_i, j)
                    if qi in index_of:
                        edges.append(index_of[qi])
                    elif self.periodic and j == self.Ny - 1:
                        edges.append(index_of[('v', right_i, j)])
                # Skip incomplete plaquettes on open boundary
                if not self.periodic and (i == self.Nx - 1 or j == self.Ny - 1):
                    continue
                if edges:
                    plaquette_ops.append(sorted(set(edges)))
        # Check
        if self.periodic:
            assert len(plaquette_ops) == self.Nx * self.Ny, "Periodic plaquette count mismatch"
        elif not self.periodic:
            assert len(plaquette_ops) == (self.Nx - 1) * (self.Ny - 1), "Non-periodic plaquette count mismatch"
        else:
            raise ValueError("Invalid boundary condition")
        return plaquette_ops

    def build_gauge_ops(self) -> list[list[int]]:
        """
        Build the list of gauge operators.  Each gauge operator consists of a single Z operator acting on a qubit.
        """
        qubits = self.list_qubits()
        gauge_ops = []
        for idx, q in enumerate(qubits):
            gauge_ops.append([idx])
        return gauge_ops

    def build_toric_code_hamiltonian(self) -> tuple[list[str], list[float]]:
        """
        Construct toric code Hamiltonian terms for given Nx, Ny, boundary, and coefficients.
        Returns (terms, coeffs), where each term is a string of 'xNzM' tokens for each qubit.
        """
        qubits = self.list_qubits()
        Nq = len(qubits)
        stars = self.build_star_ops()
        plaquettes = self.build_plaquette_ops()
        gauges = self.build_gauge_ops()
        terms = []
        coeffs = []
        # Star terms (Z on each edge in the star)
        if abs(self.c_z) > 10**-12:
            for edge_list in stars:
                word = []
                for q in range(Nq):
                    if q in edge_list:
                        word.append('x0z1')  # Z on this qubit
                    else:
                        word.append('x0z0')  # identity
                terms.append(' '.join(word))
                coeffs.append(self.c_z)
        # Plaquette terms (X on each edge in the plaquette)
        if abs(self.c_x) > 10**-12:
            for edge_list in plaquettes:
                word = []
                for q in range(Nq):
                    if q in edge_list:
                        word.append('x1z0')  # X on this qubit
                    else:
                        word.append('x0z0')  # identity
                terms.append(' '.join(word))
                coeffs.append(self.c_x)
        # Gauge terms (Z on each edge) - only if c_g != 0
        if abs(self.c_g) > 10**-12:
            for edge_list in gauges:
                word = []
                for q in range(Nq):
                    if q in edge_list:
                        word.append('x0z1')  # Z on this qubit
                    else:
                        word.append('x0z0')  # identity
                terms.append(' '.join(word))
                coeffs.append(self.c_g)
        return terms, coeffs

    def hamiltonian(self) -> PauliSum:
        ps, weights = self.build_toric_code_hamiltonian()
        return PauliSum.from_string(ps, weights=weights, dimensions=[2] * self.n_qubits)


if __name__ == "__main__":
    Nx = 2
    Ny = 2
    periodic = True
    c_x = 1.
    c_z = 1.
    c_g = 1.

    TC = ToricCode(Nx, Ny, c_x, c_z, c_g, periodic)
    hamiltonian = TC.hamiltonian()

    print(hamiltonian)
