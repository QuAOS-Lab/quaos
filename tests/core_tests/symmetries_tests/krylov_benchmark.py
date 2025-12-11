import numpy as np
import pytest

from sympleq.models.Ising import ising_chain_hamiltonian, heuristic_clifford_symmetry
from sympleq.core.paulis import PauliSum
from sympleq.core.circuits import Gate
from sympleq.core.symmetries.block_decomposition import block_decompose_optimal
from sympleq.core.symmetries.clifford import clifford_phase_decomposition
from sympleq.core.symmetries.krylov_time_evolve import (
    observable_dynamics_krylov,
)
from sympleq.core.symmetries.krylov_time_evolve import (
    krylov_observable_dynamics_plain,
)
from sympleq.core.states.state import State
from sympleq.core.circuits.gates import Hadamard as H_gate


def _ising_X0_observable(N: int) -> PauliSum:
    tableau = np.zeros((1, 2 * N), dtype=int)
    tableau[0, 0] = 1  # X on qudit 0
    return PauliSum.from_tableau(tableau, weights=[1.0], dimensions=[2] * N)


def _ising_product_state_x_plus(N: int) -> np.ndarray:
    # Reuse your existing helper if you prefer; this is a self-contained version.


    dims = [2] * N
    # Start from |0...0>
    psi0 = np.zeros(2**N, dtype=np.complex128)
    psi0[0] = 1.0
    state = State(psi0, dims)

    # Apply local Hadamard on each qubit to get |+...+>
    for q in range(N):
        state = H_gate(q, 2).act_on_state(state)

    return state.as_array()


@pytest.mark.benchmark(group="krylov_plain_vs_symmetry")
@pytest.mark.parametrize("N", [8, 10])
def test_benchmark_plain_krylov(benchmark, N: int):
    """
    Benchmark plain Krylov ⟨X_0(t)⟩ for an N-site Ising chain.
    """
    J = 1.0
    h = 0.5
    H = ising_chain_hamiltonian(N, J, h, periodic=True)
    O = _ising_X0_observable(N)
    psi0 = _ising_product_state_x_plus(N)
    times = np.linspace(0.0, 2.0, 101)

    def run():
        return krylov_observable_dynamics_plain(
            H=H,
            O=O,
            psi0=psi0,
            times=times,
            m_max=64,
        )

    result = benchmark(run)
    assert result.shape == times.shape


@pytest.mark.benchmark(group="krylov_plain_vs_symmetry")
@pytest.mark.parametrize("N", [8, 10])
def test_benchmark_symmetry_krylov_full(benchmark, N: int):
    """
    Benchmark symmetry-aware Krylov (full Hilbert) ⟨X_0(t)⟩ for an N-site Ising chain.
    Uses heuristic_clifford_symmetry + block_decompose_optimal + T/S reconstruction.
    """
    J = 1.0
    h = 0.5
    H = ising_chain_hamiltonian(N, J, h, periodic=True)
    O = _ising_X0_observable(N)
    psi0 = _ising_product_state_x_plus(N)
    times = np.linspace(0.0, 2.0, 101)

    # Symmetry construction (same as in your tests/examples)
    F = heuristic_clifford_symmetry(N)
    S_symp, T_symp = block_decompose_optimal(F.symplectic, 2)
    h_S, h_T = clifford_phase_decomposition(
        F.symplectic, F.phase_vector, S_symp, T_symp, int(H.lcm)
    )
    S_gate = Gate("S", F.qudit_indices, S_symp, F.dimensions, h_S)
    T_gate = Gate("T", F.qudit_indices, T_symp, F.dimensions, h_T)

    def run():
        return observable_dynamics_krylov(
            H=H,
            O=O,
            psi0=psi0,
            times=times,
            T_gate=T_gate,
            S_gate=S_gate,
            symmetry_eigval=None,
            use_reduced_basis=False,
            m_max=64,
        )

    result = benchmark(run)
    assert result.shape == times.shape


@pytest.mark.benchmark(group="krylov_plain_vs_symmetry")
@pytest.mark.parametrize("N", [8, 10])
def test_benchmark_symmetry_krylov_reduced(benchmark, N: int):
    """
    Benchmark symmetry-reduced Krylov ⟨X_0(t)⟩ for an N-site Ising chain.
    This uses the sector-reduced algorithm (still with dense sector basis).
    """
    J = 1.0
    h = 0.5
    H = ising_chain_hamiltonian(N, J, h, periodic=True)
    O = _ising_X0_observable(N)
    psi0 = _ising_product_state_x_plus(N)
    times = np.linspace(0.0, 2.0, 101)

    F = heuristic_clifford_symmetry(N)
    S_symp, T_symp = block_decompose_optimal(F.symplectic, 2)
    h_S, h_T = clifford_phase_decomposition(
        F.symplectic, F.phase_vector, S_symp, T_symp, int(H.lcm)
    )
    S_gate = Gate("S", F.qudit_indices, S_symp, F.dimensions, h_S)
    T_gate = Gate("T", F.qudit_indices, T_symp, F.dimensions, h_T)

    def run():
        return observable_dynamics_krylov(
            H=H,
            O=O,
            psi0=psi0,
            times=times,
            T_gate=T_gate,
            S_gate=S_gate,
            symmetry_eigval=1.0,  # choose +1 sector for benchmark
            use_reduced_basis=True,
            m_max=64,
        )

    result = benchmark(run)
    assert result.shape == times.shape
