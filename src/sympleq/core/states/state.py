import numpy as np
from typing import Sequence


class State:
    """
    Pure state in the computational basis for a mixed-radix register.

    Attributes
    ----------
    amplitudes : np.ndarray
        Complex vector of length prod(dimensions).
    dimensions : np.ndarray
        1D int array of local dimensions [d0, d1, ..., d_{n-1}].
    """
    def __init__(self, amplitudes: np.ndarray, dimensions: Sequence[int] | np.ndarray):
        dims = np.asarray(dimensions, dtype=int).reshape(-1)
        D = int(np.prod(dims))

        amps = np.asarray(amplitudes, dtype=np.complex128).reshape(-1)
        if amps.shape != (D,):
            raise ValueError(
                f"State amplitudes have shape {amps.shape}, but dimensions {dims} "
                f"imply Hilbert space dimension {D}."
            )

        self.amplitudes = amps
        self.dimensions = dims

    @classmethod
    def from_basis(cls, dimensions: Sequence[int] | np.ndarray, index: int) -> "State":
        """
        |index> in the lexicographically-ordered computational basis.
        """
        dims = np.asarray(dimensions, dtype=int).reshape(-1)
        D = int(np.prod(dims))
        if not (0 <= index < D):
            raise ValueError(f"Basis index {index} out of range for dimension {D}.")
        amps = np.zeros(D, dtype=np.complex128)
        amps[index] = 1.0
        return cls(amps, dims)

    def copy(self) -> "State":
        return State(self.amplitudes.copy(), self.dimensions.copy())

    def n_qudits(self) -> int:
        return self.dimensions.size

    def as_array(self) -> np.ndarray:
        return self.amplitudes
