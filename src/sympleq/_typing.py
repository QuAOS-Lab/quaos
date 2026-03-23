from typing import Union
import numpy as np
from numpy.typing import NDArray
import scipy.sparse as sp

ScalarType = Union[float, complex, int]

IntNDArray = NDArray[np.integer]
IntArrayLike = Union[int | list[int] | IntNDArray]

FloatNDArray = NDArray[np.floating]
FloatArrayLike = Union[float | IntArrayLike | list[float] | FloatNDArray]

ComplexNDArray = NDArray[np.complexfloating]
ComplexArrayLike = Union[ScalarType | FloatArrayLike | list[complex] | ComplexNDArray]

ComplexSparseMatrix = sp.csr_matrix
