from typing import Union
import numpy as np
from numpy.typing import NDArray
import scipy.sparse as sp

ScalarType = Union[float, complex, int]
ArrayShape = int
MatrixShape = tuple[int, int]

IntNDArray = NDArray[np.integer]
IntArrayLike = Union[int | list[int] | IntNDArray]

FloatNDArray = NDArray[np.floating]
FloatArrayLike = Union[float | IntArrayLike | list[float] | FloatNDArray]

ComplexNDArray = NDArray[np.complexfloating]
ComplexArrayLike = Union[ScalarType | FloatArrayLike | list[complex] | ComplexNDArray]

IntSparseMatrix = sp.csr_matrix
IntSparseMatrixLike = Union[
    IntSparseMatrix | sp.coo_matrix | tuple[IntArrayLike, IntArrayLike, IntArrayLike, MatrixShape]
]

FloatSparseMatrix = sp.csr_matrix
FloatSparseMatrixLike = Union[
    FloatSparseMatrix | IntSparseMatrixLike | tuple[FloatArrayLike, IntArrayLike, IntArrayLike, MatrixShape]
]

ComplexSparseMatrix = sp.csr_matrix
ComplexSparseMatrixLike = Union[
    ComplexSparseMatrix | FloatSparseMatrixLike | tuple[ComplexArrayLike, IntArrayLike, IntArrayLike, MatrixShape]
]
