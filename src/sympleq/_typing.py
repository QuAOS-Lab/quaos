from typing import Union
import numpy as np
from numpy.typing import NDArray

ScalarType = Union[float, complex, int]

IntNDArray = NDArray[np.integer]
IntArrayVariant = Union[int | list[int] | IntNDArray]

FloatNDArray = NDArray[np.floating]
FloatArrayVariant = Union[float | IntArrayVariant | list[float] | FloatNDArray]

ComplexNDArray = NDArray[np.complexfloating]
ComplexArrayVariant = Union[ScalarType | FloatArrayVariant | list[complex] | ComplexNDArray]
