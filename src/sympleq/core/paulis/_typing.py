from sympleq._typing import (
    ScalarType, IntNDArray, ComplexNDArray, IntArrayLike, ComplexArrayLike, ComplexSparseMatrix
)

__all__ = [
    'ScalarType',
    'TableauType',
    'TableauLike',
    'PhasesType',
    'PhasesLike',
    'DimensionsType',
    'DimensionsLike',
    'WeightsType',
    'WeightsLike'
]

TableauType = IntNDArray
TableauLike = IntArrayLike

PhasesType = IntNDArray
PhasesLike = IntArrayLike

DimensionsType = IntNDArray
DimensionsLike = IntArrayLike

WeightsType = ComplexNDArray
WeightsLike = ComplexArrayLike

HilbertOperator = ComplexSparseMatrix
