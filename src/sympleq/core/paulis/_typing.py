from sympleq._typing import ScalarType, IntNDArray, ComplexNDArray, IntArrayVariant, ComplexArrayVariant

__all__ = [
    'ScalarType',
    'TableauType',
    'TableauVariant',
    'PhasesType',
    'PhasesVariant',
    'DimensionsType',
    'DimensionsVariant',
    'WeightsType',
    'WeightsVariant'
]

TableauType = IntNDArray
TableauVariant = IntArrayVariant

PhasesType = IntNDArray
PhasesVariant = IntArrayVariant

DimensionsType = IntNDArray
DimensionsVariant = IntArrayVariant

WeightsType = ComplexNDArray
WeightsVariant = ComplexArrayVariant
