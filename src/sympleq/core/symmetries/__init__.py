from importlib import import_module

__all__ = ['min_qudit_clifford_symmetry', 'find_clifford_symmetries',
           'qudit_cost', 'pauli_reduce']


def __getattr__(name):
    if name in {'min_qudit_clifford_symmetry', 'find_clifford_symmetries', 'qudit_cost'}:
        clifford = import_module(".clifford", __name__)
        return getattr(clifford, name)
    if name == 'pauli_reduce':
        pauli_reduce = import_module(".pauli", __name__).pauli_reduce
        return pauli_reduce
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
