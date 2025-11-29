import os
import sys


def _setup_native_libs():
    """Ensure native libs are found on all platforms."""
    here = os.path.dirname(__file__)
    lib_dir = os.path.join(os.path.dirname(here), "sympleq.libs")

    if os.name == "nt":
        if os.path.isdir(lib_dir):
            os.add_dll_directory(lib_dir)


_setup_native_libs()

import sympleq.core.circuits as circuits
import sympleq.core.paulis as paulis
