Coding Standards
================

PEP 8 Style Guide
-----------------

The project follows **PEP 8**, the official Python style guide, to maintain readability and consistency across the codebase.

**Key Guidelines:**

- **Indentation**: 4 spaces per indentation level (no tabs)
- **Line Length**: Maximum 79 characters for code, 72 for docstrings/comments
- **Naming Conventions**:
  
  - Functions and variables: ``snake_case``
  - Classes: ``PascalCase``
  - Constants: ``UPPER_CASE``
  - Private attributes: ``_leading_underscore``

- **Imports**: Grouped and ordered (standard library, third-party, local)
- **Whitespace**: Follow PEP 8 rules for operators and commas

**Resources:**

- `PEP 8 Style Guide <https://www.python.org/dev/peps/pep-0008/>`_
- `PEP 8 Quick Reference <https://pep8.org/>`_

.. tip::
   The recommended VS Code extensions (see :doc:`extensions`) automatically detect PEP 8 violations and provide inline suggestions for fixes.


Docstring Standards
-------------------

All code must include comprehensive docstrings following the **reStructuredText (RST) style** for seamless integration with Sphinx documentation.

**Requirements:**

- All public modules, classes, methods, and functions must have docstrings
- Use RST directives (``:param:``, ``:type:``, ``:return:``, ``:raises:``)
- Include usage examples where appropriate
- Document exceptions and edge cases

See :doc:`documentation` for detailed docstring examples and formatting guidelines.


Type Hints
----------

Use Python type hints for all function signatures:

.. code-block:: python

   def process_pauli(pauli: Pauli, coefficient: complex = 1.0) -> PauliSum:
       """Process a Pauli operator with optional coefficient."""
       return PauliSum({pauli: coefficient})

**Benefits:**

- Improved code readability
- Enhanced IDE autocomplete and error detection
- Better documentation generation
- Runtime type checking support with tools like ``mypy``
