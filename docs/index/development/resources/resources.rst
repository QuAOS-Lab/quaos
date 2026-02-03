.. _development_resources:

=========
Resources
=========

This page provides a comprehensive collection of resources to help you with the development of SympleQ and the surrounding ecosystem. It covers project structure, development tools, workflows, and coding standards.


Project Structure Overview
===========================

The project is structured as shown below. Each folder serves a specific purpose which is described in the corresponding section.

.. code-block:: text

    sympleq/
    ├── .github/
    ├── .vscode/
    ├── docs/
    ├── src/
    │   └── sympleq/
    ├── scripts/
    │   ├── configs/
    │   ├── examples/
    │   ├── experiments/
    │   ├── notebooks/
    │   ├── personal/
    │   ├── tools/
    │   └── TO_BE_ORGANIZED/
    └── tests/


.github/
--------

Contains GitHub workflows for continuous integration and deployment (CI/CD). Structured following the `GitHub Actions <https://docs.github.com/en/actions/writing-workflows/quickstart>`_ convention.

**Key Components:**

- ``workflows/``: CI/CD pipeline definitions
  
  - ``check_source.yml``: Code quality validation and linting
  - ``review_code.yml``: Test execution and coverage reporting
  - ``tag_commit.yml``: Automated SemVer tagging on ``dev`` branch
  - ``publish_repo.yml``: Package building and PyPI publication

- ``actions/``: Custom reusable GitHub Actions

.. code-block:: text

    .github/
    ├── actions/
    └── workflows/
        ├── check_source.yml
        ├── review_code.yml
        ├── tag_commit.yml
        └── publish_repo.yml


.vscode/
--------

Visual Studio Code workspace configuration. Contains IDE settings, task definitions, and extension recommendations.

**Key Files:**

- ``tasks.json``: Automated development tasks (build docs, run tests, profile code)
- ``extensions.json``: Recommended VS Code extensions
- ``settings.json``: Workspace-specific editor settings

.. code-block:: text

    .vscode/
    ├── tasks.json
    ├── extensions.json
    └── settings.json

.. note::
   While typically excluded from version control, this folder is included to streamline IDE setup and provide a consistent development environment across the team.


docs/
-----

Source files for Sphinx-generated documentation.

**Structure:**

- ``conf.py``: Sphinx configuration
- ``index.rst``: Documentation entry point
- ``index/``: Documentation content organized by topic
- ``_build/``: Generated HTML output (not version controlled)
- ``_static/``: Custom CSS and static assets
- ``_templates/``: Custom Sphinx templates
- ``doc_venv/``: Dedicated virtual environment for documentation builds

.. code-block:: text

    docs/
    ├── conf.py
    ├── index.rst
    ├── doc_requirements.txt
    ├── index/
    │   ├── about/
    │   ├── development/
    │   │   ├── resources/
    │   │   └── contributing/
    │   └── guides/
    ├── _build/
    ├── _static/
    ├── _templates/
    └── doc_venv/


src/sympleq/
------------

Main source code for the SympleQ package. Contains all production code organized into logical modules.

.. code-block:: text

    src/
    └── sympleq/
        ├── ...


scripts/
--------

Development scripts, examples, experiments, and auxiliary tools.

**Sub-directories:**

``configs/``
   Configuration files for development tools:
   
   - ``dev_requirements.txt``: Development dependencies
   - ``requirements.txt``: Production dependencies
   - ``pytest.ini``: Pytest configuration
   - ``pytest_cov_config.ini``: Coverage settings
   - ``setup.cfg``: Package metadata and tool configurations
   - ``cspell.json``: Spell checker configuration

``examples/``
   Example scripts demonstrating library usage:
   
   - ``generate_basis.py``: Pauli basis generation
   - ``pauli_operations.py``: Pauli algebra operations
   - ``pauli_symmetries.py``: Symmetry detection examples
   - ``notebooks/``: Example Jupyter notebooks

``experiments/``
   Research experiments and exploratory code (not part of production package)

``notebooks/``
   Jupyter notebooks for development and analysis

``personal/``
   Personal workspace for individual developers (not version controlled)

``tools/``
   Automation scripts for development tasks:
   
   - ``python/``: Python utility scripts (profiling, benchmark listing, notebook cleaning)
   - ``shell/``: Shell scripts for Windows (``.bat``) and Linux (``.sh``)

``TO_BE_ORGANIZED/``
   Temporary storage for code awaiting organization

.. code-block:: text

    scripts/
    ├── configs/
    │   ├── dev_requirements.txt
    │   ├── requirements.txt
    │   ├── pytest.ini
    │   └── setup.cfg
    ├── examples/
    │   └── notebooks/
    ├── experiments/
    ├── notebooks/
    ├── personal/
    │   └── profiling/
    ├── tools/
    │   ├── python/
    │   │   ├── profiling_script.py
    │   │   ├── list_benchmarks.py
    │   │   └── clear_notebooks.py
    │   └── shell/
    │       ├── env.bat
    │       ├── env.sh
    │       ├── docs.bat
    │       ├── docs.sh
    │       ├── run_tests.bat
    │       ├── run_tests.sh
    │       ├── run_profiler.bat
    │       └── run_profiler.sh
    └── TO_BE_ORGANIZED/


tests/
------

Test suite for the SympleQ package. Organized to mirror the source code structure.

**Components:**

- Unit tests for each module
- Integration tests
- Benchmark tests (marked with ``@pytest.mark.benchmark``)
- Test data and fixtures in ``_test_data/``

Additional Resources
====================

External Documentation
----------------------

- `Python Documentation <https://docs.python.org/3/>`_
- `NumPy Documentation <https://numpy.org/doc/stable/>`_
- `Pytest Documentation <https://docs.pytest.org/>`_
- `Sphinx Documentation <https://www.sphinx-doc.org/>`_
- `Git Documentation <https://git-scm.com/doc>`_
- `VS Code Documentation <https://code.visualstudio.com/docs>`_

Community
---------

- GitHub Discussions: Ask questions and share ideas
- Issue Tracker: Report bugs and request features
- Code Owners: See ``.github/CODEOWNERS`` for module maintainers


Getting Help
============

If you encounter issues or have questions:

1. **Check Documentation**: Search this documentation for relevant guides
2. **Review Issues**: Check if your question has been answered in GitHub Issues
3. **Ask the Team**: Reach out to code owners for specific modules
4. **Open an Issue**: Create a new issue with detailed description and reproduction steps

.. note::
   For development environment issues, ensure you've followed all steps in :doc:`tools` and run the **Setup Dev Environment** task.


Development Guides
==================

.. toctree::
   :maxdepth: 1

   tools
   extensions
   version_control
   tasks
   documentation
   coding_standards

The following guides provide step-by-step instructions for setting up your development environment and understanding project workflows:

1. :doc:`tools` - Installation and configuration of required development tools (Python, Git, VS Code)
2. :doc:`extensions` - Essential VS Code extensions for development
3. :doc:`version_control` - Git workflow, branching strategy, and versioning
4. :doc:`tasks` - Automated development tasks (testing, documentation, profiling)
5. :doc:`documentation` - Writing and building project documentation
6. :doc:`coding_standards` - Code style, documentation standards, and best practices
