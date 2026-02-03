.. _development_tasks:

======================
Development Tasks
======================

This document describes the automated development tasks available in the project through VS Code's task runner.

Task Runner Extension
======================

The project uses the **Task Runner Code** extension (`sanaajani.taskrunnercode`) to provide a convenient UI for running tasks directly from the VS Code sidebar. This extension displays all tasks defined in `.vscode/tasks.json` in a dedicated Task Runner panel.

.. note::
   Install the recommended extension by opening the Command Palette (Ctrl+Shift+P) and selecting "Extensions: Show Recommended Extensions".

Available Tasks
===============

All tasks are defined in `.vscode/tasks.json` and have corresponding shell scripts in `scripts/tools/shell/`. Each task supports both Windows (`.bat`) and Linux (`.sh`) environments.

Build Sphinx Documentation
---------------------------

**Label:** Build Sphinx Documentation

**Purpose:** Generates HTML documentation from reStructuredText source files using Sphinx.

**Script Locations:**
   - Windows: `scripts/tools/shell/docs.bat`
   - Linux: `scripts/tools/shell/docs.sh`

**Workflow:**
   1. Creates a dedicated virtual environment in `docs/doc_venv/` if it doesn't exist
   2. Installs documentation requirements from `docs/doc_requirements.txt`
   3. Removes previous autosummary outputs
   4. Runs Sphinx build with `-E` (rebuild all) and `-b html` (HTML output)
   5. Automatically opens the generated documentation in the default browser

**Output:** HTML documentation in `docs/_build/html/`

**Environment Variables Used:**
   - ``DOC_VENV``: Path to documentation virtual environment
   - ``DOC_REQUIREMENTS``: Path to documentation requirements file
   - ``DOC_ROOT``: Documentation source directory
   - ``DOC_BUILD_DIR``: Output directory for built documentation
   - ``DOC_INDEX``: Path to index.html

Setup Dev Environment
---------------------

**Label:** Setup Dev Environment

**Purpose:** Initializes the development environment with all necessary dependencies.

**Script Locations:**
   - Windows: `scripts/tools/shell/setup_environment.bat`
   - Linux: `scripts/tools/shell/setup_environment.sh`

**Workflow:**
   1. Creates main project virtual environment in `venv/` if it doesn't exist
   2. Upgrades pip, setuptools, and setuptools-scm
   3. Installs the project package in editable mode (`pip install -e`)
   4. Creates unversioned folders (e.g., `scripts/personal/`)

**Environment Variables Used:**
   - ``SRC_VENV``: Path to source virtual environment
   - ``DEV_REQUIREMENTS``: Path to development requirements file
   - ``PYTHON_PY_SETUP``: Path to project root for editable install
   - ``PERSONAL_FOLDER``: Path to personal scripts folder

.. important::
   Run this task first before using any other development tasks.

Clear Notebooks Content
------------------------

**Label:** Clear Notebooks Content

**Purpose:** Clears output cells and execution counts from Jupyter notebooks for clean version control.

**Script Locations:**
   - Windows: `scripts/tools/shell/clear_ipynb.bat`
   - Linux: `scripts/tools/shell/clear_ipynb.sh`

**Workflow:**
   1. Activates the project virtual environment
   2. Installs required dependencies
   3. Runs the notebook cleaning script on all notebooks in ``NOTEBOOKS_ROOT_DIR``

**Environment Variables Used:**
   - ``SRC_VENV``: Path to source virtual environment
   - ``DEV_REQUIREMENTS``: Path to development requirements
   - ``CLEAR_NOTEBOOKS_SCRIPT``: Path to Python notebook cleaning script
   - ``NOTEBOOKS_ROOT_DIR``: Root directory containing notebooks to clean

.. tip::
   Use this task before committing notebooks to ensure clean diffs.

Run Tests With Coverage
------------------------

**Label:** Run Tests With Coverage

**Purpose:** Executes pytest test suite with coverage reporting.

**Script Locations:**
   - Windows: `scripts/tools/shell/run_tests.bat`
   - Linux: `scripts/tools/shell/run_tests.sh`

**Workflow:**
   1. Activates the project virtual environment
   2. Installs development requirements
   3. Runs pytest with coverage enabled
   4. Generates XML, HTML, and JUnit coverage reports
   5. Opens HTML coverage report in the default browser

**Command Line Arguments:**
   Accepts optional pytest markers as arguments to filter tests:
   
   .. code-block:: bash
   
      run_tests.bat benchmark           # Run only benchmark tests
      run_tests.bat "not acceptance"    # Skip acceptance tests

**Environment Variables Used:**
   - ``SRC_VENV``: Path to source virtual environment
   - ``DEV_REQUIREMENTS``: Path to development requirements
   - ``PYTEST_INI``: Path to pytest configuration file
   - ``PRJ_NAME``: Project name for coverage measurement
   - ``COVERAGE_REPORT_XML``: Path for XML coverage report
   - ``COVERAGE_REPORT_HTML``: Path for HTML coverage report
   - ``COVERAGE_REPORT_JUNIT``: Path for JUnit XML report

**Output:** Coverage reports in `scripts/personal/htmlcov/`

Profile Code
------------

**Label:** Profile Code

**Purpose:** Profiles benchmark tests using cProfile and visualizes results with snakeviz.

**Script Locations:**
   - Windows: `scripts/tools/shell/run_profiler.bat`
   - Linux: `scripts/tools/shell/run_profiler.sh`

**Workflow:**
   1. Activates the project virtual environment
   2. If no test name provided as argument:
      
      a. Lists all available benchmark tests using pytest collection
      b. Prompts user for interactive selection
      c. Stores selection in temporary file
   
   3. Runs the selected test directly (bypassing pytest overhead) with cProfile
   4. Generates `.prof` profiling output
   5. Launches snakeviz for interactive visualization

**Command Line Arguments:**
   Accepts optional test function signature:
   
   .. code-block:: bash
   
      run_profiler.bat tests/benchmark_test.py::test_paulisum_sum

**Environment Variables Used:**
   - ``SRC_VENV``: Path to source virtual environment
   - ``DEV_REQUIREMENTS``: Path to development requirements
   - ``PROFILING_OUTPUT_PATH``: Output directory for profiling results

**Output:** Profiling data in `scripts/personal/profiling/`

.. note::
   **Interactive Input Handling**
   
   The Profile Code task requires special handling for interactive input in VS Code's integrated terminal. The Python script ``list_benchmarks.py`` uses platform-specific methods to clear the stdin buffer before prompting for input:
   
   - **Windows:** Uses ``msvcrt.kbhit()`` and ``msvcrt.getch()``
   - **Linux/Unix:** Uses ``termios.tcflush()`` with ``TCIFLUSH``
   
   This prevents terminal echo commands (like virtual environment activation) from being inadvertently consumed as user input.

Environment Configuration
==========================

All tasks load environment variables from centralized configuration files:

- ``scripts/tools/shell/env.bat`` (Windows)
- ``scripts/tools/shell/env.sh`` (Linux/Unix)

These files define consistent paths and project settings used across all automation scripts.

Key Variables
-------------

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Variable
     - Description
   * - ``PROJECT_ROOT``
     - Absolute path to project root directory
   * - ``SRC_VENV``
     - Path to main development virtual environment
   * - ``DOC_VENV``
     - Path to documentation virtual environment
   * - ``DEV_REQUIREMENTS``
     - Path to development dependencies file
   * - ``DOC_REQUIREMENTS``
     - Path to documentation dependencies file
   * - ``PROFILING_OUTPUT_PATH``
     - Output directory for profiling results
   * - ``NOTEBOOKS_ROOT_DIR``
     - Root directory containing Jupyter notebooks
   * - ``PRJ_NAME``
     - Project name (used in coverage and imports)

Platform Support
================

.. list-table::
   :widths: 40 20 20 20
   :header-rows: 1

   * - Task
     - Windows
     - Linux
     - macOS
   * - Build Sphinx Documentation
     - ✓
     - ✓
     - ✗
   * - Setup Dev Environment
     - ✓
     - ✓
     - ✗
   * - Clear Notebooks Content
     - ✓
     - ✓
     - ✗
   * - Run Tests With Coverage
     - ✓
     - ✓
     - ✗
   * - Profile Code
     - ✓
     - ✓
     - ✗

.. note::
   macOS support could be added by creating corresponding `.sh` scripts, but is currently not implemented.

Task Execution
==============

Tasks can be executed in three ways:

1. **Task Runner Panel:** Click tasks in the Task Runner Code extension sidebar
2. **Command Palette:** Press ``Ctrl+Shift+P``, type "Tasks: Run Task", and select from the list
3. **Keyboard Shortcut:** Configure custom keyboard shortcuts in VS Code settings
4. **Terminal:** Run shell scripts directly from command line

Example: Running from Terminal
-------------------------------

.. code-block:: bash

   # Windows
   cd scripts/tools/shell
   call run_tests.bat benchmark

   # Linux
   cd scripts/tools/shell
   bash run_tests.sh benchmark

Best Practices
==============

1. **Always run Setup Dev Environment first** when cloning the repository
2. **Clear notebook outputs** before committing to version control
3. **Use test markers** to run specific test subsets during development
4. **Profile specific tests** rather than entire suites for focused optimization
5. **Check coverage reports** to identify untested code paths

Troubleshooting
===============

Virtual Environment Not Found
------------------------------

If you see "Virtual environment not found" errors, run the **Setup Dev Environment** task to initialize the development environment.

Interactive Input Issues
------------------------

If the Profile Code task doesn't accept keyboard input properly:

1. Ensure you're using a standard terminal (not a specialized shell)
2. Try running the script directly from a command prompt/terminal
3. Check that the ``VSCODE_TASK_DISABLE_ENV_CONTRIBUTION`` environment variable is respected

Documentation Build Fails
--------------------------

If Sphinx documentation fails to build:

1. Check that all required reStructuredText files are present
2. Verify docstring formatting in source code
3. Review Sphinx warnings for syntax errors
4. Ensure the documentation virtual environment is properly configured
