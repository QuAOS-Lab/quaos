# All the variables below assume PROJECT_ROOT as working directory
export PROJECT_ROOT="$(realpath "$(dirname "$0")/../../..")"
export PRJ_NAME=sympleq

# Project env variables
export CONFIGS_ROOT=scripts/configs
export SRC_ROOT=src
export DOC_ROOT=docs
export SRC_VENV=venv

# Docs env variables
export DOC_VENV=$DOC_ROOT/doc_venv
export DOC_REQUIREMENTS=$DOC_ROOT/doc_requirements.txt
export DOC_AUTOSUMMARY=$DOC_ROOT/index/guides/_autosummary
export DOC_BUILD_DIR=$DOC_ROOT/_build/html
export DOC_INDEX=$DOC_ROOT/_build/html/index.html

# Setup dev environment env variables
export PYTHON_PY_SETUP="./"
export PERSONAL_FOLDER=scripts/personal
export DIST_FOLDER=dist/
export DEV_REQUIREMENTS=scripts/configs/dev_requirements.txt

# Clear notebooks env variables
export NOTEBOOKS_ROOT_DIR=$PROJECT_ROOT
export CLEAR_NOTEBOOKS_SCRIPT=scripts/tools/python/clear_notebooks.py

# Test coverage
export PYTEST_INI=$PROJECT_ROOT/scripts/configs/pytest.ini
export COVERAGE_REPORT_XML=$PROJECT_ROOT/scripts/personal/coverage_report.xml
export COVERAGE_REPORT_HTML=$PROJECT_ROOT/scripts/personal/coverage_html_report
export COVERAGE_REPORT_JUNIT=$PROJECT_ROOT/scripts/personal/junit.xml

# Profiling
export PROFILING_OUTPUT_PATH=$PROJECT_ROOT/scripts/personal/profiling