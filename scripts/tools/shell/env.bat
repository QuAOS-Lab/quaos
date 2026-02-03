REM All the variable below assume PROJECT_ROOT as working directory
set PROJECT_ROOT=%~dp0../../..
set PRJ_NAME=sympleq

REM Project env variables
set CONFIGS_ROOT=scripts/configs
set SRC_ROOT=src
set DOC_ROOT=docs
set SRC_VENV=venv
set SRC_REQUIREMENTS=scripts/configs/requirements.txt

REM Docs env variables
set DOC_VENV=%DOC_ROOT%/doc_venv
set DOC_REQUIREMENTS=%DOC_ROOT%/doc_requirements.txt
set DOC_AUTOSUMMARY=%DOC_ROOT%/index/guides/_autosummary
set DOC_BUILD_DIR=%DOC_ROOT%/_build/html
set DOC_INDEX=%DOC_ROOT%/_build/html/index.html

REM Setup dev environment env variables
set PYTHON_PY_SETUP="./"
set PERSONAL_FOLDER=scripts/personal
set DIST_FOLDER=dist/
set DEV_REQUIREMENTS=scripts/configs/dev_requirements.txt

REM Clear notebooks env variables
set NOTEBOOKS_ROOT_DIR=%PROJECT_ROOT%
set CLEAR_NOTEBOOKS_SCRIPT=scripts/tools/python/clear_notebooks.py

REM Test coverage
set PYTEST_INI=%PROJECT_ROOT%/scripts/configs/pytest.ini
set COVERAGE_REPORT_XML=%PROJECT_ROOT%/scripts/personal/coverage_report.xml
set COVERAGE_REPORT_HTML=%PROJECT_ROOT%/scripts/personal/coverage_html_report
set COVERAGE_REPORT_JUNIT=%PROJECT_ROOT%/scripts/personal/junit.xml

REM Profiling
set PROFILING_OUTPUT_PATH=%PROJECT_ROOT%/scripts/personal/profiling