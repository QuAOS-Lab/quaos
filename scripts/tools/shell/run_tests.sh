#!/bin/bash

# Change the working directory to the script's directory
cd "$(dirname "$0")"
source env.sh
cd "$PROJECT_ROOT"

if [ ! -d "$SRC_VENV" ]; then
    echo "Virtual environment not found. Please, setup environment first."
    exit 1
fi

source "$SRC_VENV/bin/activate"

# Build pytest command with optional markers
PYTEST_CMD="pytest --override-ini pytest.ini=$PYTEST_INI --cov=$PRJ_NAME --cov-report=xml:$COVERAGE_REPORT_XML --cov-report=html:$COVERAGE_REPORT_HTML --junitxml=$COVERAGE_REPORT_JUNIT --disable-warnings -vv"

# Add markers if provided as arguments
if [ $# -gt 0 ]; then
    PYTEST_CMD="$PYTEST_CMD -m \"$*\""
fi

eval $PYTEST_CMD
xdg-open "$COVERAGE_REPORT_HTML/index.html" 2>/dev/null || open "$COVERAGE_REPORT_HTML/index.html" 2>/dev/null || start "$COVERAGE_REPORT_HTML/index.html"
