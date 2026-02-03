#!/bin/bash

# Change to the script's directory and load environment variables
cd "$(dirname "$0")"
source env.sh
cd "$PROJECT_ROOT"

if [ ! -d "$SRC_VENV" ]; then
    echo "Virtual environment not found."
    exit 1
fi

source "$SRC_VENV/bin/activate" > /dev/null 2>&1

if [ ! -f "$DEV_REQUIREMENTS" ]; then
    echo "$DEV_REQUIREMENTS not found."
    exit 1
fi

pip install -r "$DEV_REQUIREMENTS" > /dev/null 2>&1

if [ -z "$1" ]; then
    echo "Missing required argument!"
    echo "Select one from the list below or launch:"
    echo "  run_profiler.sh <benchmark_function_signature>"
    echo

    python scripts/tools/python/list_benchmarks.py "$PROJECT_ROOT" "$PROFILING_OUTPUT_PATH"
    
    if [ $? -ne 0 ]; then
        echo
        exit 1
    fi

    TEST_NAME_FILE="$PROFILING_OUTPUT_PATH/test_to_profile.txt"
    if [ -f "$TEST_NAME_FILE" ]; then
        TEST_NAME=$(cat "$TEST_NAME_FILE")
    else
        echo "No test selected."
        exit 1
    fi
else
    TEST_NAME="$1"
fi

echo
echo "Profiling test: $TEST_NAME"
echo

python scripts/tools/python/profiling_script.py "$PROJECT_ROOT" "$PROFILING_OUTPUT_PATH" "$TEST_NAME"
