#!/bin/bash

# Change to the script's directory and load environment variables
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
source env.sh
cd "$PROJECT_ROOT"

# Check if the virtual environment exists
if [ ! -d "$SRC_VENV" ]; then
    echo "Virtual environment not found."
    exit 1
fi

# Activate the virtual environment
source "$SRC_VENV/bin/activate"

# Check if accuracy is provided, default to 60
ACCURACY=${1:-60}

echo
echo "Running Vulture code analysis with $ACCURACY% confidence..."
echo

# Run Vulture to find unused code
vulture "$SRC_ROOT" --min-confidence "$ACCURACY"
