#!/usr/bin/env bash

# Change to the script's directory and load environment variables
cd "$(dirname "$0")"
source env.sh
cd "$PROJECT_ROOT"

# Check if the virtual environment exists
if [ ! -d "$SRC_VENV" ]; then
    echo "Virtual environment not found."
    exit 1
fi

echo "Activating virtual environment..."
source "$SRC_VENV/bin/activate"

echo "Clearing Jupyter notebooks..."
python3 "$CLEAR_NOTEBOOKS_SCRIPT" "$NOTEBOOKS_ROOT_DIR"

echo "Done!"
deactivate