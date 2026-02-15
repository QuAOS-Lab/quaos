#!/usr/bin/env bash

# Change to the script's directory and load environment variables
cd "$(dirname "$0")"
source env.sh
cd "$PROJECT_ROOT"

if [ ! -d "$DOC_VENV" ]; then
    echo "Creating virtual environment in $DOC_VENV..."
    python3 -m venv "$DOC_VENV"
fi

source "$DOC_VENV/bin/activate"
pip install -r "$DOC_REQUIREMENTS"
pip install -e .

rm -rf "$DOC_AUTOSUMMARY"

sphinx-build -E -b html "$DOC_ROOT" "$DOC_BUILD_DIR"
if [ $? -ne 0 ]; then
    echo "Sphinx build failed"
    exit $?
fi

# Open the built documentation index in the default browser
xdg-open $DOC_INDEX 2>/dev/null || echo "Docs built at $DOC_INDEX"
# Wait for 2 seconds to ensure the browser opens
sleep 2
deactivate