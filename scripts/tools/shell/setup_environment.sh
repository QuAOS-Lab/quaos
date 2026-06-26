#!/bin/bash

# Change the working directory to the script's
# directory and load environment variables
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
source env.sh
cd "$PROJECT_ROOT"

# Initializing virtual environment...
if [[ ! -d "$SRC_VENV" ]]; then
    echo "Creating virtual environment $SRC_VENV..."
    python -m venv "$SRC_VENV"
fi

source "$SRC_VENV/bin/activate"
python -m pip install uv
uv pip install -e "${PYTHON_PY_SETUP}[development]"

# Optional package groups
echo
read -rp "Install quantinuum packages (pytket, pytket-quantinuum, qnexus)? [Y/n]: " yn
case "$yn" in
    [Yy]*) uv pip install -e ".[quantinuum]" ;;
esac

echo
read -rp "Install RBMalgorithms packages (torch, aepsych, botorch, gpytorch)? [Y/n]: " yn
case "$yn" in
    [Yy]*) uv pip install -e ".[RBMalgorithms]" ;;
esac

deactivate

# Generating unversioned folders...
folders="$PERSONAL_FOLDER"
for folder in $folders; do
    if [[ ! -d "$folder" ]]; then
        mkdir -p "$folder"
        echo "Created folder: $folder"
    fi
done
