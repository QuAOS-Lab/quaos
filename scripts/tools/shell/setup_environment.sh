#!/usr/bin/env bash

# Change to the script's directory and load environment variables
cd "$(dirname "$0")"
source env.sh
cd "$PROJECT_ROOT"

# Initializing virtual environment...
if [ ! -d "$SRC_VENV" ]; then
    echo "Creating virtual environment $SRC_VENV..."
    python3 -m venv "$SRC_VENV"
    source "$SRC_VENV/bin/activate"
    python3 -m pip install -r "$DEV_REQUIREMENTS"
    deactivate
fi

ARCH=$(uname -m)

if [ "$ARCH" == "x86_64" ]; then
    PRESET=$CMAKE_X86_64_PRESET
elif [ "$ARCH" == "aarch64" ]; then
    echo "Architecture $SRC_VENV unsupported on Linux OS"
    exit 1
else
    echo "Architecture $SRC_VENV unsupported on Linux OS"
    exit 1
fi

echo "Detected architecture $ARCH: using preset $PRESET"

source "$SRC_VENV/bin/activate"
python3 -m pip install  --upgrade pip setuptools wheel scikit-build-core nanobind setuptools-scm
python3 -m pip install -e "$PYTHON_PY_SETUP" --config-settings=cmake.args="--preset;$PRESET"
deactivate

# Generating unversioned folders...
folders=("$PERSONAL_FOLDER")
for F in "${folders[@]}"; do
    if [ ! -d "$F" ]; then
        mkdir -p "$F"
        echo "Created folder: $F"
    fi
done
