#!/bin/bash

cd "$(dirname "$0")/../../docs"

if [ ! -d "doc_venv" ]; then
    echo "Creating virtual environment doc_venv..."
    python3 -m venv doc_venv
fi

source doc_venv/bin/activate
pip install -r doc_requirements.txt

pip install -r ../configs/requirements.txt

rm -rf "index/_autosummary/"

sphinx-build -E -b html . _build/html

if [ $? -ne 0 ]; then
    echo "Sphinx build failed"
    exit 1
fi

xdg-open _build/html/index.html || open _build/html/index.html