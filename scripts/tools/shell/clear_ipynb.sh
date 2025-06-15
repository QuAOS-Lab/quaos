#!/bin/bash

cd "$(dirname "$0")"
cd ../../..

# Check if the virtual environment exists
if [ ! -d "venv" ]; then
    echo "Virtual environment not found."
    exit 1
fi

echo "Activating virtual environment..."
source venv/bin/activate

if [ ! -f "./configs/requirements/dev_requirements.txt" ]; then
    echo "dev_requirements.txt not found."
    deactivate 2>/dev/null
    exit 1
fi

echo "Installing requirements from dev_requirements.txt..."
pip install -r ./configs/requirements/dev_requirements.txt >/dev/null 2>&1

echo "Clearing Jupyter notebooks..."
python ./scripts/tools/python/clear_notebooks.py

echo "Done!"