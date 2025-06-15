#!/bin/bash

cd "$(dirname "$0")"
cd ../..

if [ ! -d "venv" ]; then
    echo "Virtual environment not found."
    exit 1
fi

if [ ! -f "./configs/requirements.txt" ]; then
    echo "requirements.txt not found."
    exit 1
fi

source venv/bin/activate
pip freeze > ./configs/current_requirements.txt
python scripts/python/sync_requirements.py --r ./configs/requirements.txt

cd configs
mv -f updated_requirements.txt requirements.txt
rm -f current_requirements.txt
echo "Requirements have been synced and updated."
