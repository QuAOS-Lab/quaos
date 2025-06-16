#!/bin/bash

# Move to the script's directory and then go two levels up
cd "$(dirname "$0")"
cd ../..

# Initializing virtual environment
if [ ! -d "venv" ]; then
    echo "Creating virtual environment venv..."
    python3 -m venv venv
fi

# Installing dependencies if requirements.txt exists
if [ ! -f "configs/requirements.txt" ]; then
    echo "requirements.txt not found."
else
    source venv/bin/activate
    pip install -r configs/requirements.txt
    deactivate
fi

# Generating unversioned folders
folders=("data" "data/profiling" "scripts/personal")

for folder in "${folders[@]}"; do
    if [ ! -d "$folder" ]; then
        mkdir -p "$folder"
        echo "Created folder: $folder"
    fi
done

# Writing unversioned .env file
if [ ! -f "configs/.env" ]; then
    echo 'GITHUB_TOKEN="undefined"' > configs/.env
fi
echo "Please fill the .env file with the required secrets if not already done."

