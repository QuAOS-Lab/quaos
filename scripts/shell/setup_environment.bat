@echo off

cd /d %~dp0
cd ../..

REM Initializing virtual environment...

if not exist venv (
    echo Creating virtual environment venv...
    python -m venv venv
)

if not exist "configs/requirements.txt" (
    echo requirements.txt not found.
) else (
    call venv/Scripts/activate
    call pip install -r ./configs/requirements.txt
    call deactivate
)


REM Generating unversioned folders...

set "folders=data data/profiling scripts/personal"

for %%F in (%folders%) do (
    if not exist "%%F" (
        mkdir "%%F"
        echo Created folder: %%F
    )
)


REM Writing unversioned .env file...

if not exist "./configs/.env" (
    echo GITHUB_TOKEN="undefined" > ./configs/.env
)
echo "Please fill the .env file with the required secrets if not already done."