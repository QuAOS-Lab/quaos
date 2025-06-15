@echo off

cd /d %~dp0
cd ../../..

REM Check if the virtual environment exists
if not exist "venv" (
    echo Virtual environment not found.
    exit /b 1
)

echo Activating virtual environment...
call venv/Scripts/activate

if not exist "./configs/requirements/dev_requirements.txt" (
    echo dev_requirements.txt not found.
    exit /b 1
)

echo Installing requirements from dev_requirements.txt...
call pip install -r ./configs/requirements/dev_requirements.txt >nul 2>&1

echo Clearing Jupyter notebooks...
call python ./scripts/tools/python/clear_notebooks.py

echo Done!