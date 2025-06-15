@echo off

cd /d %~dp0
cd ../../docs

if not exist doc_venv (
    echo Creating virtual environment doc_venv...
    python -m venv doc_venv
)

call doc_venv\Scripts\activate
pip install -r doc_requirements.txt

pip install -r ../configs/requirements.txt

rmdir /s /q "index/_autosummary/"

sphinx-build -E -b html . _build/html

if %errorlevel% neq 0 (
    echo Sphinx build failed
    exit /b %errorlevel%
)

start "" "%cd%\_build\html\index.html"
