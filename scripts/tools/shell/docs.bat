@echo off

REM Change the working directory to the script's
REM directory and load environment variables
cd /d %~dp0
call env.bat
cd %PROJECT_ROOT%

if not exist %DOC_VENV% (
    echo Creating virtual environment in %DOC_VENV%...
    python -m venv %DOC_VENV%
)

call %DOC_VENV%\Scripts\activate
pip install -r %DOC_REQUIREMENTS%
pip install -e .

rmdir /s /q "%DOC_AUTOSUMMARY%"

sphinx-build -E -b html %DOC_ROOT% %DOC_BUILD_DIR%

if %errorlevel% neq 0 (
    echo Sphinx build failed
    exit /b %errorlevel%
)

start "" "%DOC_INDEX%"
call deactivate