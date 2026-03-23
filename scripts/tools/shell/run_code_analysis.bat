@echo off

REM Change the working directory to the script's
REM directory and load environment variables
cd /d %~dp0
call env.bat
cd %PROJECT_ROOT%

REM Check if the virtual environment exists
if not exist "%SRC_VENV%" (
    echo Virtual environment not found.
    exit /b 1
)

call %SRC_VENV%/Scripts/activate >nul 2>&1

if not exist "%DEV_REQUIREMENTS%" (
    echo %DEV_REQUIREMENTS% not found.
    exit /b 1
)

call pip install -r %DEV_REQUIREMENTS% >nul 2>&1

REM Check if accuracy is provided, default to 60
if "%~1"=="" (
    set ACCURACY=60
) else (
    set ACCURACY=%~1
)

echo.
echo Running Vulture code analysis with %ACCURACY% confidence...
echo.

REM Run Vulture to find unused code
vulture "%SRC_ROOT%" --min-confidence %ACCURACY%
