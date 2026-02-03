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

REM Check if test function name is provided
if "%~1"=="" (
    echo Missing required argument! 
    echo Select one from the list below or launch:
    echo   run_profiler.bat ^<benchmark_function_signature^>
    echo.
    python scripts/tools/python/list_benchmarks.py %PROJECT_ROOT% %PROFILING_OUTPUT_PATH%
    
    if errorlevel 1 (
        echo.
        exit /b 1
    )

    set /p TEST_NAME=< "%PROFILING_OUTPUT_PATH%\test_to_profile.txt"
) else (
    set TEST_NAME=%~1
)

echo.
echo Profiling test: %TEST_NAME%
echo.

REM Run the profiling script
call python scripts/tools/python/profiling_script.py %PROJECT_ROOT% %PROFILING_OUTPUT_PATH% %TEST_NAME%
