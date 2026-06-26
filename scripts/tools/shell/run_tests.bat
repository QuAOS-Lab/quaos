@echo off

REM Change the working directory to the script's
REM directory and load environment variables
cd /d %~dp0
call env.bat
cd %PROJECT_ROOT%

if not exist %SRC_VENV% (
    echo Virtual environment not found. Please, setup environment first.
    exit /b 1
)

call %SRC_VENV%/Scripts/activate

REM Build pytest command with optional markers
set PYTEST_CMD=pytest --override-ini "pytest.ini=%PYTEST_INI%" --cov=%PRJ_NAME% --cov-report=xml:%COVERAGE_REPORT_XML% --cov-report=html:%COVERAGE_REPORT_HTML% --junitxml=%COVERAGE_REPORT_JUNIT% --disable-warnings -vv

REM Add markers if provided as arguments
if not "%~1"=="" (
    set PYTEST_CMD=%PYTEST_CMD% -m "%*"
)

%PYTEST_CMD%
start "" "%COVERAGE_REPORT_HTML%\index.html"