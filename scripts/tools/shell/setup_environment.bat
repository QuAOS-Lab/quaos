@echo off

REM Change the working directory to the script's
REM directory and load environment variables
cd /d %~dp0
call env.bat
cd %PROJECT_ROOT%

CALL %VCVARSALL_PATH% x64

REM Initializing virtual environment...
if not exist %SRC_VENV% (
    echo Creating virtual environment %SRC_VENV%...
    python -m venv %SRC_VENV%
    call %SRC_VENV%/Scripts/activate
    call python -m pip install -r %DEV_REQUIREMENTS%
    call deactivate
)

set ARCH=%PROCESSOR_ARCHITECTURE%
REM On 64-bit Windows with 32-bit CMD, check PROCESSOR_ARCHITEW6432
if "%ARCH%"=="x86" (
    if not "%PROCESSOR_ARCHITEW6432%"=="" (
        set ARCH=%PROCESSOR_ARCHITEW6432%
    )
)
if /I "%ARCH%"=="AMD64" (
    set PRESET=%CMAKE_X86_64_PRESET%
) else if /I "%ARCH%"=="ARM64" (
    echo Architecture %SRC_VENV% unsupported on Windows os
    exit -1
) else (
    echo Architecture %SRC_VENV% unsupported on Windows os
    exit -1
)
echo Detected architecture %ARCH%: using preset %PRESET%

call %SRC_VENV%/Scripts/activate
call python -m pip install --upgrade pip setuptools wheel scikit-build-core nanobind setuptools-scm
call python -m pip install -e %PYTHON_PY_SETUP% --config-settings=cmake.args="--preset;%PRESET%"
call deactivate

REM Generating unversioned folders...
set "folders=%PERSONAL_FOLDER%"
for %%F in (%folders%) do (
    if not exist "%%F" (
        mkdir "%%F"
        echo Created folder: %%F
    )
)
