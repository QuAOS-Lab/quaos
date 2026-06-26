@echo off

REM Change the working directory to the script's
REM directory and load environment variables
cd /d %~dp0
call env.bat
cd %PROJECT_ROOT%

REM Initializing virtual environment...
if exist %SRC_VENV% (rmdir /s /q %SRC_VENV% 2>nul)
    echo Creating virtual environment %SRC_VENV%...
    powershell -NonInteractive -Command "py -m venv '%SRC_VENV%'"
    if not exist %SRC_VENV%/Scripts/activate.bat (
        echo ERROR: Failed to create virtual environment.
        echo Ensure Python is installed, then run this task again.
        exit /b 1
    )
)

call %SRC_VENV%/Scripts/activate.bat
call python -m pip install uv
call uv pip install -e ".[development]"

REM Optional package groups
echo.
choice /c YN /m "Install quantinuum packages (pytket, pytket-quantinuum, qnexus)?"
if errorlevel 2 goto skip_quantinuum
call uv pip install -e ".[quantinuum]"
:skip_quantinuum

echo.
choice /c YN /m "Install RBMalgorithms packages (torch, aepsych, botorch, gpytorch)?"
if errorlevel 2 goto skip_rbmalgorithms
call uv pip install -e ".[RBMalgorithms]"
:skip_rbmalgorithms

call deactivate

REM Generating unversioned folders...
set "folders=%PERSONAL_FOLDER%"
for %%F in (%folders%) do (
    if not exist "%%F" (
        mkdir "%%F"
        echo Created folder: %%F
    )
)
