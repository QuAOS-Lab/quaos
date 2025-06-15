@echo off

cd /d %~dp0
cd ../..

if not exist "venv" (
    echo Virtual environment not found.
    exit /b 1
)

if not exist "./configs/requirements.txt" (
    echo requirements.txt not found.
    exit /b 1
)

call venv/Scripts/activate
call pip freeze > ./configs/current_requirements.txt
call python scripts/python/sync_requirements.py  --r ./configs/requirements.txt

cd configs
move /y updated_requirements.txt requirements.txt
del current_requirements.txt
echo Requirements have been synced and updated.