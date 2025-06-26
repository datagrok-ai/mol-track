@echo off
setlocal enabledelayedexpansion

set IMAGE_NAME=moltrack
set ENV_DIR=.moltrack-env

echo Building Docker image: %IMAGE_NAME% from Dockerfile
docker build -t %IMAGE_NAME% .

echo Looking for a free port
set "HOST_PORT="
for /L %%p in (5432,1,5500) do (
    netstat -an | findstr ":%%p " | findstr "LISTENING" >nul
    if errorlevel 1 (
        set "HOST_PORT=%%p"
        echo Found free port: %%p
        goto :found_port
    ) else (
        echo Port %%p is in use.
    )
)

echo No free port found
exit /b 1

:found_port
set DB_PORT=%HOST_PORT%
echo Set DB_PORT environment variable to %DB_PORT%

docker ps -aq -f "name=%IMAGE_NAME%" > temp_containers.txt
set /p EXISTING_CONTAINER=<temp_containers.txt
del temp_containers.txt

if not "%EXISTING_CONTAINER%"=="" (
    echo Stopping and removing existing container: %IMAGE_NAME%
    docker stop %IMAGE_NAME% >nul 2>&1
    docker rm %IMAGE_NAME% >nul 2>&1
)

echo Running Docker container: %IMAGE_NAME% on port %HOST_PORT%
docker run -d --name %IMAGE_NAME% -p %HOST_PORT%:5432 %IMAGE_NAME%

echo Creating Python virtual environment at %ENV_DIR%
python -m venv %ENV_DIR%

echo Activating virtual environment...
call %ENV_DIR%\Scripts\activate.bat

echo Installing uv...
pip install uv

echo Creating uv virtual environment...
uv venv

echo Syncing dependencies with uv...
uv sync

echo Running Uvicorn app...
uv run --active uvicorn app.main:app --reload

endlocal