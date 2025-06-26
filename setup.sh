#!/bin/bash
set -e

IMAGE_NAME="moltrack"
ENV_DIR=".moltrack-env"

echo "Building Docker image: $IMAGE_NAME from Dockerfile"
docker build -t $IMAGE_NAME .

function find_free_port() {
    for port in $(seq 5432 5500); do
        if ! lsof -i:$port >/dev/null 2>&1; then
            echo $port
            return
        fi
    done
    echo "No free port found" >&2
    exit 1
}

HOST_PORT=$(find_free_port)
echo "Found free host port: $HOST_PORT"

export DB_PORT=$HOST_PORT
echo "Set DB_PORT environment variable to $DB_PORT"

if [ "$(docker ps -aq -f name=$IMAGE_NAME)" ]; then
    echo "Stopping and removing existing container: $IMAGE_NAME"
    docker stop $IMAGE_NAME || true
    docker rm $IMAGE_NAME
fi

echo "Running Docker container: $IMAGE_NAME"
docker run -d --name $IMAGE_NAME -p $HOST_PORT:5432 $IMAGE_NAME

echo "Creating Python virtual environment at $ENV_DIR"
python3 -m venv $ENV_DIR

echo "Activating virtual environment..."
source $ENV_DIR/bin/activate

if ! command -v uv &> /dev/null; then
    echo "Installing uv..."
    pip install uv
fi

echo "Creating uv virtual environment..."
uv venv || true

echo "Syncing dependencies with uv..."
uv sync

echo "Running Uvicorn app..."
uv run --active uvicorn app.main:app --reload