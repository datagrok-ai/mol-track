import json
from pathlib import Path
from typing import Any, Dict


def load_config(path: Path) -> Dict[str, Any]:
    """Load configuration settings for the performance tests."""

    with open(path / "config.json") as f:
        config = json.load(f)
        return config


def format_time(seconds):
    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    seconds = seconds % 60
    return f"{int(hours)}:{int(minutes):02}:{int(seconds):02}"
