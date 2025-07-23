import pandas as pd
from pathlib import Path
import json
from typing import Any
import logging

logger = logging.getLogger(__name__)


class FileWriter:
    def __init__(self, output_dir: Path, csv_headers: list[str] = None):
        self.output_dir = output_dir
        self.csv_headers = csv_headers
        self._is_first_csv_chunk = True

    def write_csv_chunk(self, filename: str, rows: list[dict]) -> None:
        file_path = self.output_dir / filename
        df = pd.DataFrame(rows, columns=self.csv_headers)
        mode = "w" if self._is_first_csv_chunk else "a"
        try:
            df.to_csv(file_path, mode=mode, header=self._is_first_csv_chunk, index=False)
            self._is_first_csv_chunk = False
        except Exception as e:
            logger.error(f"Failed to write CSV chunk to {file_path}: {e}")
            raise

    def write_json(self, filename: str, data: dict[str, Any]) -> None:
        file_path = self.output_dir / filename
        try:
            with open(file_path, "w", encoding="utf-8") as f:
                json.dump(data, f, ensure_ascii=False, indent=4)
        except Exception as e:
            logger.error(f"Failed to write JSON to {file_path}: {e}")
            raise
