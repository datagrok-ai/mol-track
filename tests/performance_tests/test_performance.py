from pathlib import Path
from functools import wraps
import time
import json
import logging
import httpx
from httpx import HTTPError, RequestError
from tests.performance_tests.data_generator import DataGenerator
from tests.performance_tests.enums import DatasetSize
from app.utils.enums import ScopeClass
from typing import Dict, Any
from tests.performance_tests.utils import load_config, format_time

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

DATA_DIR = Path("./tests/performance_tests/data")


def time_it(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        elapsed = time.time() - start
        return result, elapsed

    return wrapper


class PerformanceTesterError(Exception):
    """Custom exception for test performance errors."""

    pass


class PerformanceTester:
    def __init__(self, config: Dict[str, Any], base_url: str = "http://localhost:8000"):
        self.scopes = {
            ScopeClass.COMPOUND: "/v1/compounds/",
            ScopeClass.BATCH: "/v1/batches/",
            ScopeClass.ASSAY: "/v1/assays",
            ScopeClass.ASSAY_RUN: "/v1/assay_runs/",
            ScopeClass.ASSAY_RESULT: "/v1/assay_results/",
        }
        self.times = {}
        self.base_url = base_url
        size = DatasetSize(config.get("dataset_size", "SMALL"))
        logger.info(f"Test dataset size: {size.value}")
        self.client = httpx.Client(base_url=base_url, timeout=100_000)
        self.config = config

    def __enter__(self):
        return self

    def __exit__(self, _exc_type, _exc_value, _traceback):
        self.client.close()

    def run_test_suite(self):
        logger.info("Starting performance test suite...")
        _, gen_time = self.generate_data()
        logger.info(f"Data generation completed in {format_time(gen_time)}")
        _, pop_time = self.populate_database()

        scope_performance = {
            scope: f"Total time for adding {scope}s to database: "
            f"{format_time(self.times[scope])} for a total of {self.total_records[scope]} records"
            for scope in self.times
        }

        return {
            "generate_data_time": f"{format_time(gen_time)}",
            "populate_database_time": f"{format_time(pop_time)}",
            "time_per_scope": scope_performance,
            "total": f"Total time for all operations: {format_time(gen_time + pop_time)}",
        }

    @time_it
    def generate_data(self):
        dg = DataGenerator(self.config, output_dir=str(DATA_DIR))
        self.total_records = dg.generate_all()

    @time_it
    def populate_database(self):
        try:
            # Add schemas
            for scope in self.scopes.keys():
                self._post_json_file(file_path=DATA_DIR / f"schemas/{scope.value}_schema.json", url="/v1/schema/")

            # Add assays
            self._post_json_file(file_path=DATA_DIR / "ASSAY.json", url="/v1/assays")

            # Add compounds, batches, assay runs, and assay results
            for scope, url in self.scopes.items():
                if scope != ScopeClass.ASSAY:
                    _, measured_time = self._post_csv_with_mapping(
                        file_path=DATA_DIR / f"{scope.value}.csv",
                        url=url,
                        mapping_file=DATA_DIR / f"{scope.value}_mapping.json",
                    )
                    self.times[scope] = measured_time
                    logger.info(
                        f"Populated database for scope {scope.value} "
                        f"({self.total_records[scope]} records) in {format_time(measured_time)}"
                    )
        except PerformanceTesterError as e:
            logger.error(f"Failed to add data for scope '{scope.value}': {e}")
            raise

    def _read_file(self, path: Path, binary=False):
        mode = "rb" if binary else "r"
        with path.open(mode, encoding=None if binary else "utf-8") as f:
            return f.read() if binary else json.load(f)

    @time_it
    def _post_json_file(self, file_path: Path, url: str):
        data = self._read_file(file_path)
        return self._safe_post(url, json=data)

    @time_it
    def _post_csv_with_mapping(self, file_path: Path, mapping_file: Path, url: str):
        files = {"csv_file": (file_path.name, file_path.open("rb"), "text/csv")}
        data = {
            "error_handling": "reject_row",
            "mapping": json.dumps(self._read_file(mapping_file)),
            "output_format": "json",
        }
        return self._safe_post(url, files=files, data=data)

    def _safe_post(self, url, **kwargs):
        headers = {"accept": "application/json"}
        try:
            response = self.client.post(url, headers=headers, **kwargs)
            response.raise_for_status()
            return response
        except (HTTPError, RequestError) as err:
            raise PerformanceTesterError(f"HTTP error posting to {url}: {err}")
        except Exception as exc:
            raise PerformanceTesterError(f"Unexpected error posting to {url}: {exc}")


if __name__ == "__main__":
    config = load_config(DATA_DIR)
    with PerformanceTester(config, base_url="http://localhost:8000") as tp:
        results = tp.run_test_suite()
        logger.info(json.dumps(results, indent=4))
