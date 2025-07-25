from pathlib import Path
from functools import wraps
import time
import json
import logging
import httpx
from httpx import HTTPError, RequestError
from data_generator import DataGenerator

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
    def __init__(self, size: str, base_url="http://localhost:8000"):
        self.scopes = {
            "compound": "/v1/compounds/",
            "batch": "/v1/batches/",
            "assay": "/v1/assays",
            "assay_run": "/v1/assay_runs/",
            "assay_result": "/v1/assay_results/",
        }
        self.times = {}
        self.base_url = base_url
        self.size = size
        self.client = httpx.Client(base_url=base_url, timeout=100_000)

    def __enter__(self):
        return self

    def __exit__(self, _exc_type, _exc_value, _traceback):
        self.client.close()

    def run_test_suite(self):
        _, gen_time = self.generate_data()
        _, pop_time = self.populate_database()

        scope_performance = {
            scope: f"Total time for adding {scope}s to database: "
            f"{self.times[scope]:.2f} seconds for a total of {self.total_records[scope]} records"
            for scope in self.times
        }

        return {
            "generate_data_time": f"{gen_time:.2f} seconds",
            "populate_database_time": f"{pop_time:.2f} seconds",
            "time_per_scope": scope_performance,
            "total": f"Total time for all operations: {(gen_time + pop_time):.2f} seconds",
        }

    @time_it
    def generate_data(self):
        dg = DataGenerator(size=self.size, output_dir=str(DATA_DIR))
        self.total_records = dg.generate_all()

    @time_it
    def populate_database(self):
        try:
            # Add schemas
            for scope in self.scopes:
                self._post_json_file(file_path=DATA_DIR / f"schemas/{scope}_schema.json", url="/v1/schema/")

            # Add assays
            self._post_json_file(file_path=DATA_DIR / "assays.json", url="/v1/assays")

            # Add compounds, batches, assay runs, and assay results
            for scope, url in self.scopes.items():
                if scope != "assay":
                    _, measured_time = self._post_csv_with_mapping(
                        file_path=DATA_DIR / f"{scope}.csv", url=url, mapping_file=DATA_DIR / f"{scope}_mapping.json"
                    )
                    self.times[scope] = measured_time
        except PerformanceTesterError as e:
            logger.error(f"Failed to add data for scope '{scope}': {e}")
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
    with PerformanceTester("xsmall", base_url="http://localhost:8000") as tp:
        results = tp.run_test_suite()
        logger.info(json.dumps(results, indent=4))
