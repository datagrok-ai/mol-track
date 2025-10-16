import json
from pathlib import Path
from typing import Optional
import pytest
import os
import uuid
import sys
from fastapi.testclient import TestClient
from sqlalchemy import create_engine, text
from sqlalchemy.orm import sessionmaker
from app.utils.logging_utils import logger
from app.utils import enums

# Set the DB_SCHEMA environment variable
os.environ["DB_SCHEMA"] = "moltrack"

# Add the parent directory to sys.path to allow imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# Now import from the project directly
from app.main import app, get_db
from app.setup.database import SQLALCHEMY_DATABASE_URL

# Create a unique test database name
test_db_suffix = str(uuid.uuid4())[:8]
test_db_name = f"test_moltrack_{test_db_suffix}"

# Extract connection details from the main database URL
# Assuming format: postgresql://username:password@host:port/dbname
db_url_parts = SQLALCHEMY_DATABASE_URL.split("/")
base_url = "/".join(db_url_parts[:-1])
admin_db = db_url_parts[-1].split("?")[0]  # Get the database name without query parameters
ADMIN_DATABASE_URL = f"{base_url}/{admin_db}"
TEST_DATABASE_URL = f"{base_url}/{test_db_name}"

# Set the DATABASE_URL environment variable for testing
os.environ["DATABASE_URL"] = TEST_DATABASE_URL

# Create engine for admin connection (to create/drop test database)
admin_engine = create_engine(ADMIN_DATABASE_URL, isolation_level="AUTOCOMMIT")

# Get the schema name from environment variable or use default
DB_SCHEMA = os.environ.get("DB_SCHEMA", "moltrack")

DATA_DIR = Path(__file__).parent.parent / "data"
BLACK_DIR = DATA_DIR / "black"
SIMPLE_DIR = DATA_DIR / "simple"
EXCLUDE_TABLES = ["users", "settings", "semantic_types", "properties"]
KEEP_PROPERTIES_KEYS = ["corporate_compound_id", "corporate_batch_id"]
SCHEMA_FILES = [
    "batches_schema.json",
    "compounds_schema.json",
    "assays_schema.json",
    "assay_runs_schema.json",
    "assay_results_schema.json",
]
DATA_PATHS = {
    "compounds": ("compounds.csv", "compounds_mapping.json"),
    "batches": ("batches.csv", "batches_mapping.json"),
    "assay_runs": ("assay_runs.csv", "assay_runs_mapping.json"),
    "assay_results": ("assay_results.csv", "assay_results_mapping.json"),
}


def truncate_all_except(db, schema, exclude_tables, keep_properties_keys=None):
    tables_query = text("""
        SELECT tablename FROM pg_tables
        WHERE schemaname = :schema
          AND tablename <> ALL(:exclude_tables)
    """)
    result = db.execute(tables_query, {"schema": schema, "exclude_tables": exclude_tables})
    tables = [row[0] for row in result]

    if tables:
        db.execute(text(f"TRUNCATE {', '.join(f'{schema}.{t}' for t in tables)} RESTART IDENTITY CASCADE;"))

    if "properties" in exclude_tables and keep_properties_keys:
        db.execute(
            text(f"""
                DELETE FROM {schema}.properties
                WHERE name <> ALL(:keep_keys)
            """),
            {"keep_keys": keep_properties_keys},
        )

    db.commit()


@pytest.fixture(scope="module")
def setup_test_db():
    """
    Set up a test PostgreSQL database using SQL commands through SQLAlchemy
    """
    # Create the test database
    with admin_engine.connect() as conn:
        # Disconnect all active connections to the test database if it exists
        conn.execute(
            text(f"""
            SELECT pg_terminate_backend(pg_stat_activity.pid)
            FROM pg_stat_activity
            WHERE pg_stat_activity.datname = '{test_db_name}'
            AND pid <> pg_backend_pid();
        """)
        )

        # Drop the database if it exists
        try:
            conn.execute(text(f"DROP DATABASE IF EXISTS {test_db_name}"))
        except Exception as e:
            logger.error(f"Error dropping database: {e}")

        # Create the test database
        conn.execute(text(f"CREATE DATABASE {test_db_name}"))

        # Grant privileges
        conn.execute(text(f"GRANT ALL PRIVILEGES ON DATABASE {test_db_name} TO CURRENT_USER"))

    # Create engine for the test database
    test_engine = create_engine(TEST_DATABASE_URL)

    # Read the schema from db folder
    schema_paths = [
        os.path.join(os.path.dirname(os.path.dirname(__file__)), "db", "schema.sql"),
        os.path.join(os.path.dirname(os.path.dirname(__file__)), "db", "schema_rdkit.sql"),  # new file
        os.path.join(os.path.dirname(os.path.dirname(__file__)), "db", "setup.sql"),
    ]

    # Apply the schemas to the test database
    with test_engine.connect() as conn:
        for schema_path in schema_paths:
            with open(schema_path, "r") as f:
                schema_sql = f.read()
            schema_sql = schema_sql.replace(":LOGIN", "postgres")

            # Execute each statement in the schema
            for statement in schema_sql.split(";"):
                statement = statement.replace("^", ";")  # had to replace it for trigger function
                statement = statement.replace(":", "\\:")  # escape colons for SQLAlchemy
                if statement.strip():
                    conn.execute(text(statement))
            conn.execute(text("COMMIT"))

    yield

    # Drop the test database after tests
    with admin_engine.connect() as conn:
        # Disconnect all active connections to the test database
        conn.execute(
            text(f"""
            SELECT pg_terminate_backend(pg_stat_activity.pid)
            FROM pg_stat_activity
            WHERE pg_stat_activity.datname = '{test_db_name}'
            AND pid <> pg_backend_pid();
        """)
        )
        conn.execute(text("COMMIT"))

        # Drop the test database
        conn.execute(text(f"DROP DATABASE IF EXISTS {test_db_name}"))
        conn.execute(text("COMMIT"))


# Create engine and session factory for the test database
@pytest.fixture
def test_engine(setup_test_db):
    """Create an engine for the test database"""
    engine = create_engine(TEST_DATABASE_URL)
    yield engine
    engine.dispose()


@pytest.fixture
def test_db(test_engine):
    """Create a fresh database session for each test"""
    TestingSessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=test_engine)
    db = TestingSessionLocal()
    try:
        yield db
    finally:
        # Clean up data after each test
        truncate_all_except(db, DB_SCHEMA, EXCLUDE_TABLES, KEEP_PROPERTIES_KEYS)


@pytest.fixture
def client(test_db):
    """Create a test client with the test database"""

    def override_get_db():
        try:
            yield test_db
        finally:
            pass

    app.dependency_overrides[get_db] = override_get_db

    with TestClient(app) as c:
        yield c

    app.dependency_overrides.clear()


def read_json(file_path: Path) -> dict:
    with file_path.open("r", encoding="utf-8") as f:
        return json.load(f)


def read_csv(file_path: Path) -> str:
    return file_path.read_text(encoding="utf-8")


@pytest.fixture
def preload_schema(client):
    for schema_file in SCHEMA_FILES:
        data = read_json(BLACK_DIR / schema_file)
        client.post("/v1/schema/", json=data)


@pytest.fixture
def preload_additions(client):
    file_path = DATA_DIR / "additions.csv"
    files = {"csv_file": (str(file_path), read_csv(file_path), "text/csv")}
    client.post("/v1/additions/", files=files)


def preload_entity(
    client,
    endpoint: str,
    file_path: str,
    mapping_path: Optional[Path] = None,
    error_handling=enums.ErrorHandlingOptions.reject_row,
    mime_type: str = "text/csv",
    output_format: str = enums.OutputFormat.json,
):
    files = {"file": (str(file_path), open(file_path, "rb"), mime_type)}
    data = {"error_handling": error_handling.value, "output_format": output_format.value}

    if mapping_path is not None:
        data["mapping"] = json.dumps(read_json(mapping_path))

    return client.post(endpoint, files=files, data=data)


def _preload_compounds(
    client, csv_path: str, mapping_path: Optional[str] = None, error_handling=enums.ErrorHandlingOptions.reject_row
):
    return preload_entity(client, "/v1/compounds/", csv_path, mapping_path, error_handling)


@pytest.fixture
def preload_compounds(client):
    return _preload_compounds(client, BLACK_DIR / "compounds.csv", BLACK_DIR / "compounds_mapping.json")


def _preload_batches(
    client, csv_path: str, mapping_path: Optional[str] = None, error_handling=enums.ErrorHandlingOptions.reject_row
):
    return preload_entity(client, "/v1/batches/", csv_path, mapping_path, error_handling)


@pytest.fixture
def preload_batches(client):
    return _preload_batches(client, BLACK_DIR / "batches.csv", BLACK_DIR / "batches_mapping.json")


def _preload_assays(client, json_path: str):
    data = read_json(json_path)
    return client.post("/v1/assays", json=data)


@pytest.fixture
def preload_assays(client):
    return _preload_assays(client, BLACK_DIR / "assays.json")


def _preload_assay_runs(
    client, csv_path: str, mapping_path: Optional[str] = None, error_handling=enums.ErrorHandlingOptions.reject_row
):
    return preload_entity(client, "/v1/assay_runs/", csv_path, mapping_path, error_handling)


@pytest.fixture
def preload_assay_runs(client):
    return _preload_assay_runs(client, BLACK_DIR / "assay_runs.csv", BLACK_DIR / "assay_runs_mapping.json")


def _preload_assay_results(
    client, csv_path: str, mapping_path: Optional[str] = None, error_handling=enums.ErrorHandlingOptions.reject_row
):
    return preload_entity(client, "/v1/assay_results/", csv_path, mapping_path, error_handling)


@pytest.fixture
def preload_assay_results(client):
    return _preload_assay_results(client, BLACK_DIR / "assay_results.csv", BLACK_DIR / "assay_results_mapping.json")


def preload_data_from_dir(client, data_dir: str, use_mapping_files: bool = True):
    for schema_file in SCHEMA_FILES:
        schema_data = read_json(data_dir / schema_file)
        client.post("/v1/schema/", json=schema_data)

    file_path = DATA_DIR / "additions.csv"
    files = {"file": (str(file_path), read_csv(file_path), "text/csv")}
    client.post("/v1/additions/", files=files)

    assays_data = read_json(data_dir / "assays.json")
    client.post("/v1/assays", json=assays_data)

    for endpoint, (csv_file, mapping_file) in {
        "/v1/compounds/": DATA_PATHS["compounds"],
        "/v1/batches/": DATA_PATHS["batches"],
        "/v1/assay_runs/": DATA_PATHS["assay_runs"],
        "/v1/assay_results/": DATA_PATHS["assay_results"],
    }.items():
        csv_path = data_dir / csv_file
        files = {"file": (str(csv_path), read_csv(csv_path), "text/csv")}
        data = {
            "error_handling": enums.ErrorHandlingOptions.reject_row.value,
        }
        if use_mapping_files:
            mapping_path = data_dir / mapping_file
            data["mapping"] = json.dumps(read_json(mapping_path))

        client.post(endpoint, files=files, data=data)


@pytest.fixture
def preload_black_data(client):
    preload_data_from_dir(client, BLACK_DIR, use_mapping_files=True)


@pytest.fixture
def preload_simple_data(client):
    preload_data_from_dir(client, SIMPLE_DIR, use_mapping_files=False)


# Common test data
aspirin_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
aspirin_smiles_noncanonical = "CC(Oc1c(C(O)=O)cccc1)=O"
caffeine_smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
