import pytest
import os
import uuid
import sys
from fastapi.testclient import TestClient
from sqlalchemy import create_engine, text
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import ProgrammingError

# Set the DB_SCHEMA environment variable
os.environ["DB_SCHEMA"] = "moltrack"

# Add the parent directory to sys.path to allow imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Now import from the project directly
from main import app, get_db
from database import Base, SQLALCHEMY_DATABASE_URL

# Create a unique test database name
test_db_suffix = str(uuid.uuid4())[:8]
test_db_name = f"test_moltrack_{test_db_suffix}"

# Extract connection details from the main database URL
# Assuming format: postgresql://username:password@host:port/dbname
db_url_parts = SQLALCHEMY_DATABASE_URL.split('/')
base_url = '/'.join(db_url_parts[:-1])
admin_db = db_url_parts[-1].split('?')[0]  # Get the database name without query parameters
ADMIN_DATABASE_URL = f"{base_url}/{admin_db}"
TEST_DATABASE_URL = f"{base_url}/{test_db_name}"

# Set the DATABASE_URL environment variable for testing
os.environ["DATABASE_URL"] = TEST_DATABASE_URL

# Create engine for admin connection (to create/drop test database)
admin_engine = create_engine(ADMIN_DATABASE_URL, 
        isolation_level="AUTOCOMMIT")

# Get the schema name from environment variable or use default
DB_SCHEMA = os.environ.get("DB_SCHEMA", "moltrack")

@pytest.fixture(scope="module")
def setup_test_db():
    """
    Set up a test PostgreSQL database using SQL commands through SQLAlchemy
    """
    # Create the test database
    with admin_engine.connect() as conn:
        # Disconnect all active connections to the test database if it exists
        conn.execute(text(f"""
            SELECT pg_terminate_backend(pg_stat_activity.pid)
            FROM pg_stat_activity
            WHERE pg_stat_activity.datname = '{test_db_name}'
            AND pid <> pg_backend_pid();
        """))
        
        # Drop the database if it exists
        try:
            conn.execute(text(f"DROP DATABASE IF EXISTS {test_db_name}"))
        except Exception as e:
            print(f"Error dropping database: {e}")
        
        # Create the test database
        conn.execute(text(f"CREATE DATABASE {test_db_name}"))

        # Grant privileges
        conn.execute(text(f"GRANT ALL PRIVILEGES ON DATABASE {test_db_name} TO CURRENT_USER"))
    
    # Create engine for the test database
    test_engine = create_engine(TEST_DATABASE_URL)
    
    # Read the schema from db.sql
    db_sql_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'db', 'schema.sql')
    with open(db_sql_path, 'r') as f:
        schema_sql = f.read()
    schema_sql = schema_sql.replace(':LOGIN', 'postgres')
    
    # Apply the schema to the test database
    with test_engine.connect() as conn:
        # Execute each statement in the schema
        for statement in schema_sql.split(';'):
            statement = statement.replace('^', ';')   # had to replace it for trigger function 
            if statement.strip():
                conn.execute(text(statement))
        conn.execute(text("COMMIT"))
    
    yield
    
    # Drop the test database after tests
    with admin_engine.connect() as conn:
        # Disconnect all active connections to the test database
        conn.execute(text(f"""
            SELECT pg_terminate_backend(pg_stat_activity.pid)
            FROM pg_stat_activity
            WHERE pg_stat_activity.datname = '{test_db_name}'
            AND pid <> pg_backend_pid();
        """))
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
        db.execute(text(f"TRUNCATE {DB_SCHEMA}.compounds, {DB_SCHEMA}.batches, {DB_SCHEMA}.properties RESTART IDENTITY CASCADE;"))
        db.commit()

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

# Common test data
aspirin_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
aspirin_smiles_noncanonical = "CC(Oc1c(C(O)=O)cccc1)=O"
caffeine_smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" 