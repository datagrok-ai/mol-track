# MolTrack Server

A lightweight, flexible and extendable FastAPI-based server for managing chemical compounds, batches, and properties, with the RDKit cartridge-enabled Postgres for chemical intelligence. Ideal for labs, startups, and small- to medium-sized biotech companies.

See also [user stories](./docs/user-stories.md)

## Features

Work in progress:

1. Compound Registration
    * [x] **Unique Identifiers**: Automatically assigns unique identifiers (e.g., registration numbers, UUIDs) to new compounds.
    * [x] **Duplicate Detection**: Prevents registration of duplicates.
    * [x] **Structure Validation**: Valence checking, standardization of stereochemistry, etc.
    * [ ] **Structure Standardization**: Converts entered structures into a consistent format, handling tautomerization, salts, and stereo conventions.
2. Metadata
    * [ ] **Custom Attributes**: Supports capturing custom metadata (e.g., biological data, physicochemical properties, origin information) and ties it to the appropriate entity (compound, batch/lot).
    * [ ] **Attachment Management**: Allows attaching documents (NMR spectra, mass spectrometry data, analytical certificates).
3. Batches and Lots
    * [ ] **Batch Registration**: Manages registration of multiple batches or lots for a single compound.
    * [ ] **Duplicate Detection**: Prevents the registration of duplicates
    * [ ] **Purity and Inventory Tracking**: Tracks batch-specific details such as purity, quantity, storage location, supplier, and expiration dates.
4. Protocols and Assay Results
    * [ ] **Protocols**: Define assay types used to measure batches.
    * [ ] **Assay Results**: Register and query assay results.
5. Search
    * [ ] **Structure-based Search**: Supports exact, substructure, similarity, and Markush searches.
    * [ ] **Metadata Search**: Enables querying by metadata fields such as IDs, names, properties, and batch information.
6. Audit and Compliance
    * [ ] **Audit Trails**: Records detailed logs of registration, editing, and deletion activities for compliance and traceability.
    * [ ] **Role-based Access Control**: Implements security controls to ensure sensitive data is accessible only by authorized users.
7. Integration and APIs
    * [x] **API Access**: Provides RESTful APIs to facilitate integration with other lab informatics systems (ELNs, LIMS, inventory management systems).
9. User Interface
    * [ ] **Chemical Drawing Integration**: Allows users to input structures directly using chemical drawing tools (e.g., MarvinJS, ChemDraw, Ketcher).
    * [ ] **Custom Reports**: Generates reports on compound libraries, registration statistics, and inventory statuses.
    * [ ] **Visualization Tools**: Includes dashboards and data visualization features for quick analysis and decision-making.

## Automated setup

To simplify and speed up the installation and launch process, we provide two automated setup scripts:

* `setup.bat` for **Windows**
* `setup.sh` for **macOS/Linux**

[Manual setup](#manual-setup) can be time-consuming and error-prone, requiring multiple steps such as building Docker images, running containers, configuring virtual environments, and starting the server. These scripts handle all of that automatically, so you can get your environment ready with a single command.

Both scripts accept an optional `--run_server` flag (if specified, the Uvicorn server is started).

> **Note:** Docker must be installed and running on your machine before running these scripts.

### On Windows

1. Open CMD in the project directory.
2. Run:

```cmd
setup.bat # Run setup only, skip starting server
setup.bat --run_server # Run setup and start server
```

### On macOS/Linux

1. Open a terminal in the project directory.
2. Make the script executable and run:

```bash
chmod +x setup.sh
./setup.sh # Run setup only, skip starting server
./setup.sh --run_server # Run setup and start server
```

## Manual setup

### 1. Create and activate a virtual environment

Create a new virtual environment:

```bash
python3 -m venv .venv
```

Activate the environment:

* **Windows (CMD):**

  ```cmd
  .venv\Scripts\activate
  ```
* **macOS/Linux:**

  ```bash
  source .venv/bin/activate
  ```

### 2. Install `uv`

Install `uv` package using pip:

```bash
pip install uv
```

*For alternative installation options, see the [official docs](https://docs.astral.sh/uv/guides/install-python/#getting-started).*

### 3. Initialize the project environment with `uv`

Set up the virtual environment and dependencies using `uv` commands:

```bash
uv venv
uv sync
```

* `uv venv` creates the virtual environment.
* `uv sync` installs all required dependencies.

### 4. Configure the database connection

Edit the `database.py` file and update the `SQLALCHEMY_DATABASE_URL` variable with your PostgreSQL connection string:

```python
SQLALCHEMY_DATABASE_URL = "postgresql://user:password@host:port/database"
```

Make sure your database server is running and accessible.

### 5. Run the server

Start the FastAPI server with:

```bash
uv run --active uvicorn app.main:app --reload
```

You can now access the API at [http://localhost:8000](http://localhost:8000).

## Setting up pytest in VS Code

To configure pytest in VS Code, follow these steps:

1. Install the **Python** extension

   * Open the **Extensions** view (`Ctrl+Shift+X` on Windows/Linux or `Cmd+Shift+X` on macOS).
   * Search for **Python** and install the official extension by Microsoft.

2. Click the **Testing** icon (beaker icon) in the **Activity bar**.

3. Configure python tests

   * Click on **Configure Python Tests** button.
   * When prompted, select:

     * **Test framework**: `pytest`
     * **Test directory**: folder containing the tests (important: ensure it contains an `__init__.py` file — this is required for test discovery to work properly)

Your tests should now be detected and listed in the **Testing panel**.


## API Documentation

Once the server is running, you can access:
- Interactive API documentation: http://localhost:8000/docs
- Alternative API documentation: http://localhost:8000/redoc


## API Endpoints

### Compounds
- `GET /compounds/` - List all compounds
- `POST /compounds/` - Create a new compound
- `GET /compounds/{compound_id}` - Get a specific compound
- `PUT /compounds/{compound_id}` - Update a compound
- `DELETE /compounds/{compound_id}` - Delete a compound

### Batches
- `GET /batches/` - List all batches
- `POST /batches/` - Create a new batch
- `GET /batches/{batch_id}` - Get a specific batch
- `PUT /batches/{batch_id}` - Update a batch
- `DELETE /batches/{batch_id}` - Delete a batch

### Properties
- `POST /properties/` - Create a new property
- `GET /compounds/{compound_id}/properties/` - Get properties for a compound 

