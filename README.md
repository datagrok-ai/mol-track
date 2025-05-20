# MolTrack Server

A lightweight, flexible and extendable FastAPI-based server for managing chemical compounds, batches, and properties, with the RDKit cartridge-enabled Postgres for chemical intelligence. Ideal for labs, startups, and small- to medium-sized biotech companies.

See also [user stories](user-stories.md)

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



## Setup

1. Install `uv`:
   ```
   pip install uv
   ```
   *For other installation methods, refer to the [official docs](https://docs.astral.sh/uv/guides/install-python/#getting-started).*

2. Create a virtual environment:
   ```
   uv venv
   ```

3. Activate the virtual environment:
   * **Windows**: `.venv\Scripts\activate`
   * **macOS/Linux**: `source .venv/bin/activate`

4. Install dependencies:
   ```
   uv sync
   ```

5. Configure the database connection:
   - Open `database.py` and update the `SQLALCHEMY_DATABASE_URL` with your PostgreSQL connection details.

## Running the Server

Start the server with:

```
uvicorn main:app --reload
```

The API will be available at http://localhost:8000

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

