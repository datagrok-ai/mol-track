# MolTrack Server

A lightweight, flexible and extendable FastAPI-based server for managing chemical compounds, batches, and properties, with the RDKit cartridge-enabled Postgres for chemical intelligence. Ideal for labs, startups, and small- to medium-sized biotech companies.

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

1. Create a virtual environment:
   ```
   python -m venv venv
   ```

2. Activate the virtual environment:
   - Windows: `venv\Scripts\activate`
   - Linux/Mac: `source venv/bin/activate`

3. Install dependencies:
   ```
   pip install -r requirements.txt
   ```

4. Configure the database connection:
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

## User Stories
Here we build a list of user stories that we should support as we progress.

- As an administrator, I can add a chemist to the user management system, so they may register their batches.
- As a chemist, I want to add a batch of a compound that I have synthesized and receive an appropriate batch identifier and compound identifier so that the biologist can generate identifiable unique assay results.
- As a chemist, I want to edit a batch that I have registered to correct or add additional information or update the structure.
- As a chemist, I should only be allowed to update the compound structure if I own all the associated batches, so that I do not perturb the registration number for batches/compounds made by others.
- As a chemist, the molecule coordinates that I use to register should be saved to preserve an appropriate orientation.
- 
- As a registrar, I want to add batches (and compounds) in bulk to support the case when we might acquire a set of batches and compounds for research purposes.
- As a registrar, I want to edit the batches (and compounds) that others have registered to correct errors.
- As a biologist, I want to view all the information for a specific batch so that I might prepare my experiment correctly.
- As an administrator, I want to configure the business rules for batch uniquenss verification.
- As an administrator, I want to configure the business rules for compound uniqueness verification.
- As an administrator, I want to configure the business rules for no-structure uniqueness.
- As a chemist, I want to register a batch with a no-structure compound and receive appropriate batch and compound identifiers.
- As a chemist or registrar, I want to annotate my registration with aliases, so that I can capture external identifiers.
- As chemist, I want to query compounds using substructure, similarity, exact, no-stereo exact or tautomeric methods, so that I might retrieve the exact subset that is scientifically interesting.
- As a scientist, I want to query batches by properties and/or aliases to find the relevant batches.
- As a scientist, I want to query compounds by properties and/or aliases to find the relevant compounds.
- As a peptide chemist, I want to register my compounds using a helm representation, 
- As an administrator, I want to configure my system so that certain properties must adhere to controlled vocabularies, in order to maintain data cleanliness.
- As the system, I should be able to provide a preferred rendering of the compound as an image or svg, so that I can used in downstream systems that don't have native molecular rendering.
- As a chemist or registrar, I should be able to associate a preferred rendering for compound, so that consumers will be able to retrieve and benefit from the preferred rendering.


## Persona / Roles
- Administrator has all the appropriate rights to administer the system
- Chemist (or editor) has the ability to add batches and compounds and update those that created.
- Registar has the ability to add batches and compounds in bulk and edit all compounds and batches.
- Viewer has the ability to search and retrieve compounds, batches/lots and appropriate properties

## Concepts
- As a compound registration system, meaning registering and tracking the parent structure and obtaining a unique corporate id for the structure.  
- As a batch registration system, where the batch entity is the primary entity and compound represents a structural annotation for the batch.
- Salts and/or solvents are an additional annotation on the batch, but all serve as a point of aggregation - 'find all salt/solvent form for batches of this structure.'
- The inital focus is small molecules and no-structs, but the system should be evolvable to handle other modalities: peptides, siRNA, ASO, biologics.