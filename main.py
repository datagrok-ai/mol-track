from typing import Any
from fastapi import FastAPI, Depends, HTTPException, Body
from sqlalchemy.orm import Session
from sqlalchemy.sql import text
from rdkit import Chem

import models
from chemistry_utils import (
    calculate_no_stereo_smiles_hash,
    calculate_tautomer_hash,
    standardize_mol,
)
from logging_setup import logger

# Handle both package imports and direct execution
try:
    # When imported as a package (for tests)
    from . import models, crud
    from .database import SessionLocal
except ImportError:
    # When run directly
    import models
    import crud
    from database import SessionLocal

# models.Base.metadata.create_all(bind=engine)

app = FastAPI(title="MolTrack API", description="API for managing chemical compounds and batches")

admin_user_id: str | None = None


# Dependency
def get_db():
    db = SessionLocal()
    logger.debug(db.bind.url)
    logger.debug("Database connection successful")
    try:
        yield db
    finally:
        db.close()


def get_admin_user(db: Session):
    admin = db.query(models.User).filter(models.User.first_name == "Admin").first()
    if not admin:
        raise Exception("Admin user not found.")

    global admin_user_id
    admin_user_id = admin.id


@app.on_event("startup")
def on_startup():
    db = SessionLocal()
    try:
        get_admin_user(db)
    finally:
        db.close()


# Compounds endpoints
@app.post("/compounds/", response_model=models.CompoundResponse)
def create_compound(compound: models.CompoundCreate, db: Session = Depends(get_db)):
    return crud.create_compound(db=db, compound=compound)



@app.post("/compounds/batch/", response_model=list[models.CompoundResponse])
def create_compounds_batch(compounds: list[str] = Body(..., embed=True), db: Session = Depends(get_db)):
    """
    Create multiple compounds from a list of SMILES strings.

    All SMILES must be valid and not already exist in the database.
    If any SMILES is invalid or already exists, the entire batch will fail.
    """
    return crud.create_compounds_batch(db=db, smiles_list=compounds)


# Think of removing the schema at all and use as params
@app.get("/compounds/", response_model=list[models.CompoundResponse])
def read_compounds(query: models.CompoundQueryParams = Depends(), db: Session = Depends(get_db)):
    """
    Get a list of compounds with optional filtering by substructure.

    - **substructure**: Optional SMILES pattern to search for substructures
    - **skip**: Number of records to skip (for pagination)
    - **limit**: Maximum number of records to return (for pagination)
    """
    compounds = crud.get_compounds_ex(db, query_params=query)
    return compounds


@app.get("/compounds/{compound_id}", response_model=models.CompoundResponse)
def read_compound(compound_id: int, db: Session = Depends(get_db)):
    db_compound = crud.get_compound(db, compound_id=compound_id)
    if db_compound is None:
        raise HTTPException(status_code=404, detail="Compound not found")
    return db_compound


# Batches endpoints
@app.post("/batches/", response_model=models.BatchResponse)
def create_batch(batch: models.BatchBase, db: Session = Depends(get_db)):
    return crud.create_batch(db=db, batch=batch)


@app.get("/batches/", response_model=list[models.BatchResponse])
def read_batches(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    batches = crud.get_batches(db, skip=skip, limit=limit)
    return batches


@app.get("/batches/{batch_id}", response_model=models.BatchResponse)
def read_batch(batch_id: int, db: Session = Depends(get_db)):
    db_batch = crud.get_batch(db, batch_id=batch_id)
    if db_batch is None:
        raise HTTPException(status_code=404, detail="Batch not found")
    return db_batch


@app.post("/semantic-types/", response_model=models.SemanticType)
def create_semantic_type_endpoint(semantic_type: models.SemanticTypeBase, db: Session = Depends(get_db)):
    return crud.create_semantic_type(db=db, semantic_type=semantic_type)


# Properties endpoints
@app.post("/properties/", response_model=models.PropertyResponse)
def create_property(property: models.PropertyBase, db: Session = Depends(get_db)):
    return crud.create_property(db=db, property=property)


@app.get("/properties/", response_model=list[models.PropertyResponse])
def read_properties(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    properties = crud.get_properties(db, skip=skip, limit=limit)
    return properties


@app.get("/properties/{property_id}", response_model=models.PropertyResponse)
def read_property(property_id: int, db: Session = Depends(get_db)):
    db_property = crud.get_property(db, property_id=property_id)
    if db_property is None:
        raise HTTPException(status_code=404, detail="Property not found")
    return db_property


# AssayType endpoints
@app.post("/assay-types/", response_model=models.AssayTypeResponse)
def create_assay_type(assay_type: models.AssayTypeCreate, db: Session = Depends(get_db)):
    # Validate that all property IDs exist
    if assay_type.property_ids:
        for property_id in assay_type.property_ids:
            db_property = crud.get_property(db, property_id=property_id)
            if db_property is None:
                raise HTTPException(status_code=404, detail=f"Property with ID {property_id} not found")

    return crud.create_assay_type(db=db, assay_type=assay_type)


@app.get("/assay-types/", response_model=list[models.AssayTypeResponse])
def read_assay_types(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    assay_types = crud.get_assay_types(db, skip=skip, limit=limit)
    return assay_types


@app.get("/assay-types/{assay_type_id}", response_model=models.AssayTypeResponse)
def read_assay_type(assay_type_id: int, db: Session = Depends(get_db)):
    db_assay_type = crud.get_assay_type(db, assay_type_id=assay_type_id)
    if db_assay_type is None:
        raise HTTPException(status_code=404, detail="Assay type not found")
    return db_assay_type


# Assay endpoints
@app.post("/assays/", response_model=models.AssayResponse)
def create_assay(assay: models.AssayCreate, db: Session = Depends(get_db)):
    # Validate that the assay type exists
    db_assay_type = crud.get_assay_type(db, assay_type_id=assay.assay_type_id)
    if db_assay_type is None:
        raise HTTPException(status_code=404, detail=f"Assay type with ID {assay.assay_type_id} not found")

    # Validate that all property IDs exist
    if assay.property_ids:
        for property_id in assay.property_ids:
            db_property = crud.get_property(db, property_id=property_id)
            if db_property is None:
                raise HTTPException(status_code=404, detail=f"Property with ID {property_id} not found")

    return crud.create_assay(db=db, assay=assay)


@app.get("/assays/", response_model=list[models.AssayResponse])
def read_assays(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    assays = crud.get_assays(db, skip=skip, limit=limit)
    return assays


@app.get("/assays/{assay_id}", response_model=models.AssayResponse)
def read_assay(assay_id: int, db: Session = Depends(get_db)):
    db_assay = crud.get_assay(db, assay_id=assay_id)
    if db_assay is None:
        raise HTTPException(status_code=404, detail="Assay not found")
    return db_assay


# AssayResult endpoints
# a bit different response
@app.post("/assay-results/", response_model=models.AssayResultResponse)
def create_assay_result(assay_result: models.AssayResultBase, db: Session = Depends(get_db)):
    """
    Create a single assay result entry for a specific property.

    - **assay_id**: The ID of the assay
    - **batch_id**: The ID of the batch
    - **property_id**: The ID of the property
    - **value_num/value_string/value_bool**: The value of the measurement (use the appropriate field based on property type)
    """
    return crud.create_assay_result(db=db, assay_result=assay_result)


@app.post("/batch-assay-results/", response_model=models.BatchAssayResultsResponse)
def create_batch_assay_results(batch_results: models.BatchAssayResultsCreate, db: Session = Depends(get_db)):
    """
    Register multiple measurements for a batch against an assay at once.

    - **assay_id**: The ID of the assay
    - **batch_id**: The ID of the batch
    - **measurements**: A dictionary mapping property names to their values
    """
    return crud.create_batch_assay_results(db=db, batch_results=batch_results)


@app.get("/assay-results/", response_model=list[models.AssayResultResponse])
def read_assay_results(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    """
    Get a list of all individual assay results.
    """
    return crud.get_assay_results(db, skip=skip, limit=limit)


@app.get("/assay-results/{assay_result_id}", response_model=models.AssayResultResponse)
def read_assay_result(assay_result_id: int, db: Session = Depends(get_db)):
    """
    Get a specific assay result by ID.
    """
    db_assay_result = crud.get_assay_result(db, assay_result_id=assay_result_id)
    if db_assay_result is None:
        raise HTTPException(status_code=404, detail="Assay result not found")
    return db_assay_result


@app.get("/batches/{batch_id}/assay-results", response_model=list[models.BatchAssayResultsResponse])
def read_batch_assay_results(batch_id: int, db: Session = Depends(get_db)):
    """
    Get all assay results for a specific batch, grouped by assay.
    """
    return crud.get_batch_assay_results(db, batch_id=batch_id)


# BatchDetail endpoints
@app.post("/batch-details/", response_model=models.BatchDetailResponse)
def create_batch_detail(batch_detail: models.BatchDetailBase, db: Session = Depends(get_db)):
    """
    Create a batch detail entry.

    - **batch_id**: The ID of the batch
    - **property_id**: The ID of the property
    - **result_value**: The value for this property
    """
    # Validate batch exists
    db_batch = crud.get_batch(db, batch_id=batch_detail.batch_id)
    if db_batch is None:
        raise HTTPException(status_code=404, detail=f"Batch with ID {batch_detail.batch_id} not found")

    # Validate property exists
    db_property = crud.get_property(db, property_id=batch_detail.property_id)
    if db_property is None:
        raise HTTPException(status_code=404, detail=f"Property with ID {batch_detail.property_id} not found")

    return crud.create_batch_detail(db=db, batch_detail=batch_detail)


@app.get("/batch-details/{batch_detail_id}", response_model=models.BatchDetailResponse)
def read_batch_detail(batch_detail_id: int, db: Session = Depends(get_db)):
    """
    Get a specific batch detail by ID.

    - **batch_detail_id**: The ID of the batch detail
    """
    db_batch_detail = crud.get_batch_detail(db, batch_detail_id=batch_detail_id)
    if db_batch_detail is None:
        raise HTTPException(status_code=404, detail="Batch detail not found")
    return db_batch_detail


@app.get("/batches/{batch_id}/details", response_model=list[models.BatchDetailResponse])
def read_batch_details(batch_id: int, skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    """
    Get all details for a specific batch.

    - **batch_id**: The ID of the batch
    - **skip**: Number of records to skip
    - **limit**: Maximum number of records to return
    """
    # Validate batch exists
    db_batch = crud.get_batch(db, batch_id=batch_id)
    if db_batch is None:
        raise HTTPException(status_code=404, detail=f"Batch with ID {batch_id} not found")

    return crud.get_batch_details_by_batch(db, batch_id=batch_id, skip=skip, limit=limit)


@app.post("/search/compounds/exact", response_model=list[models.Compound])
def search_compounds_exact(request: models.ExactSearchModel, db: Session = Depends(get_db)):
    """
    Endpoint for exact compound search.
    """
    try:
        # Validate and generate hash_mol using the Pydantic model
        exact_params = models.ExactSearchModel(**request.model_dump())

        # Use the generated hash_mol to query the database
        return crud.get_compound_by_hash(db=db, hash_mol=exact_params.hash_mol)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


# Less precise search
@app.post("/search/compounds/structure", response_model=list[models.Compound])
def search_compound_structure(request: models.SearchCompoundStructure, db: Session = Depends(get_db)):
    """
    Perform a dynamic structure-based search for compounds.

    - **search_type**: Type of structure search (e.g., "substructure", "tautomer", "stereo", "similarity", "connectivity").
    - **query_smiles**: SMILES string for the structure search.
    - **search_parameters**: Additional parameters for the search.
    """
    try:
        query_smiles = request.query_smiles

        search_functions = {
            "substructure": crud.search_compounds_substructure,
            "tautomer": crud.search_compounds_tautomer,
            "stereo": crud.search_compounds_stereo,
            "connectivity": crud.search_compounds_connectivity,
            "similarity": crud.search_compounds_similarity,
        }

        search_func = search_functions.get(request.search_type)
        if not search_func:
            raise HTTPException(status_code=400, detail="Invalid search type")

        return search_func(db=db, query_smiles=query_smiles, search_parameters=request.search_parameters)

    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.post("/search/complex", response_model=list[dict[str, Any]])
def search_complex_query(request: models.ComplexQueryRequest, db: Session = Depends(get_db)):
    """
    Perform a complex query across multiple tables (e.g., compounds, batch, assays).
    """
    try:
        # Validate and build the query
        query_parts = []
        query_params = {}  # Dictionary to hold query parameters
        selected_columns = {}  # Dictionary to track selected columns for each table

        for idx, condition in enumerate(request.conditions):
            table = condition.table
            field = condition.field
            operator = condition.operator
            value = condition.value
            if table not in ["compounds", "batch", "assays"]:
                raise HTTPException(status_code=400, detail=f"Invalid table name: {table}")

            # Handle SMILES-based queries substructure, tautomer, similarity etc
            if condition.query_smiles:
                mol = Chem.MolFromSmiles(condition.query_smiles)
                if mol is None:
                    raise HTTPException(
                        status_code=400,
                        detail=f"Invalid SMILES string: {condition.query_smiles}",
                    )
                standardized_mol = standardize_mol(mol)

                # Calculate the appropriate hash based on the field
                if field == "hash_tautomer":
                    value = calculate_tautomer_hash(standardized_mol)
                elif field == "hash_no_stereo_smiles":
                    value = calculate_no_stereo_smiles_hash(standardized_mol)
                else:
                    raise HTTPException(status_code=400, detail=f"Unsupported hash field: {field}")

            # If the field is a UUID field, inline the UUID value directly into the query
            if field in ["hash_tautomer", "hash_no_stereo_smiles"]:
                query_parts.append(f"{table}.{field} {operator} '{value}'")
            else:
                # For non-UUID fields, use parameterized queries
                param_name = f"param_{idx}"  # Unique parameter name
                query_parts.append(f"{table}.{field} {operator} :{param_name}")
                query_params[param_name] = value  # Store the parameter value

            if condition.columns:
                selected_columns[table] = condition.columns
            elif table == "compounds":
                # Default columns to return for the compounds table
                selected_columns[table] = ["id", "canonical_smiles"]

        # Combine conditions with the specified logic
        combined_conditions = f" {request.logic} ".join(query_parts)

        select_clauses = []
        for table, columns in selected_columns.items():
            for column in columns:
                select_clauses.append(f"{table}.{column}")
        select_clause = ", ".join(select_clauses)

        # Build the final SQL query
        sql_query = f"""
            SELECT {select_clause}
            FROM {", ".join(set([cond.table for cond in request.conditions]))}
            WHERE {combined_conditions}
        """

        # Execute the query
        result = db.execute(text(sql_query), query_params)
        column_names = [col[0] for col in result.cursor.description]

        # Convert the result into a list of dictionaries
        json_result = [dict(zip(column_names, row)) for row in result]
        return json_result

    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Error building or executing query: {str(e)}")
