from fastapi import FastAPI, Depends, HTTPException, Body
import re
from sqlalchemy.orm import Session
from typing import List
import models

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
    print(db.bind.url)
    print("Database connection successful")
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


@app.post("/compounds/batch/", response_model=List[models.CompoundResponse])
def create_compounds_batch(compounds: List[str] = Body(..., embed=True), db: Session = Depends(get_db)):
    """
    Create multiple compounds from a list of SMILES strings.

    All SMILES must be valid and not already exist in the database.
    If any SMILES is invalid or already exists, the entire batch will fail.
    """
    return crud.create_compounds_batch(db=db, smiles_list=compounds)


# Think of removing the schema at all and use as params
@app.get("/compounds/", response_model=List[models.CompoundResponse])
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


@app.get("/batches/", response_model=List[models.BatchResponse])
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


@app.get("/properties/", response_model=List[models.PropertyResponse])
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


@app.get("/assay-types/", response_model=List[models.AssayTypeResponse])
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


@app.get("/assays/", response_model=List[models.AssayResponse])
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


@app.get("/assay-results/", response_model=List[models.AssayResultResponse])
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


@app.get("/batches/{batch_id}/assay-results", response_model=List[models.BatchAssayResultsResponse])
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


@app.get("/batches/{batch_id}/details", response_model=List[models.BatchDetailResponse])
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

# SynonymType endpoints
@app.post("/synonym-types/", response_model=models.SynonymTypeResponse)
def create_synonym_type(synonym_type: models.SynonymTypeBase, db: Session = Depends(get_db)):
    return crud.create_synonym_type(db=db, synonym_type=synonym_type)

@app.get("/synonym-types/", response_model=List[models.SynonymTypeResponse])
def read_synonym_types(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    return crud.get_synonym_types(db, skip=skip, limit=limit)

@app.get("/synonym-types/{type_id}", response_model=models.SynonymTypeResponse)
def read_synonym_type(type_id: int, db: Session = Depends(get_db)):
    synonym_type = crud.get_synonym_type(db, type_id=type_id)
    if not synonym_type:
        raise HTTPException(status_code=404, detail="Synonym type not found")
    return synonym_type

# Synonym validation function
def validate_synonym_value(db: Session, synonym_type_id: int, synonym_value: str):
    """
    Validate the synonym value against the pattern defined in the synonym type.

    Args:
        db (Session): Database session.
        synonym_type_id (int): ID of the synonym type.
        synonym_value (str): Synonym value to validate.

    Raises:
        HTTPException: If the synonym type is not found or the value does not match the pattern.
    """
    # Retrieve the synonym type and its pattern
    synonym_type = db.query(models.SynonymType).filter(models.SynonymType.id == synonym_type_id).first()
    if not synonym_type:
        raise HTTPException(status_code=404, detail="Synonym type not found")

    # Validate the synonym_value against the pattern
    if synonym_type.pattern and not re.match(synonym_type.pattern, synonym_value):
        raise HTTPException(
            status_code=422,
            detail=f"Synonym value '{synonym_value}' does not match the required pattern: {synonym_type.pattern}"
        )

# Compound synonym endpoints
@app.post("/compound-synonyms/", response_model=models.CompoundSynonymResponse)
def create_compound_synonym(compound_synonym: models.CompoundSynonymBase, db: Session = Depends(get_db)):
    validate_synonym_value(db, compound_synonym.synonym_type_id, compound_synonym.synonym_value)
    return crud.create_compound_synonym(db=db, synonym=compound_synonym)

@app.get("/compound-synonyms/", response_model=List[models.CompoundSynonymResponse])
def read_compound_synonyms(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    return crud.get_compound_synonyms(db, skip=skip, limit=limit)

@app.get("/compound-synonyms/{synonym_id}", response_model=models.CompoundSynonymResponse)
def read_compound_synonym(synonym_id: int, db: Session = Depends(get_db)):
    synonym = crud.get_compound_synonym(db, synonym_id=synonym_id)
    if not synonym:
        raise HTTPException(status_code=404, detail="Compound synonym not found")
    return synonym

# @app.put("/compound-synonyms/{synonym_id}", response_model=models.CompoundSynonymResponse)
# def update_compound_synonym(synonym_id: int, compound_synonym: models.CompoundSynonym, db: Session = Depends(get_db)):
#     validate_synonym_value(db, compound_synonym.synonym_type_id, compound_synonym.synonym_value)
#     return crud.update_compound_synonym(db=db, compound_synonym_id=synonym_id, compound_synonym=compound_synonym)

# Batch synonym endpoints
@app.post("/batch-synonyms/", response_model=models.BatchSynonymResponse)
def create_batch_synonym(batch_synonym: models.BatchSynonymBase, db: Session = Depends(get_db)):
    validate_synonym_value(db, batch_synonym.synonym_type_id, batch_synonym.synonym_value)
    return crud.create_batch_synonym(db=db, synonym=batch_synonym)

@app.get("/batch-synonyms/", response_model=List[models.BatchSynonymResponse])
def read_batch_synonyms(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    return crud.get_batch_synonyms(db, skip=skip, limit=limit)

@app.get("/batch-synonyms/{synonym_id}", response_model=models.BatchSynonymResponse)
def read_batch_synonym(synonym_id: int, db: Session = Depends(get_db)):
    synonym = crud.get_batch_synonym(db, synonym_id=synonym_id)
    if not synonym:
        raise HTTPException(status_code=404, detail="Batch synonym not found")
    return synonym

# @app.put("/batch-synonyms/{synonym_id}", response_model=models.BatchSynonymResponse)
# def update_batch_synonym(synonym_id: int, batch_synonym: models.BatchSynonym, db: Session = Depends(get_db)):
#     validate_synonym_value(db, batch_synonym.synonym_type_id, batch_synonym.synonym_value)
#     return crud.update_batch_synonym(db=db, batch_synonym_id=synonym_id, batch_synonym=batch_synonym)

# Synonym Search
@app.get("/synonym-search/compounds", response_model=List[models.CompoundResponse])
def search_compounds_by_synonym(synonym_value: str, skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    """
    Search for compounds by their synonyms.

    - **synonym_value**: The synonym value to search for
    - **skip**: Number of records to skip (for pagination)
    - **limit**: Maximum number of records to return (for pagination)
    """
    return crud.search_compounds_by_synonym(db, synonym_value, skip, limit)

@app.get("/synonym-search/batches", response_model=List[models.BatchResponse])
def search_batches_by_synonym(synonym_value: str, skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    """
    Search for batches by their synonyms.

    - **synonym_value**: The synonym value to search for
    - **skip**: Number of records to skip (for pagination)
    - **limit**: Maximum number of records to return (for pagination)
    """
    return crud.search_batches_by_synonym(db, synonym_value, skip, limit)

# Synonym retrieval by compound or batch ID
@app.get("/compounds/{compound_id}/synonyms", response_model=List[models.CompoundSynonymResponse])
def read_compound_synonyms_by_compound_id(compound_id: int, db: Session = Depends(get_db)):
    """
    Get all synonyms for a specific compound.

    - **compound_id**: The ID of the compound
    """
    # Query the database for synonyms related to the compound
    synonyms = db.query(models.CompoundSynonym).filter(models.CompoundSynonym.compound_id == compound_id).all()
    if not synonyms:
        raise HTTPException(status_code=404, detail="No synonyms found for the specified compound")
    return synonyms

@app.get("/batches/{batch_id}/synonyms", response_model=List[models.BatchSynonymResponse])
def read_batch_synonyms_by_batch_id(batch_id: int, db: Session = Depends(get_db)):
    """
    Get all synonyms for a specific batch.

    - **batch_id**: The ID of the batch
    """
    # Query the database for synonyms related to the batch
    synonyms = db.query(models.BatchSynonym).filter(models.BatchSynonym.batch_id == batch_id).all()
    if not synonyms:
        raise HTTPException(status_code=404, detail="No synonyms found for the specified batch")
    return synonyms
