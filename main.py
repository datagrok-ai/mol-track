import csv
import io
import json
from fastapi import FastAPI, Depends, File, Form, HTTPException, Body, UploadFile
from sqlalchemy.orm import Session
from typing import List, Optional
from batch_registrar import BatchRegistrar
from compound_registrar import CompoundRegistrar
import models
import enums


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


# TODO: Move to the utils or a separate module
def create_if_not_exists(model_cls, base_cls, values, create_fn, created_list, db):
    for item in values:
        if not isinstance(item, dict):
            item = item.dict()

        # TODO: Think how to better handle empty values
        # TODO: Think abouut how to be with column names that are not the same as model fields
        item = {k: v for k, v in item.items() if v not in (None, "", [], {}, ())}
        # TODO: Improve lookup for existing records (if exists we will need an update method, TBD)
        if not db.query(model_cls).filter_by(**item).first():
            try:
                instance = create_fn(db, base_cls(**item))
                created_list.append({"id": instance.id, "name": instance.name})
            except Exception as e:
                db.rollback()
                raise HTTPException(status_code=500, detail=f"Error creating {model_cls.__name__}: {str(e)}")


@app.post("/schema/")
def preload_schema(payload: models.SchemaPayload, db: Session = Depends(get_db)):
    created_synonyms = []
    created_properties = []

    create_if_not_exists(
        model_cls=models.SynonymType,
        base_cls=models.SynonymTypeBase,
        values=payload.synonym_types,
        create_fn=crud.create_synonym_type,
        created_list=created_synonyms,
        db=db,
    )

    create_if_not_exists(
        model_cls=models.Property,
        base_cls=models.PropertyBase,
        values=payload.properties,
        create_fn=crud.create_property,
        created_list=created_properties,
        db=db,
    )

    return {
        "status": "success",
        "created": {
            "synonym_types": created_synonyms,
            "property_types": created_properties,
        },
    }


@app.post("/v1/compounds/")
def register_compounds_v1(
    csv_file: UploadFile = File(...),
    mapping: Optional[str] = Form(None),
    error_handling: enums.ErrorHandlingOptions = Form(enums.ErrorHandlingOptions.reject_all),
    db: Session = Depends(get_db),
):
    csv_content = csv_file.file.read().decode("utf-8")
    registrar = CompoundRegistrar(db=db, mapping=mapping, error_handling=error_handling)
    rows = registrar.process_csv(csv_content)
    registrar.register_all(rows)
    return registrar.result()


@app.get("/v1/compounds/", response_model=List[models.CompoundResponse])
def read_compounds_v1(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    compounds = crud.get_compounds_v1(db, skip=skip, limit=limit)
    return compounds


@app.get("/v1/compounds/{compound_id}", response_model=models.CompoundResponse)
def read_compound_v1(compound_id: int, db: Session = Depends(get_db)):
    db_compound = crud.get_compound(db, compound_id=compound_id)
    if db_compound is None:
        raise HTTPException(status_code=404, detail="Compound not found")
    return db_compound


@app.get("/v1/compounds/{compound_id}/synonyms", response_model=List[models.CompoundSynonym])
def read_compound_synonyms_v1(compound_id: int, db: Session = Depends(get_db)):
    compound = crud.get_compound(db, compound_id=compound_id)
    if not compound:
        raise HTTPException(status_code=404, detail="Compound not found")
    return compound.compound_synonyms


@app.get("/v1/compounds/{compound_id}/properties", response_model=List[models.Property])
def read_compound_properties_v1(compound_id: int, db: Session = Depends(get_db)):
    compound = crud.get_compound(db, compound_id=compound_id)
    if not compound:
        raise HTTPException(status_code=404, detail="Compound not found")
    return compound.properties


@app.post("/v1/additions/")
def create_additions(
    additions: Optional[str] = Form(None), file: Optional[UploadFile] = File(None), db: Session = Depends(get_db)
):
    created_additions = []
    input_addtions = []

    if file:
        if file.content_type != "text/csv":
            raise HTTPException(status_code=400, detail="Only CSV files are accepted")
        csv_content = file.file.read().decode("utf-8")
        input_addtions = list(csv.DictReader(io.StringIO(csv_content)))

    elif additions:
        input_addtions = json.loads(additions)["additions"]

    create_if_not_exists(
        model_cls=models.Addition,
        base_cls=models.AdditionBase,
        values=input_addtions,
        create_fn=crud.create_addition,
        created_list=created_additions,
        db=db,
    )

    return {
        "status": "success",
        "created": {
            "additions": created_additions,
        },
    }


@app.get("/v1/additions/", response_model=List[models.AdditionResponse])
def read_additions_v1(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    additions = crud.get_additions_v1(db, skip=skip, limit=limit)
    return additions


@app.get("/v1/additions/salts", response_model=List[models.AdditionResponse])
def read_additions_salts_v1(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    salts = crud.get_additions_v1(db, skip=skip, limit=limit, role=enums.AdditionsRole.SALT)
    return salts


@app.get("/v1/additions/solvates", response_model=List[models.AdditionResponse])
def read_additions_solvates_v1(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    solvents = crud.get_additions_v1(db, skip=skip, limit=limit, role=enums.AdditionsRole.SOLVATE)
    return solvents


@app.get("/v1/additions/{addition_id}", response_model=models.AdditionBase)
def read_addition_v1(addition_id: int, db: Session = Depends(get_db)):
    db_addition = crud.get_addition_v1(db, addition_id=addition_id)
    return db_addition


@app.post("/v1/batches/")
def register_batches_v1(
    csv_file: UploadFile = File(...),
    mapping: Optional[str] = Form(None),
    error_handling: enums.ErrorHandlingOptions = Form(enums.ErrorHandlingOptions.reject_all),
    db: Session = Depends(get_db),
):
    csv_content = csv_file.file.read().decode("utf-8")
    registrar = BatchRegistrar(db=db, mapping=mapping, error_handling=error_handling)
    rows = registrar.process_csv(csv_content)
    registrar.register_all(rows)
    return registrar.result()


@app.get("/v1/batches/", response_model=List[models.BatchResponse])
def read_batches_v1(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    batches = crud.get_batches(db, skip=skip, limit=limit)
    return batches


@app.get("/v1/batches/{batch_id}", response_model=models.BatchResponse)
def read_batch_v1(batch_id: int, db: Session = Depends(get_db)):
    db_batch = crud.get_batch(db, batch_id=batch_id)
    if db_batch is None:
        raise HTTPException(status_code=404, detail="Batch not found")
    return db_batch


@app.get("/v1/batches/{batch_id}/synonyms", response_model=List[models.BatchSynonym])
def read_batch_synonyms_v1(batch_id: int, db: Session = Depends(get_db)):
    batch = crud.get_batch(db, batch_id=batch_id)
    if not batch:
        raise HTTPException(status_code=404, detail="Batch not found")
    return batch.batch_synonyms
