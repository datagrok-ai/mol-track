import csv
import io
from fastapi import APIRouter, FastAPI, Depends, File, Form, HTTPException, UploadFile
from sqlalchemy import insert
from sqlalchemy.orm import Session
from typing import List, Optional, Type
from registration.assay_result_registrar import AssayResultsRegistrar
from registration.assay_run_registrar import AssayRunRegistrar
from registration.batch_registrar import BatchRegistrar
from registration.compound_registrar import CompoundRegistrar
import models
import enums
from services.property_service import PropertyService

from typing import Any
from sqlalchemy.sql import text
from rdkit import Chem

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
router = APIRouter(prefix="/v1")

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


def get_or_raise_exception(get_func, db, id, not_found_msg):
    item = get_func(db, id)
    if not item:
        raise HTTPException(status_code=404, detail=not_found_msg)
    return item


# === Schema endpoints for supplementary data like properties and synonyms ===
# https://github.com/datagrok-ai/mol-track/blob/main/api_design.md#schema---wip
@router.post("/schema/")
def preload_schema(payload: models.SchemaPayload, db: Session = Depends(get_db)):
    created_synonyms = crud.create_properties(db, payload.synonym_types)
    created_properties = crud.create_properties(db, payload.properties)

    return {
        "status": "success",
        "created": {
            "synonym_types": created_synonyms,
            "property_types": created_properties,
        },
    }


@router.get("/schema/compounds", response_model=List[models.PropertyBase])
def get_schema_compounds(db: Session = Depends(get_db)):
    return crud.get_entities_by_scope(db, enums.ScopeClass.COMPOUND)


@router.get("/schema/compounds/synonyms", response_model=List[models.SynonymTypeBase])
def get_schema_compound_synonyms(db: Session = Depends(get_db)):
    return crud.get_entities_by_scope(db, enums.ScopeClass.COMPOUND, crud.get_synonym_id(db))


def fetch_additions(db: Session):
    return (
        db.query(models.Addition)
        .join(models.BatchAddition, models.Addition.id == models.BatchAddition.addition_id)
        .distinct()
        .all()
    )


@router.get("/schema/batches", response_model=models.SchemaBatchResponse)
def get_schema_batches(db: Session = Depends(get_db)):
    properties = crud.get_entities_by_scope(db, enums.ScopeClass.BATCH)
    additions = fetch_additions(db)
    return models.SchemaBatchResponse(properties=properties, additions=additions)


@router.get("/schema/batches/synonyms", response_model=models.SchemaBatchResponse)
def get_schema_batch_synonyms(db: Session = Depends(get_db)):
    synonym_types = crud.get_entities_by_scope(db, enums.ScopeClass.BATCH, crud.get_synonym_id(db))
    additions = fetch_additions(db)
    return models.SchemaBatchResponse(synonym_types=synonym_types, additions=additions)


# === Compounds endpoints ===
# https://github.com/datagrok-ai/mol-track/blob/main/api_design.md#register-virtual-compounds
def process_registration(
    registrar_class: Type,
    csv_file: UploadFile,
    mapping: Optional[str],
    error_handling,
    output_format,
    db: Session,
):
    csv_content = csv_file.file.read().decode("utf-8")
    registrar = registrar_class(db=db, mapping=mapping, error_handling=error_handling)
    rows = registrar.process_csv(csv_content)
    registrar.register_all(rows)
    return registrar.result(output_format=output_format)


@router.post("/compounds/")
def register_compounds(
    csv_file: UploadFile = File(...),
    mapping: Optional[str] = Form(None),
    error_handling: enums.ErrorHandlingOptions = Form(enums.ErrorHandlingOptions.reject_all),
    output_format: enums.OutputFormat = Form(enums.OutputFormat.json),
    db: Session = Depends(get_db),
):
    return process_registration(CompoundRegistrar, csv_file, mapping, error_handling, output_format, db)


@router.get("/compounds/", response_model=List[models.CompoundResponse])
def get_compounds(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    compounds = crud.read_compounds(db, skip=skip, limit=limit)
    return compounds


@router.get("/compounds/{compound_id}", response_model=models.CompoundResponse)
def get_compound_by_id(compound_id: int, db: Session = Depends(get_db)):
    return get_or_raise_exception(crud.get_compound_by_id, db, compound_id, "Compound not found")


@router.get("/compounds/{compound_id}/synonyms", response_model=List[models.PropertyWithValue])
def get_compound_synonyms(compound_id: int, db: Session = Depends(get_db)):
    compound = get_or_raise_exception(crud.get_compound_by_id, db, compound_id, "Compound not found")
    return [prop for prop in compound.properties if prop.semantic_type_id == crud.get_synonym_id(db)]


@router.get("/compounds/{compound_id}/properties", response_model=List[models.PropertyWithValue])
def get_compound_properties(compound_id: int, db: Session = Depends(get_db)):
    compound = get_or_raise_exception(crud.get_compound_by_id, db, compound_id, "Compound not found")
    return compound.properties


@router.delete("/compounds/{compound_id}", response_model=models.Compound)
def delete_compound_by_id(compound_id: int, db: Session = Depends(get_db)):
    batches = crud.get_batches_by_compound(db, compound_id=compound_id)
    if batches:
        raise HTTPException(status_code=400, detail="Compound has dependent batches")
    return crud.delete_compound(db, compound_id=compound_id)


# TODO: Create the utils module and move there
def clean_empty_values(d: dict) -> dict:
    return {k: (None if isinstance(v, str) and v.strip() == "" else v) for k, v in d.items()}


# === Additions endpoints ===
# https://github.com/datagrok-ai/mol-track/blob/main/api_design.md#additions
@router.post("/additions/")
def create_additions(file: Optional[UploadFile] = File(None), db: Session = Depends(get_db)):
    if not file:
        raise HTTPException(status_code=400, detail="CSV file is required.")

    if file.content_type != "text/csv":
        raise HTTPException(status_code=400, detail="Only CSV files are accepted.")

    try:
        content = file.file.read().decode("utf-8")
        reader = csv.DictReader(io.StringIO(content))
        input_additions = [models.AdditionBase.model_validate(clean_empty_values(row)) for row in reader]
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Error parsing CSV: {str(e)}")

    created_additions = crud.create_additions(db=db, additions=input_additions)

    return {
        "status": "success",
        "created": {
            "additions": created_additions,
        },
    }


@router.get("/additions/", response_model=List[models.Addition])
def read_additions_v1(db: Session = Depends(get_db)):
    return crud.get_additions(db)


@router.get("/additions/salts", response_model=List[models.Addition])
def read_additions_salts_v1(db: Session = Depends(get_db)):
    return crud.get_additions(db, role=enums.AdditionsRole.SALT)


@router.get("/additions/solvates", response_model=List[models.Addition])
def read_additions_solvates_v1(db: Session = Depends(get_db)):
    return crud.get_additions(db, role=enums.AdditionsRole.SOLVATE)


@router.get("/additions/{addition_id}", response_model=models.Addition)
def read_addition_v1(addition_id: int, db: Session = Depends(get_db)):
    return crud.get_addition_by_id(db, addition_id=addition_id)


@router.put("/additions/{addition_id}", response_model=models.Addition)
def update_addition_v1(addition_id: int, addition_update: models.AdditionUpdate, db: Session = Depends(get_db)):
    return crud.update_addition_by_id(db, addition_id, addition_update)


@router.delete("/additions/{addition_id}", response_model=models.Addition)
def delete_addition(addition_id: int, db: Session = Depends(get_db)):
    return crud.delete_addition_by_id(db, addition_id=addition_id)


# === Batches endpoints ===
# https://github.com/datagrok-ai/mol-track/blob/main/api_design.md#register-batches
@router.post("/batches/")
def register_batches_v1(
    csv_file: UploadFile = File(...),
    mapping: Optional[str] = Form(None),
    error_handling: enums.ErrorHandlingOptions = Form(enums.ErrorHandlingOptions.reject_all),
    output_format: enums.OutputFormat = Form(enums.OutputFormat.json),
    db: Session = Depends(get_db),
):
    return process_registration(BatchRegistrar, csv_file, mapping, error_handling, output_format, db)


@router.get("/batches/", response_model=List[models.BatchResponse])
def read_batches_v1(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    batches = crud.get_batches(db, skip=skip, limit=limit)
    return batches


@router.get("/batches/{batch_id}", response_model=models.BatchResponse)
def read_batch_v1(batch_id: int, db: Session = Depends(get_db)):
    return get_or_raise_exception(crud.get_batch, db, batch_id, "Batch not found")


@router.get("/batches/{batch_id}/properties", response_model=List[models.PropertyWithValue])
def read_batch_properties_v1(batch_id: int, db: Session = Depends(get_db)):
    batch = get_or_raise_exception(crud.get_batch, db, batch_id, "Batch not found")
    return batch.properties


@router.get("/batches/{batch_id}/synonyms", response_model=List[models.PropertyWithValue])
def read_batch_synonyms_v1(batch_id: int, db: Session = Depends(get_db)):
    batch = get_or_raise_exception(crud.get_batch, db, batch_id, "Batch not found")
    return [prop for prop in batch.properties if prop.semantic_type_id == crud.get_synonym_id(db)]


@router.get("/batches/{batch_id}/additions", response_model=List[models.BatchAddition])
def read_batch_additions_v1(batch_id: int, db: Session = Depends(get_db)):
    batch = get_or_raise_exception(crud.get_batch, db, batch_id, "Batch not found")
    return batch.batch_additions


# === Assay data endpoints ===
# https://github.com/datagrok-ai/mol-track/blob/main/api_design.md#assay-data-domain
@router.post("/assays")
def create_assays(payload: List[models.AssayCreateBase], db: Session = Depends(get_db)):
    all_properties = {p.name: p for p in crud.get_properties(db)}
    property_service = PropertyService(all_properties)

    assays_to_insert = [
        {"name": assay.name, "created_by": admin_user_id, "updated_by": admin_user_id} for assay in payload
    ]

    stmt = insert(models.Assay).returning(models.Assay.id)
    inserted_ids = [row[0] for row in db.execute(stmt.values(assays_to_insert)).fetchall()]

    detail_records = []
    property_records = []

    for assay_id, assay in zip(inserted_ids, payload):
        entity_ids = {"assay_id": assay_id}
        inserted, updated = property_service.build_details_records(
            models.AssayDetail,
            properties=assay.extra_fields,
            entity_ids=entity_ids,
            scope=enums.ScopeClass.ASSAY,
            include_user_fields=False,
        )
        detail_records.extend(inserted)

        for prop_data in assay.assay_result_properties:
            prop_info = property_service.get_property_info(prop_data.name, enums.ScopeClass.ASSAY_RESULT)
            property_records.append(
                {
                    "assay_id": assay_id,
                    "property_id": prop_info["property"].id,
                    "required": prop_data.required,
                }
            )

    if detail_records:
        db.execute(insert(models.AssayDetail).values(detail_records))
    if property_records:
        db.execute(insert(models.AssayProperty).values(property_records))

    db.commit()
    return {"status": "success", "created": assays_to_insert}


@router.get("/assays/", response_model=list[models.AssayResponse])
def get_assays(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    assays = crud.get_assays(db, skip=skip, limit=limit)
    return assays


@router.get("/assays/{assay_id}", response_model=models.AssayResponse)
def get_assay_by_id(assay_id: int, db: Session = Depends(get_db)):
    db_assay = crud.get_assay(db, assay_id=assay_id)
    if db_assay is None:
        raise HTTPException(status_code=404, detail="Assay not found")
    return db_assay


@router.post("/assay_runs/")
def create_assay_runs(
    csv_file: UploadFile = File(...),
    mapping: Optional[str] = Form(None),
    error_handling: enums.ErrorHandlingOptions = Form(enums.ErrorHandlingOptions.reject_all),
    output_format: enums.OutputFormat = Form(enums.OutputFormat.json),
    db: Session = Depends(get_db),
):
    return process_registration(AssayRunRegistrar, csv_file, mapping, error_handling, output_format, db)


@router.get("/assay_runs/", response_model=list[models.AssayRunResponse])
def get_assay_runs(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    assay_runs = crud.get_assay_runs(db, skip=skip, limit=limit)
    return assay_runs


@router.get("/assay_runs/{assay_run_id}", response_model=models.AssayRunResponse)
def get_assay_run_by_id(assay_run_id: int, db: Session = Depends(get_db)):
    db_assay_run = crud.get_assay_run(db, assay_run_id=assay_run_id)
    if db_assay_run is None:
        raise HTTPException(status_code=404, detail="Assay run found")
    return db_assay_run


@router.post("/assay_results/")
def create_assay_results(
    csv_file: UploadFile = File(...),
    mapping: Optional[str] = Form(None),
    error_handling: enums.ErrorHandlingOptions = Form(enums.ErrorHandlingOptions.reject_all),
    output_format: enums.OutputFormat = Form(enums.OutputFormat.json),
    db: Session = Depends(get_db),
):
    return process_registration(AssayResultsRegistrar, csv_file, mapping, error_handling, output_format, db)


# === Search endpoints ===
# https://github.com/datagrok-ai/mol-track/blob/main/api_design.md#search---wip
@router.post("/search/compounds/exact", response_model=list[models.Compound])
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
@router.post("/search/compounds/structure", response_model=list[models.Compound])
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


@router.post("/search/complex", response_model=list[dict[str, Any]])
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


app.include_router(router)
