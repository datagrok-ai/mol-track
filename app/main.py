from contextlib import asynccontextmanager
import csv
import io
from fastapi import APIRouter, FastAPI, Depends, File, Form, HTTPException, UploadFile
from sqlalchemy import insert
from sqlalchemy.orm import Session
from typing import List, Optional, Type
from app.services.registrars.assay_result_registrar import AssayResultsRegistrar
from app.services.registrars.assay_run_registrar import AssayRunRegistrar
from app.services.registrars.batch_registrar import BatchRegistrar
from app.services.registrars.compound_registrar import CompoundRegistrar
from app import models
from app import crud
from app.utils import enums
from app.services.property_service import PropertyService
from app.services.search.engine import SearchEngine
from app.services.search.search_filter_builder import SearchFilterBuilder

from sqlalchemy.sql import text

from app.utils.logging_utils import logger


# Handle both package imports and direct execution
try:
    # When imported as a package (for tests)
    from . import models
    from .setup.database import SessionLocal
except ImportError:
    # When run directly
    import app.models as models
    from app.setup.database import SessionLocal

# models.Base.metadata.create_all(bind=engine)


@asynccontextmanager
async def lifespan(app: FastAPI):
    db = SessionLocal()
    try:
        get_admin_user(db)
    finally:
        db.close()

    yield


app = FastAPI(title="MolTrack API", description="API for managing chemical compounds and batches", lifespan=lifespan)
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
    admin = db.query(models.User).filter(models.User.email == "admin@datagrok.ai").first()
    if not admin:
        raise Exception("Admin user not found.")
    global admin_user_id
    admin_user_id = admin.id


def get_or_raise_exception(get_func, db, id, not_found_msg):
    item = get_func(db, id)
    if not item:
        raise HTTPException(status_code=404, detail=not_found_msg)
    return item


# === Schema endpoints for supplementary data like properties and synonyms ===
# https://github.com/datagrok-ai/mol-track/blob/main/api_design.md#schema---wip
@router.post("/schema/")
def preload_schema(payload: models.SchemaPayload, db: Session = Depends(get_db)):
    try:
        created_synonyms = crud.create_properties(db, payload.synonym_types)
        created_properties = crud.create_properties(db, payload.properties)
        return {
            "status": "success",
            "synonym_types": created_synonyms,
            "property_types": created_properties,
        }
    except Exception as e:
        db.rollback()
        return {"status": "failed", "error": str(e)}


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
def create_additions(csv_file: Optional[UploadFile] = File(None), db: Session = Depends(get_db)):
    if not csv_file:
        raise HTTPException(status_code=400, detail="CSV file is required.")

    if csv_file.content_type != "text/csv":
        raise HTTPException(status_code=400, detail="Only CSV files are accepted.")

    try:
        content = csv_file.file.read().decode("utf-8")
        reader = csv.DictReader(io.StringIO(content))
        input_additions = [models.AdditionBase.model_validate(clean_empty_values(row)) for row in reader]
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Error parsing CSV: {str(e)}")

    try:
        created_additions = crud.create_additions(db=db, additions=input_additions)
        return {
            "status": "success",
            "additions": created_additions,
        }
    except Exception as e:
        db.rollback()
        return {"status": "failed", "error": str(e)}


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
    dependent_batch_addition = crud.get_batch_addition_for_addition(db, addition_id)
    if dependent_batch_addition is not None:
        raise HTTPException(status_code=400, detail="Addition has dependent batches")
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


@router.delete("/batches/{batch_id}", response_model=models.Batch)
def delete_batch_by_id(batch_id: int, db: Session = Depends(get_db)):
    assay_results = crud.get_all_assay_results_for_batch(db, batch_id)
    if assay_results:
        raise HTTPException(status_code=400, detail="Batch has dependent assay results")
    return crud.delete_batch(db, batch_id)


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


@router.patch("/admin/compound-matching-rule")
def update_compound_matching_rule(
    rule: enums.CompoundMatchingRule = Form(enums.CompoundMatchingRule.ALL_LAYERS), db: Session = Depends(get_db)
):
    """
    Update the compound matching rule.
    """
    try:
        old_value_query = db.execute(text("SELECT value FROM moltrack.settings WHERE name = 'Compound Matching Rule'"))
        old_value = old_value_query.scalar()

        if old_value == rule.value:
            return {"status": "success", "message": f"Compound matching rule is already set to {rule.value}"}

        db.execute(
            text("UPDATE moltrack.settings SET value = :rule WHERE name = 'Compound Matching Rule'"),
            {"rule": rule.value},
        )
        db.commit()
        return {"status": "success", "message": f"Compound matching rule updated from {old_value} to {rule.value}"}
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=f"Error updating compound matching rule: {str(e)}")


@router.patch("/admin/institution-id-pattern")
def update_institution_id_pattern(
    scope: enums.ScopeClassReduced = Form(enums.ScopeClassReduced.BATCH),
    pattern: str = Form(default="DG-{:05d}"),
    db: Session = Depends(get_db),
):
    """
    Update the pattern for generating corporate IDs for compounds or batches.
    """

    EXPECTED_PATTERN = r"^.{0,10}\{\:0?[1-9]d\}.{0,10}$"
    import re

    if not pattern or not re.match(EXPECTED_PATTERN, pattern):
        raise HTTPException(
            status_code=400,
            detail="""Invalid pattern format. 
                    The pattern must contain '{:d}'.
                    You can also use '{:0Nd}' for zero-padded numbers (numbers will be padded with zeros to N digits).,
                    Pattern can also have prefix and postfix, meant for identification of institution.
                    Example: 'DG-{:05d}' for ids in format 'DG-00001', 'DG-00002' etc.""",
        )

    setting_name = "corporate_batch_id" if scope == "BATCH" else "corporate_compound_id"

    try:
        db.execute(
            text("UPDATE moltrack.properties SET pattern = :pattern WHERE name = :setting"),
            {"setting": setting_name, "pattern": pattern},
        )
        db.commit()
        return {
            "status": "success",
            "message": f"Corporate ID pattern for {scope} updated to {pattern}, ids will be looking like {pattern.format(1)}",
        }
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=f"Error updating corporate ID pattern: {str(e)}")


@router.patch("/admin/molregno-sequence-start")
def set_molregno_sequence_start(start_value: int = Form(...), db: Session = Depends(get_db)):
    return seq_start_update(start_value, "moltrack.molregno_seq", db)


@router.patch("/admin/batchregno-sequence-start")
def set_batchregno_sequence_start(start_value: int = Form(...), db: Session = Depends(get_db)):
    return seq_start_update(start_value, "moltrack.batch_regno_seq", db)


def seq_start_update(start_value: int, seq_name, db: Session):
    """
    Set the starting value for the sequence.
    """
    if start_value < 1:
        raise HTTPException(status_code=400, detail="Start value must be greater than 0")

    max_batchregno = db.execute(text(f"SELECT last_value FROM {seq_name}")).scalar_one()

    if start_value <= max_batchregno:
        raise HTTPException(
            status_code=400,
            detail=f"Start value {start_value} must be greater than the current max {seq_name} {max_batchregno}",
        )

    try:
        db.execute(text("SELECT setval(:seq_name, :start_value)"), {"seq_name": seq_name, "start_value": start_value})
        db.commit()
        return {"status": "success", "message": f"The {seq_name} set to {start_value}"}
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=f"Error setting {seq_name}: {str(e)}")


# === Search endpoints ===
# TODO: Maybe we should move this to a separate module?
def advanced_search(request: models.SearchRequest, db: Session = Depends(get_db)):
    """
    Advanced multi-level search endpoint supporting compounds, batches, and assay_results.
    """
    try:
        engine = SearchEngine(db)
        return engine.search(request)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/search/compounds", response_model=models.SearchResponse)
def search_compounds_advanced(output: List[str], filter: Optional[models.Filter] = None, db: Session = Depends(get_db)):
    """
    Endpoint for compound-level searches.

    Automatically sets level to 'compounds' and accepts filter parameters directly.
    """
    request = models.SearchRequest(level="compounds", output=output, filter=filter)
    return advanced_search(request, db)


@router.post("/search/batches", response_model=models.SearchResponse)
def search_batches_advanced(output: List[str], filter: Optional[models.Filter] = None, db: Session = Depends(get_db)):
    """
    Endpoint for batch-level searches.

    Automatically sets level to 'batches' and accepts filter parameters directly.
    """
    request = models.SearchRequest(level="batches", output=output, filter=filter)
    return advanced_search(request, db)


@router.post("/search/assay-results", response_model=models.SearchResponse)
def search_assay_results_advanced(
    output: List[str], filter: Optional[models.Filter] = None, db: Session = Depends(get_db)
):
    """
    Endpoint for assay result-level searches.

    Automatically sets level to 'assay_results' and accepts filter parameters directly.
    """
    request = models.SearchRequest(level="assay_results", output=output, filter=filter)
    return advanced_search(request, db)


@router.post("/search/generate-filter")
def generate_search_filter(expression: str, db: Session = Depends(get_db)):
    try:
        builder = SearchFilterBuilder(db)
        filter = builder.build_filter(expression)

        return filter
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


app.include_router(router)
