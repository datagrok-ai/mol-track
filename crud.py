from sqlalchemy.orm import Session
from fastapi import HTTPException
from rdkit import Chem
from typing import List, Optional
from sqlalchemy import insert, text
from sqlalchemy.orm import selectinload
from datetime import datetime
import models as models
from rdkit.Chem import Descriptors, rdMolDescriptors
import main
import enums


from typing import Type, Dict, Any
from rdkit.Chem.RegistrationHash import HashLayer
from chemistry_utils import standardize_mol, generate_hash_layers, generate_uuid_from_string

# Handle both package imports and direct execution
try:
    # When imported as a package (for tests)
    from . import models
except ImportError:
    # When run directly
    import models


# === Compound-related operations ===
def get_compound_by_hash(db: Session, hash_mol: str):
    """
    Search for compounds in the database by hash_mol.
    """
    if not isinstance(hash_mol, str) or len(hash_mol) != 40:
        raise HTTPException(status_code=400, detail=f"Invalid hash_mol format: {hash_mol}")

    return db.query(models.Compound).filter(models.Compound.hash_mol == hash_mol).all()


def enrich_compound(compound: models.Compound) -> models.CompoundResponse:
    return models.CompoundResponse(
        **compound.dict(), properties=enrich_properties(compound, "compound_details", "compound_id")
    )


def read_compounds(db: Session, skip: int = 0, limit: int = 100):
    compounds = db.query(models.Compound).offset(skip).limit(limit).all()
    return [enrich_compound(c) for c in compounds]


def get_compound_by_id(db: Session, compound_id: int):
    compound = (
        db.query(models.Compound)
        .options(selectinload(models.Compound.properties).selectinload(models.Property.compound_details))
        .filter(models.Compound.id == compound_id)
        .first()
    )

    if not compound:
        return None

    return enrich_compound(compound)


def delete_compound(db: Session, compound_id: int):
    db_compound = db.get(models.Compound, compound_id)
    if db_compound is None:
        raise HTTPException(status_code=404, detail="Compound not found")
    db_compound.deleted_at = datetime.now()

    db.commit()
    db.refresh(db_compound)
    return db_compound


def get_compounds_ex(db: Session, query_params: models.CompoundQueryParams):
    """
    Get compounds with optional filtering parameters.

    Args:
        db: Database session
        query_params: Query parameters including substructure, skip, and limit

    Returns:
        List of compounds matching the query parameters
    """
    # If substructure is provided, use substructure search
    if query_params.substructure:
        # Validate the substructure SMILES
        mol = Chem.MolFromSmiles(query_params.substructure)
        if mol is None:
            raise HTTPException(status_code=400, detail="Invalid substructure SMILES string")

        # SQL query using the RDKit cartridge substructure operator '@>'
        sql = text(f"""
            SELECT c.* FROM {models.DB_SCHEMA}.compounds c
            JOIN rdk.mols ON rdk.mols.id = c.id
            WHERE rdk.mols.m@>'{query_params.substructure}'
            ORDER BY c.id
            OFFSET :skip LIMIT :limit
        """)

        result = db.execute(sql, {"skip": query_params.skip, "limit": query_params.limit})
        compounds = []
        for row in result:
            compound = models.Compound()
            for column, value in row._mapping.items():
                setattr(compound, column, value)
            compounds.append(compound)

        return compounds

    else:
        # If no substructure provided, use regular get_compounds function
        return read_compounds(db, skip=query_params.skip, limit=query_params.limit)


# === Batch-related operations ===
def enrich_batch(batch: models.Batch) -> models.BatchResponse:
    return models.BatchResponse(**batch.dict(), properties=enrich_properties(batch, "batch_details", "batch_id"))


def get_batch(db: Session, batch_id: int):
    batch = (
        db.query(models.Batch)
        .options(selectinload(models.Batch.properties).selectinload(models.Property.batch_details))
        .filter(models.Batch.id == batch_id)
        .first()
    )
    if not batch:
        return None
    return enrich_batch(batch)


def get_batches(db: Session, skip: int = 0, limit: int = 100):
    batches = db.query(models.Batch).offset(skip).limit(limit).all()
    return [enrich_batch(batch) for batch in batches]


def get_batches_by_compound(db: Session, compound_id: int, skip: int = 0, limit: int = 100):
    return db.query(models.Batch).filter(models.Batch.compound_id == compound_id).offset(skip).limit(limit).all()


# === Property-related operations ===
def create_properties(db: Session, properties: list[models.PropertyBase]) -> list[dict]:
    return bulk_create_if_not_exists(db, models.Property, models.PropertyBase, properties)


def get_properties(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.Property).offset(skip).limit(limit).all()


# === Assay-related operations ===
def get_assay(db: Session, assay_id: int):
    return db.query(models.Assay).filter(models.Assay.id == assay_id).first()


def get_assays(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.Assay).offset(skip).limit(limit).all()


# === AssayRun-related operations ===
def get_assay_run(db: Session, assay_run_id: int):
    # Get the assay
    assay = db.query(models.AssayRun).filter(models.AssayRun.id == assay_run_id).first()

    if assay:
        # Get properties associated with this assay through assay details
        assay_details = db.query(models.AssayRunDetail).filter(models.AssayRunDetail.assay_run_id == assay_run_id).all()
        property_ids = [detail.property_id for detail in assay_details]

        # If no properties from assay details, get them from the assay type
        if not property_ids and assay.assay:
            assay_properties = (
                db.query(models.AssayProperty).filter(models.AssayProperty.assay_id == assay.assay_id).all()
            )
            property_ids = [prop.property_id for prop in assay_properties]

        # Get the property objects
        if property_ids:
            properties = db.query(models.Property).filter(models.Property.id.in_(property_ids)).all()
            # Add properties to assay
            assay.properties = properties
        else:
            assay.properties = []

    return assay


def get_assay_runs(db: Session, skip: int = 0, limit: int = 100):
    # Get assay_runs with pagination
    assay_runs = db.query(models.AssayRun).offset(skip).limit(limit).all()

    # For each assay, add its properties
    for assay_run in assay_runs:
        # Get properties associated with this assay through assay details
        assay_details = db.query(models.AssayRunDetail).filter(models.AssayRunDetail.assay_run_id == assay_run.id).all()
        property_ids = [detail.property_id for detail in assay_details]

        # If no properties from assay details, get them from the assay type
        if not property_ids and assay_run.assay:
            assay_properties = (
                db.query(models.AssayProperty).filter(models.AssayProperty.assay_id == assay_run.assay_id).all()
            )
            property_ids = [prop.property_id for prop in assay_properties]

        # Get the property objects
        if property_ids:
            properties = db.query(models.Property).filter(models.Property.id.in_(property_ids)).all()
            # Add properties to assay
            assay_run.properties = properties
        else:
            assay_run.properties = []

    return assay_runs


# === AssayResult-related operations ===
def get_assay_result(db: Session, assay_result_id: int):
    return db.query(models.AssayResult).filter(models.AssayResult.id == assay_result_id).first()


def get_assay_results(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.AssayResult).offset(skip).limit(limit).all()


def get_batch_assay_results(db: Session, batch_id: int):
    """Get all assay results for a specific batch"""
    results = db.query(models.AssayResult).filter(models.AssayResult.batch_id == batch_id).all()

    # Group results by assay_run_id
    grouped_results = {}
    for result in results:
        assay_run_id = result.assay_run_id
        if assay_run_id not in grouped_results:
            # Get the assay name
            assay_run = db.query(models.AssayRun).filter(models.AssayRun.id == assay_run_id).first()
            assay = db.query(models.Assay).filter(models.Assay.id == assay_run.assay_id).first()
            assay_name = assay.name if assay else "Unknown Assay"

            grouped_results[assay_run_id] = {
                "assay_run_id": assay_run_id,
                "batch_id": batch_id,
                "assay_name": assay_name,
                "measurements": {},
            }

        # Get property name and type
        property = db.query(models.Property).filter(models.Property.id == result.property_id).first()
        property_name = property.name if property else f"Property-{result.property_id}"
        property_type = property.value_type if property else "double"

        # Get value based on property type
        value = None
        if property_type in ("int", "double"):
            value = result.value_num
        elif property_type == "string":
            value = result.value_string
        elif property_type == "bool":
            value = result.value_bool

        # If we have a qualifier other than "=" (0), include it in the result
        if result.value_qualifier != 0:
            grouped_results[assay_run_id]["measurements"][property_name] = {
                "qualifier": result.value_qualifier,
                "value": value,
            }
        else:
            grouped_results[assay_run_id]["measurements"][property_name] = value

    return list(grouped_results.values())


# === Addition-related operations ===
def enrich_addition(add: models.AdditionBase) -> models.AdditionBase:
    smiles, molfile, formula, mw = add.smiles, add.molfile, add.formula, add.molecular_weight

    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if not mol and molfile:
        try:
            mol = Chem.MolFromMolBlock(molfile, sanitize=True)
        except Exception:
            mol = None

    if not mol:
        return add

    add.smiles = smiles or Chem.MolToSmiles(mol)
    add.molfile = molfile or Chem.MolToMolBlock(mol)
    add.formula = formula or rdMolDescriptors.CalcMolFormula(mol)
    add.molecular_weight = mw or Descriptors.MolWt(mol)

    return add


def create_additions(db: Session, additions: list[models.AdditionBase]) -> list[dict]:
    enriched_additions = [enrich_addition(add) for add in additions]
    return bulk_create_if_not_exists(db, models.Addition, models.AdditionBase, enriched_additions, validate=False)


def get_additions(db: Session, role: enums.AdditionsRole | None = None) -> List[models.AdditionBase]:
    query = db.query(models.Addition)
    if role is None:
        return query.all()
    return query.filter_by(role=role).all()


def get_addition_by_id(db: Session, addition_id: int) -> models.Addition:
    db_addition = db.get(models.Addition, addition_id)
    if db_addition is None:
        raise HTTPException(status_code=404, detail="Addition not found")
    return db_addition


def update_addition_by_id(db: Session, addition_id: int, addition_update: models.AdditionUpdate):
    db_addition = get_addition_by_id(db, addition_id=addition_id)
    update_data = addition_update.dict(exclude_unset=True)
    for key, value in update_data.items():
        setattr(db_addition, key, value)

    db_addition.updated_at = datetime.now()
    db.add(db_addition)
    db.commit()
    db.refresh(db_addition)
    return db_addition


def delete_addition_by_id(db: Session, addition_id: int):
    db_addition = get_addition_by_id(db, addition_id=addition_id)
    db_addition.deleted_at = datetime.now()
    db_addition.is_active = False
    db.commit()
    db.refresh(db_addition)
    return db_addition


# === Helper functions ===
def bulk_create_if_not_exists(
    db: Session,
    model_cls: Type,
    base_model_cls: Type,
    items: List[Any],
    *,
    name_attr: str = "name",
    validate: bool = True,
) -> List[Dict[str, Any]]:
    """
    Bulk insert records into the database for the given SQLModel class, only if records with the same
    unique identifier (specified by `name_attr`) do not already exist.

    This function is designed for efficient batch creation of SQLModel instances,
    validating input items against a base SQLModel class before insertion, and
    adding audit fields such as created_by and updated_by automatically.

    Args:
        db (Session): SQLAlchemy session used to query and insert data.
        model_cls (Type[SQLModel]): The SQLModel ORM model class representing the target table.
        base_model_cls (Type[SQLModel]): The SQLModel base class used for validation and serialization.
        items (List[SQLModel]): List of SQLModel instances to insert.
        name_attr (str, optional): Attribute name used to check uniqueness (default: "name").
        validate (bool, optional): Whether to validate each item using `base_model_cls` before insert (default: True).

    Returns:
        List[Dict[str, Any]]: List of inserted records.
    """
    input_names = [getattr(item, name_attr) for item in items]
    existing_names = {
        name
        for (name,) in db.query(getattr(model_cls, name_attr))
        .filter(getattr(model_cls, name_attr).in_(input_names))
        .all()
    }

    to_insert = []
    for item in items:
        item_name = getattr(item, name_attr)
        if item_name not in existing_names:
            validated = base_model_cls.model_validate(item) if validate else item
            data = validated.model_dump()
            data.update(
                {
                    "created_by": main.admin_user_id,
                    "updated_by": main.admin_user_id,
                }
            )
            to_insert.append(data)

    if not to_insert:
        return []

    stmt = insert(model_cls).values(to_insert).returning(model_cls)
    result = db.execute(stmt).fetchall()
    db.commit()
    return [model_cls.model_validate(row[0]) for row in result]


def get_synonym_id(db: Session) -> int:
    result = db.query(models.SemanticType.id).filter(models.SemanticType.name == "Synonym").scalar()
    if result is None:
        raise HTTPException(status_code=400, detail="Semantic type 'Synonym' not found.")
    return result


def get_entities_by_scope(db: Session, scope: enums.ScopeClass, semantic_type_id: Optional[int] = None):
    query = db.query(models.Property).filter(models.Property.scope == scope)
    if semantic_type_id is not None:
        query = query.filter(models.Property.semantic_type_id == semantic_type_id)
    return query.all()


def enrich_properties(owner, detail_attr: str, id_attr: str) -> list[models.PropertyWithValue]:
    enriched = []
    owner_id = getattr(owner, "id")
    for prop in owner.properties:
        details = getattr(prop, detail_attr, [])
        detail = next((d for d in details if getattr(d, id_attr, None) == owner_id), None)
        enriched.append(
            models.PropertyWithValue(
                **prop.dict(),
                value_num=getattr(detail, "value_num", None),
                value_string=getattr(detail, "value_string", None),
                value_datetime=getattr(detail, "value_datetime", None),
                value_uuid=getattr(detail, "value_uuid", None),
            )
        )
    return enriched


# === Search-related operations ===
def search_compounds_exact(query_smiles: str, search_parameters: models.ExactSearchParameters, db: Session):
    """
    Perform an exact search for compounds.

    - **query_smiles**: The SMILES string to search against.
    - **search_parameters**: Parameters for the exact search.
    """
    if not search_parameters:
        raise HTTPException(status_code=400, detail="Search parameters are required for exact search")
    exact_params = models.ExactSearchParameters(**search_parameters)
    return db.query(db=db, query_smiles=query_smiles, fields=exact_params.field)


def search_compounds_substructure(db: Session, query_smiles: str, search_parameters: Dict[str, Any]):
    """
    Placeholder for substructure search.
    """
    raise NotImplementedError("Substructure search is not implemented yet.")


def search_compounds_similarity(db: Session, query_smiles: str, search_parameters: Dict[str, Any]):
    """
    Placeholder for similarity search.
    """
    raise NotImplementedError("Similarity search is not implemented yet.")


def get_standardized_mol_and_layers(query_smiles: str, http_errors: bool = False) -> dict:
    missing_msg = "Query SMILES is required"
    invalid_msg = "Invalid SMILES string"

    if not query_smiles:
        if http_errors:
            raise HTTPException(status_code=400, detail=missing_msg)
        else:
            raise ValueError(missing_msg)

    mol = Chem.MolFromSmiles(query_smiles)
    if mol is None:
        if http_errors:
            raise HTTPException(status_code=400, detail=invalid_msg)
        else:
            raise ValueError(invalid_msg)

    standardized_mol = standardize_mol(mol)
    mol_layers = generate_hash_layers(standardized_mol)
    return mol_layers


def search_compounds_by_hash(
    db: Session, query_smiles: str, hash_layer: Any, hash_attr_name: str
) -> List[models.Compound]:
    mol_layers = get_standardized_mol_and_layers(query_smiles, http_errors=True)
    hash_value = generate_uuid_from_string(mol_layers[hash_layer])
    return db.query(models.Compound).filter(getattr(models.Compound, hash_attr_name) == hash_value).all()


def search_compounds_tautomer(db: Session, query_smiles: str, search_parameters: Dict[str, Any]):
    """
    Perform a tautomer search for compounds. For a given molecule, find all compounds that have the same tautomer hash regardless of their stereochemistry.
    """
    return search_compounds_by_hash(db, query_smiles, HashLayer.TAUTOMER_HASH, "hash_tautomer")


def search_compounds_stereo(db: Session, query_smiles: str, search_parameters: Dict[str, Any]):
    """
    Perform a stereo search for compounds regardless tautomeric state."""
    return search_compounds_by_hash(db, query_smiles, HashLayer.NO_STEREO_SMILES, "hash_no_stereo_smiles")


def search_compounds_connectivity(db: Session, query_smiles: str, search_parameters: Dict[str, Any]):
    """
    Perform a connectivity search for compounds.
    This will find compounds that have the same connectivity regardless of their stereochemistry or tautomeric state.
    """
    return search_compounds_by_hash(db, query_smiles, HashLayer.NO_STEREO_TAUTOMER_HASH, "hash_no_stereo_tautomer")
