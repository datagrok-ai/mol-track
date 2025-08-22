from typing import Any, Dict, List, NamedTuple, Optional, Union, Literal, get_args
from pydantic import ConfigDict, field_validator, model_validator
from sqlalchemy import Column, DateTime, Enum, CheckConstraint
from sqlmodel import SQLModel, Field, Relationship
from sqlalchemy.sql import func
from app.utils import enums
import os
from datetime import datetime
import uuid
# import crud
# # Handle both package imports and direct execution
# try:
#     # When imported as a package (for tests)
#     from .database import Base
# except ImportError:
#     # When run directly
#     from database import Base

# Get the schema name from environment variable or use default
DB_SCHEMA = os.environ.get("DB_SCHEMA", "moltrack")


class User(SQLModel, table=True):
    __tablename__ = "users"
    __table_args__ = {"schema": DB_SCHEMA}

    id: uuid.UUID = Field(primary_key=True, nullable=False, default_factory=uuid.uuid4)
    email: str = Field(nullable=False)
    first_name: str = Field(nullable=False)
    last_name: str = Field(nullable=False)
    has_password: bool = Field(nullable=False)
    is_active: bool = Field(default=False, nullable=False)
    is_service_account: bool = Field(default=False, nullable=False)
    created_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now(), nullable=False))
    updated_at: datetime = Field(
        sa_column=Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False)
    )
    created_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    updated_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)


class Settings(SQLModel, table=True):
    __tablename__ = "settings"
    __table_args__ = {"schema": DB_SCHEMA}

    id: int = Field(primary_key=True, nullable=False)
    name: str = Field(nullable=False)
    value: str = Field(nullable=False)
    description: str = Field(nullable=False)


class CompoundQueryParams(SQLModel):
    substructure: Optional[str] = None
    skip: int = 0
    limit: int = 100


class CompoundBase(SQLModel):
    original_molfile: Optional[str] = Field(default=None)  # as sketched by the chemist
    is_archived: Optional[bool] = Field(default=False)


class CompoundCreate(CompoundBase):
    smiles: str


class CompoundResponseBase(CompoundBase):
    id: int = Field(primary_key=True, index=True)
    canonical_smiles: str = Field(nullable=False)
    inchi: str = Field(nullable=False)
    inchikey: str = Field(nullable=False, unique=True)
    created_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now()))
    updated_at: datetime = Field(
        sa_column=Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now())
    )


class CompoundResponse(CompoundResponseBase):
    batches: List["Batch"] = []
    properties: Optional[List["PropertyWithValue"]] = []


class CompoundDetailBase(SQLModel):
    compound_id: int = Field(foreign_key="moltrack.compounds.id", nullable=False)
    property_id: int = Field(foreign_key="moltrack.properties.id", nullable=False)


class CompoundDetailCreate(CompoundDetailBase):
    value: Any


class CompoundDetail(CompoundDetailBase, table=True):
    __tablename__ = "compound_details"
    __table_args__ = {"schema": DB_SCHEMA}

    id: int = Field(primary_key=True, index=True)
    created_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now(), nullable=False))
    updated_at: datetime = Field(
        sa_column=Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False)
    )
    created_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    updated_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)

    value_datetime: Optional[datetime] = Field(sa_column=Column(DateTime(timezone=True)))
    value_uuid: Optional[uuid.UUID] = Field(default_factory=uuid.uuid4)
    value_num: Optional[float]
    value_string: Optional[str]

    compound: "Compound" = Relationship(back_populates="compound_details")
    property: "Property" = Relationship(back_populates="compound_details")


class Compound(CompoundResponseBase, table=True):
    __tablename__ = "compounds"
    __table_args__ = {"schema": DB_SCHEMA}

    molregno: int = Field(nullable=False)
    formula: str = Field(nullable=False)
    hash_mol: str = Field(nullable=False)
    hash_tautomer: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    hash_canonical_smiles: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    hash_no_stereo_smiles: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    hash_no_stereo_tautomer: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    created_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    updated_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    deleted_at: Optional[datetime] = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now()))
    deleted_by: Optional[uuid.UUID] = Field(default=None)

    batches: List["Batch"] = Relationship(back_populates="compound")
    compound_details: List["CompoundDetail"] = Relationship(back_populates="compound")
    properties: List["Property"] = Relationship(
        back_populates="compounds", link_model=CompoundDetail, sa_relationship_kwargs={"viewonly": True}
    )


class CompoundDetailUpdate(SQLModel):
    property_id: int
    value: Any


class CompoundUpdate(SQLModel):
    original_molfile: Optional[str] = None
    is_archived: Optional[bool] = None
    canonical_smiles: Optional[str] = None
    properties: Optional[List[CompoundDetailUpdate]] = None


class BatchDetailBase(SQLModel):
    batch_id: int = Field(foreign_key=f"{DB_SCHEMA}.batches.id", nullable=False)
    property_id: int = Field(foreign_key=f"{DB_SCHEMA}.properties.id", nullable=False)

    value_qualifier: Optional[int] = Field(default=0, nullable=False)  # 0 for "=", 1 for "<", 2 for ">"
    value_datetime: Optional[datetime] = Field(default=None, sa_column=Column(DateTime(timezone=True)))
    value_uuid: Optional[uuid.UUID] = Field(default_factory=uuid.uuid4)
    value_num: Optional[float] = Field(default=None)
    value_string: Optional[str] = Field(default=None)


class BatchDetailResponse(BatchDetailBase):
    id: int = Field(primary_key=True, index=True)


class BatchDetail(BatchDetailResponse, table=True):
    __tablename__ = "batch_details"
    __table_args__ = {"schema": DB_SCHEMA}

    created_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now(), nullable=False))
    updated_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now(), nullable=False))
    created_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    updated_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)

    batch: "Batch" = Relationship(back_populates="batch_details")
    property: "Property" = Relationship(back_populates="batch_details")


class BatchBase(SQLModel):
    notes: Optional[str] = Field(default=None)
    compound_id: int = Field(foreign_key=f"{DB_SCHEMA}.compounds.id", nullable=False)


class BatchResponseBase(BatchBase):
    id: int = Field(primary_key=True, index=True)
    created_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now(), nullable=False))


class BatchResponse(BatchResponseBase):
    batch_additions: List["BatchAddition"] = []
    properties: Optional[List["PropertyWithValue"]] = []


class Batch(BatchResponseBase, table=True):
    __tablename__ = "batches"
    __table_args__ = {"schema": DB_SCHEMA}

    updated_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now(), nullable=False))
    created_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    updated_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    batch_regno: int = Field(nullable=False)

    compound: "Compound" = Relationship(back_populates="batches")
    assay_results: List["AssayResult"] = Relationship(back_populates="batch")
    batch_details: List["BatchDetail"] = Relationship(back_populates="batch")
    batch_additions: List["BatchAddition"] = Relationship(back_populates="batch")
    properties: List["Property"] = Relationship(
        back_populates="batches", link_model=BatchDetail, sa_relationship_kwargs={"viewonly": True}
    )


class SemanticTypeBase(SQLModel):
    name: str = Field(nullable=False)
    description: Optional[str] = Field(default=None, nullable=True)


class SemanticType(SemanticTypeBase, table=True):
    __tablename__ = "semantic_types"
    __table_args__ = {"schema": DB_SCHEMA}

    id: int = Field(primary_key=True, nullable=False)


class PropertyBase(SQLModel):
    name: str = Field(nullable=False)
    description: Optional[str] = Field(default=None)
    value_type: enums.ValueType = Field(sa_column=Column(Enum(enums.ValueType), nullable=False))
    property_class: enums.PropertyClass = Field(sa_column=Column(Enum(enums.PropertyClass), nullable=False))
    unit: Optional[str] = Field(default=None)
    entity_type: enums.EntityType = Field(sa_column=Column(Enum(enums.EntityType), nullable=False))
    semantic_type_id: Optional[int] = Field(foreign_key=f"{DB_SCHEMA}.semantic_types.id", nullable=True, default=None)
    pattern: Optional[str] = Field(default=None)
    friendly_name: Optional[str] = Field(default=None)


class PropertyWithValue(PropertyBase):
    value_datetime: Optional[datetime] = None
    value_uuid: Optional[uuid.UUID] = None
    value_num: Optional[float] = None
    value_string: Optional[str] = None


class SynonymTypeBase(PropertyBase):
    """A specialized type of PropertyBase for synonyms with fixed constraints"""

    value_type: enums.ValueType = Field(default=enums.ValueType.string)
    property_class: enums.PropertyClass = Field(default=enums.PropertyClass.DECLARED)
    unit: Optional[str] = Field(default="")
    semantic_type_id: Optional[int] = Field(default=1)

    model_config = ConfigDict(populate_by_name=True)

    @field_validator("value_type")
    def validate_value_type(cls, v):
        if v != enums.ValueType.string:
            raise ValueError("SynonymType must have value_type of string")
        return v

    @field_validator("property_class")
    def validate_property_class(cls, v):
        if v != enums.PropertyClass.DECLARED:
            raise ValueError("SynonymType must have property_class of DECLARED")
        return v

    @field_validator("unit")
    def validate_unit(cls, v):
        if v not in (None, ""):
            raise ValueError("SynonymType must have empty unit")
        return v


class PropertyResponse(PropertyBase):
    id: int = Field(primary_key=True, index=True)
    created_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now()))


class AssayProperty(SQLModel, table=True):
    __tablename__ = "assay_properties"
    __table_args__ = {"schema": DB_SCHEMA}

    assay_id: int = Field(foreign_key=f"{DB_SCHEMA}.assays.id", primary_key=True)
    property_id: int = Field(foreign_key=f"{DB_SCHEMA}.properties.id", primary_key=True)
    required: bool = Field(default=False)

    assay: "Assay" = Relationship(back_populates="property_requirements")
    property: "Property" = Relationship()


class AssayBase(SQLModel):
    name: str = Field(nullable=False)
    description: Optional[str] = Field(default=None)


class AssayCreate(AssayBase):
    property_ids: List[int] = []  # List of property IDs to associate with this assay type
    property_requirements: List[Dict[str, Any]] = []  # List of property requirements
    property_details: List[Dict[str, Any]] = []  # List of property metadata


class AssayResponseBase(AssayBase):
    id: int = Field(primary_key=True, index=True)
    created_at: datetime = Field(sa_column=Column(DateTime, server_default=func.now()))
    updated_at: datetime = Field(sa_column=Column(DateTime, server_default=func.now(), onupdate=func.now()))


class AssayResponse(AssayResponseBase):
    properties: List["Property"] = []
    assay_details: List["AssayDetail"] = []
    property_requirements: List["AssayProperty"] = []


class Assay(AssayResponseBase, table=True):
    __tablename__ = "assays"
    __table_args__ = {"schema": DB_SCHEMA}

    created_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    updated_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)

    properties: List["Property"] = Relationship(
        back_populates="assays", link_model=AssayProperty, sa_relationship_kwargs={"viewonly": True}
    )

    assay_runs: List["AssayRun"] = Relationship(back_populates="assay")
    assay_details: List["AssayDetail"] = Relationship(back_populates="assay")
    property_requirements: List["AssayProperty"] = Relationship(back_populates="assay")


class AssayRunBase(SQLModel):
    name: str = Field(nullable=False)
    description: Optional[str] = Field(default=None)
    assay_id: int = Field(foreign_key=f"{DB_SCHEMA}.assays.id", nullable=False)


class AssayRunResponseBase(AssayRunBase):
    id: int = Field(primary_key=True, index=True)
    created_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now()))


class AssayRunResponse(AssayRunResponseBase):
    assay: Optional["Assay"] = None
    assay_run_details: List["AssayRunDetail"] = []
    properties: List["Property"] = []


class AssayRunCreate(AssayRunBase):
    property_ids: List[int] = []  # List of property IDs to associate with this assay


class AssayRunDetail(SQLModel, table=True):
    __tablename__ = "assay_run_details"
    __table_args__ = {"schema": DB_SCHEMA}

    assay_run_id: int = Field(foreign_key=f"{DB_SCHEMA}.assay_runs.id", primary_key=True)
    property_id: int = Field(foreign_key=f"{DB_SCHEMA}.properties.id", primary_key=True)

    value_datetime: Optional[datetime] = Field(default=None, sa_column=Column(DateTime(timezone=True)))
    value_uuid: Optional[uuid.UUID] = Field(default_factory=uuid.uuid4)
    value_num: Optional[float] = Field(default=None)
    value_string: Optional[str] = Field(default=None)

    assay_run: "AssayRun" = Relationship(back_populates="assay_run_details")
    property: Optional["Property"] = Relationship()


class AssayRun(AssayRunResponseBase, table=True):
    __tablename__ = "assay_runs"
    __table_args__ = {"schema": DB_SCHEMA}

    updated_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now()))
    created_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    updated_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)

    # Relationships - use assay_properties via assay to get list of expected properties
    # No direct properties relationship as assay_properties table no longer exists
    assay: "Assay" = Relationship(back_populates="assay_runs")
    assay_results: List["AssayResult"] = Relationship(back_populates="assay_run")
    assay_run_details: List["AssayRunDetail"] = Relationship(back_populates="assay_run")

    properties: List["Property"] = Relationship(
        back_populates="assay_runs",
        link_model=AssayRunDetail,
        sa_relationship_kwargs={"lazy": "joined", "viewonly": True},
    )


class AssayResultDetail(SQLModel, table=True):
    __tablename__ = "assay_result_details"
    __table_args__ = {"schema": DB_SCHEMA}

    assay_result_id: int = Field(foreign_key=f"{DB_SCHEMA}.assay_results.id", primary_key=True)
    property_id: int = Field(foreign_key=f"{DB_SCHEMA}.properties.id", primary_key=True)

    value_qualifier: Optional[int] = Field(default=0, nullable=False)  # 0 for "=", 1 for "<", 2 for ">"
    value_num: Optional[float] = Field(default=None)
    value_string: Optional[str] = Field(default=None)
    value_bool: Optional[bool] = Field(default=None)

    assay_result: "AssayResult" = Relationship(back_populates="assay_result_details")
    property: Optional["Property"] = Relationship()


class Property(PropertyResponse, table=True):
    __tablename__ = "properties"
    __table_args__ = (
        CheckConstraint(
            "value_type IN ('int', 'double', 'bool', 'datetime', 'string', 'uuid')", name="properties_value_type_check"
        ),
        CheckConstraint(
            "property_class IN ('DECLARED', 'CALCULATED', 'MEASURED', 'PREDICTED')",
            name="properties_property_class_check",
        ),
        CheckConstraint(
            "entity_type IN ('BATCH', 'COMPOUND', 'ASSAY', 'ASSAY_TYPES', 'ASSAY_RESULT', 'SYSTEM')",
            name="properties_entity_type_check",
        ),
        {"schema": DB_SCHEMA},
    )

    updated_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now()))
    created_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    updated_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)

    batch_details: List["BatchDetail"] = Relationship(back_populates="property")
    compound_details: List["CompoundDetail"] = Relationship(back_populates="property")

    assays: List["Assay"] = Relationship(link_model=AssayProperty, sa_relationship_kwargs={"viewonly": True})
    assay_runs: List["AssayRun"] = Relationship(
        back_populates="properties", link_model=AssayRunDetail, sa_relationship_kwargs={"viewonly": True}
    )
    assay_results: List["AssayResult"] = Relationship(
        back_populates="properties", link_model=AssayResultDetail, sa_relationship_kwargs={"viewonly": True}
    )
    compounds: List["Compound"] = Relationship(
        back_populates="properties", link_model=CompoundDetail, sa_relationship_kwargs={"viewonly": True}
    )
    batches: List["Batch"] = Relationship(
        back_populates="properties", link_model=BatchDetail, sa_relationship_kwargs={"viewonly": True}
    )


class AssayResultBase(SQLModel):
    batch_id: int = Field(foreign_key=f"{DB_SCHEMA}.batches.id", nullable=False)
    assay_run_id: int = Field(foreign_key=f"{DB_SCHEMA}.assay_runs.id", nullable=False)
    updated_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now()))
    created_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    updated_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)


class AssayResultResponseBase(AssayResultBase):
    id: int = Field(primary_key=True, index=True)


# Extended response model with backward compatibility field
class AssayResultResponse(AssayResultResponseBase):
    # Add a computed field for backward compatibility
    result_value: Optional[Union[float, str, bool]] = None

    @field_validator("result_value")
    def compute_result_value(cls, v, values):
        """Compute result_value from the appropriate typed value field"""
        if "value_num" in values and values["value_num"] is not None:
            return values["value_num"]
        elif "value_string" in values and values["value_string"] is not None:
            return values["value_string"]
        elif "value_bool" in values and values["value_bool"] is not None:
            return values["value_bool"]
        return None


class AssayResult(AssayResultResponseBase, table=True):
    __tablename__ = "assay_results"
    __table_args__ = {"schema": DB_SCHEMA}

    batch: "Batch" = Relationship(back_populates="assay_results")
    assay_run: "AssayRun" = Relationship(back_populates="assay_results")
    assay_result_details: List["AssayResultDetail"] = Relationship(back_populates="assay_result")

    properties: List["Property"] = Relationship(
        back_populates="assay_results",
        link_model=AssayResultDetail,
        sa_relationship_kwargs={"lazy": "joined", "viewonly": True},
    )


class AssayDetail(SQLModel, table=True):
    __tablename__ = "assay_details"
    __table_args__ = {"schema": DB_SCHEMA}

    assay_id: int = Field(foreign_key=f"{DB_SCHEMA}.assays.id", primary_key=True)
    property_id: int = Field(foreign_key=f"{DB_SCHEMA}.properties.id", primary_key=True)

    value_datetime: Optional[datetime] = Field(default=None, sa_column=Column(DateTime(timezone=True)))
    value_uuid: Optional[uuid.UUID] = Field(default_factory=uuid.uuid4)
    value_num: Optional[float] = Field(default=None)
    value_string: Optional[str] = Field(default=None)

    assay: "Assay" = Relationship(back_populates="assay_details")
    property: "Property" = Relationship()


class BatchAssayResultsCreate(SQLModel):
    assay_run_id: int
    batch_id: int
    measurements: Dict[str, Union[float, str, bool, Dict[str, Any]]]


class BatchAssayResultsResponse(SQLModel):
    assay_run_id: int
    batch_id: int
    assay_name: str
    measurements: Dict[str, Union[float, str, bool, Dict[str, Any]]]


class AdditionFields(SQLModel):
    name: Optional[str] = Field(default=None, unique=True)
    description: Optional[str] = None
    code: Optional[str] = None
    smiles: Optional[str] = None
    role: Optional[enums.AdditionsRole] = Field(default=None, sa_column=Column(Enum(enums.AdditionsRole)))
    molfile: Optional[str] = None
    formula: Optional[str] = None
    # Issue related to alias not working properly in SQLModel: https://github.com/fastapi/sqlmodel/issues/374
    molecular_weight: Optional[float] = Field(
        default=None,
        schema_extra={"validation_alias": "molecular weight"},
    )

    model_config = ConfigDict(populate_by_name=True)


class AdditionBase(AdditionFields):
    name: str = Field(nullable=False, unique=True)
    role: enums.AdditionsRole = Field(sa_column=Column(Enum(enums.AdditionsRole)))


class AdditionUpdate(AdditionFields):
    pass


class Addition(AdditionBase, table=True):
    __tablename__ = "additions"
    __table_args__ = (
        CheckConstraint("role IN ('SALT', 'SOLVATE')", name="additions_role_check"),
        {"schema": DB_SCHEMA},
    )

    id: int = Field(primary_key=True, index=True)
    created_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now(), nullable=False))
    updated_at: datetime = Field(
        sa_column=Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False)
    )
    created_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    updated_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    is_active: bool = Field(default=True)
    is_archived: bool = Field(default=False)
    deleted_at: Optional[datetime] = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now()))
    deleted_by: Optional[uuid.UUID] = Field(default=None)


class BatchAdditionBase(SQLModel):
    batch_id: int = Field(foreign_key="moltrack.batches.id", nullable=False, unique=True)
    addition_id: int = Field(foreign_key="moltrack.additions.id", nullable=False, unique=True)
    addition_equivalent: float = Field(default=1)


class BatchAddition(BatchAdditionBase, table=True):
    __tablename__ = "batch_additions"
    __table_args__ = {"schema": DB_SCHEMA}

    id: int = Field(primary_key=True, index=True)
    created_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now(), nullable=False))
    updated_at: datetime = Field(
        sa_column=Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False)
    )
    created_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    updated_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)

    batch: "Batch" = Relationship(back_populates="batch_additions")


class SchemaPayload(SQLModel):
    properties: List["PropertyBase"] = Field(default_factory=list)
    synonym_types: List["SynonymTypeBase"] = Field(default_factory=list)


class AssayPropertiesPayload(SQLModel):
    assay_type_details_properties: List["PropertyBase"] = Field(default_factory=list)
    assay_details_properties: List["PropertyBase"] = Field(default_factory=list)
    assay_type_properties: List["PropertyBase"] = Field(default_factory=list)


class AdditionsPayload(SQLModel):
    additions: List["AdditionBase"] = Field(default_factory=list)


class SchemaBatchResponse(SQLModel):
    properties: Optional[List["PropertyBase"]] = Field(default_factory=list)
    synonym_types: Optional[List["SynonymTypeBase"]] = Field(default_factory=list)
    additions: List["AdditionBase"] = Field(default_factory=list)


class AssayTypeCreateBase(SQLModel):
    name: str
    description: Optional[str] = None
    details: Dict[str, Any]


class AssayTypeCreate(SQLModel):
    assay_type: AssayTypeCreateBase


class AssayResultProperty(SQLModel):
    name: str
    required: bool


class AssayCreateBase(AssayBase):
    assay_result_properties: List[AssayResultProperty]
    extra_fields: Dict[str, Any] = Field(default_factory=dict)

    model_config = {"extra": "allow"}

    @model_validator(mode="before")
    def collect_extra_fields(cls, values):
        known_keys = {"name", "assay_result_properties"}
        extra = {k: v for k, v in values.items() if k not in known_keys}
        values["extra_fields"] = extra
        return values


class UpdateCheckResult(NamedTuple):
    action: str
    update_data: Optional[Dict[str, Any]] = None


def validate_field(v: str) -> str:
    # Basic field format validation
    if not v or not isinstance(v, str):
        raise ValueError("Field must be a non-empty string")

    # Check for valid field format (table.field or table.details.property)
    parts = v.split(".")
    if len(parts) < 2 or len(parts) > 3 or (len(parts) == 3 and parts[1] != "details"):
        raise ValueError("Field must be in format 'table.field' or 'table.details.property'")

    valid_tables = get_args(Level)
    if parts[0] not in valid_tables:
        allowed = ", ".join(valid_tables)
        raise ValueError(f"Invalid table: {parts[0]}. Must be one of {allowed}")

    return v


# Advanced Search Models - New Recursive Structure
class AtomicCondition(SQLModel):
    """Individual atomic search condition with field, operator, and value"""

    field: str  # e.g., "compounds.canonical_smiles", "compounds.details.chembl"
    operator: enums.CompareOp
    value: Any
    threshold: Optional[float] = None  # For similarity searches (e.g., molecular similarity)

    @field_validator("field")
    def validate_field_format(cls, v):
        return validate_field(v)

    @model_validator(mode="before")
    def validate_threshold(cls, values):
        # Validate threshold is provided for operators that require it
        if isinstance(values, AtomicCondition):
            operator = values.get("operator")
            threshold = values.get("threshold")
            if operator == enums.CompareOp.IS_SIMILAR:
                if threshold is None:
                    raise ValueError(f"Operator {operator.value} requires a threshold value")
            elif threshold is not None:
                raise ValueError(f"Threshold not supported for operator: {operator}")
        return values


class LogicalNode(SQLModel):
    """Logical node combining multiple filters with AND/OR operator"""

    operator: enums.LogicOp
    conditions: List["Filter"]

    @field_validator("conditions")
    def validate_conditions(cls, v):
        if not v or len(v) == 0:
            raise ValueError("LogicalNode must have at least one condition")
        if len(v) == 1:
            raise ValueError("LogicalNode with single condition is redundant, use AtomicCondition directly")
        return v


# Filter is a union type - can be either AtomicCondition or LogicalNode
Filter = Union[AtomicCondition, LogicalNode]


Level = Literal["compounds", "batches", "assay_results", "assay_runs", "assays"]


AggregationOp = Union[enums.AggregationNumericOp, enums.AggregationStringOp]


class Aggregation(SQLModel):
    field: str
    operation: AggregationOp

    @field_validator("field", mode="after")
    def validate_field_format(cls, v):
        return validate_field(v)


class SearchRequest(SQLModel):
    """Main search request model with recursive filter structure"""

    level: Level
    output: List[str]  # Columns to return
    aggregations: Optional[List[Aggregation]] = Field(default_factory=list)
    filter: Optional[Filter] = None
    output_format: enums.SearchOutputFormat

    @field_validator("output")
    def validate_output(cls, v):
        if not v:
            raise ValueError("Output must specify at least one column")
        return v


class SearchResponse(SQLModel):
    """Search response model"""

    status: str
    data: List[Dict[str, Any]]
    total_count: int
    level: str
    columns: List[str]


# Update forward references for recursive types
LogicalNode.model_rebuild()


class Token(NamedTuple):
    type: str
    value: Union[str, float, bool, None]
