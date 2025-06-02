from typing import Any, Dict, List, Optional, Union
from pydantic import validator
from sqlalchemy import Column, DateTime, Enum, CheckConstraint
from sqlmodel import SQLModel, Field, Relationship
from sqlalchemy.sql import func
import enums
import os
from datetime import datetime
import uuid

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
    compound_synonyms: List["CompoundSynonym"] = []
    properties: List["Property"] = []


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


class Compound(CompoundResponseBase, table=True):
    __tablename__ = "compounds"
    __table_args__ = {"schema": DB_SCHEMA}

    molregno: int = Field(nullable=False)
    formula: str = Field(nullable=False)
    hash_mol: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    hash_tautomer: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    hash_canonical_smiles: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    hash_no_stereo_smiles: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    hash_no_stereo_tautomer: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    created_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    updated_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    deleted_at: Optional[datetime] = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now()))
    deleted_by: Optional[uuid.UUID] = Field(default=None)

    batches: List["Batch"] = Relationship(back_populates="compound")
    compound_synonyms: List["CompoundSynonym"] = Relationship(back_populates="compound")
    properties: List["Property"] = Relationship(
        back_populates="compounds",
        link_model=CompoundDetail,
    )


class BatchBase(SQLModel):
    notes: Optional[str] = Field(default=None)
    compound_id: int = Field(foreign_key=f"{DB_SCHEMA}.compounds.id", nullable=False)


class BatchResponseBase(BatchBase):
    id: int = Field(primary_key=True, index=True)
    created_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now(), nullable=False))


class BatchResponse(BatchResponseBase):
    batch_details: List["BatchDetail"] = []
    batch_synonyms: List["BatchSynonym"] = []


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
    batch_synonyms: List["BatchSynonym"] = Relationship(back_populates="batch")


class SemanticTypeBase(SQLModel):
    name: str = Field(nullable=False)
    description: Optional[str] = Field(default=None, nullable=True)


class SemanticType(SemanticTypeBase, table=True):
    __tablename__ = "semantic_types"
    __table_args__ = {"schema": DB_SCHEMA}

    id: int = Field(primary_key=True, nullable=False)


class PropertyBase(SQLModel):
    name: str = Field(nullable=False)
    value_type: enums.ValueType = Field(sa_column=Column(Enum(enums.ValueType), nullable=False))
    property_class: enums.PropertyClass = Field(sa_column=Column(Enum(enums.PropertyClass), nullable=False))
    unit: Optional[str] = Field(default=None)
    scope: enums.ScopeClass = Field(sa_column=Column(Enum(enums.ScopeClass), nullable=False))
    semantic_type_id: Optional[int] = Field(foreign_key=f"{DB_SCHEMA}.semantic_types.id", nullable=True, default=None)


class PropertyResponse(PropertyBase):
    id: int = Field(primary_key=True, index=True)
    created_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now()))


class AssayTypeProperty(SQLModel, table=True):
    __tablename__ = "assay_type_properties"
    __table_args__ = {"schema": DB_SCHEMA}

    assay_type_id: int = Field(foreign_key=f"{DB_SCHEMA}.assay_types.id", primary_key=True)
    property_id: int = Field(foreign_key=f"{DB_SCHEMA}.properties.id", primary_key=True)
    required: bool = Field(default=False)

    assay_type: "AssayType" = Relationship(back_populates="property_requirements")
    property: "Property" = Relationship()


class AssayTypeBase(SQLModel):
    name: str = Field(nullable=False)
    description: Optional[str] = Field(default=None)


class AssayTypeCreate(AssayTypeBase):
    property_ids: List[int] = []  # List of property IDs to associate with this assay type
    property_requirements: List[Dict[str, Any]] = []  # List of property requirements
    property_details: List[Dict[str, Any]] = []  # List of property metadata


class AssayTypeResponseBase(AssayTypeBase):
    id: int = Field(primary_key=True, index=True)
    created_at: datetime = Field(sa_column=Column(DateTime, server_default=func.now()))
    updated_at: datetime = Field(sa_column=Column(DateTime, server_default=func.now(), onupdate=func.now()))


class AssayTypeResponse(AssayTypeResponseBase):
    properties: List["Property"] = []
    assay_type_details: List["AssayTypeDetail"] = []
    property_requirements: List["AssayTypeProperty"] = []


class AssayType(AssayTypeResponseBase, table=True):
    __tablename__ = "assay_types"
    __table_args__ = {"schema": DB_SCHEMA}

    created_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    updated_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)

    properties: List["Property"] = Relationship(back_populates="assay_types", link_model=AssayTypeProperty)

    assays: List["Assay"] = Relationship(back_populates="assay_type")
    assay_type_details: List["AssayTypeDetail"] = Relationship(back_populates="assay_type")
    property_requirements: List["AssayTypeProperty"] = Relationship(back_populates="assay_type")


class AssayBase(SQLModel):
    name: str = Field(nullable=False)
    description: Optional[str] = Field(default=None)
    assay_type_id: int = Field(foreign_key=f"{DB_SCHEMA}.assay_types.id", nullable=False)


class AssayResponseBase(AssayBase):
    id: int = Field(primary_key=True, index=True)
    created_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now()))


class AssayResponse(AssayResponseBase):
    assay_type: Optional["AssayType"] = None
    assay_details: List["AssayDetail"] = []
    properties: List["Property"] = []


class AssayCreate(AssayBase):
    property_ids: List[int] = []  # List of property IDs to associate with this assay


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
    property: Optional["Property"] = Relationship()


class Assay(AssayResponseBase, table=True):
    __tablename__ = "assays"
    __table_args__ = {"schema": DB_SCHEMA}

    updated_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now()))
    created_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    updated_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)

    # Relationships - use assay_type_properties via assay_type to get list of expected properties
    # No direct properties relationship as assay_properties table no longer exists
    assay_type: "AssayType" = Relationship(back_populates="assays")
    assay_results: List["AssayResult"] = Relationship(back_populates="assay")
    assay_details: List["AssayDetail"] = Relationship(back_populates="assay")

    properties: List["Property"] = Relationship(
        back_populates="assays", link_model=AssayDetail, sa_relationship_kwargs={"lazy": "joined"}
    )


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


class Property(PropertyResponse, table=True):
    __tablename__ = "properties"
    __table_args__ = (
        CheckConstraint(
            "value_type IN ('int', 'double', 'bool', 'datetime', 'string')", name="properties_value_type_check"
        ),
        CheckConstraint(
            "property_class IN ('CALCULATED', 'MEASURED', 'PREDICTED')", name="properties_property_class_check"
        ),
        CheckConstraint("scope IN ('BATCH', 'COMPOUND', 'ASSAY', 'SYSTEM')", name="properties_scope_check"),
        {"schema": DB_SCHEMA},
    )

    updated_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now()))
    created_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    updated_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)

    assay_results: List["AssayResult"] = Relationship(back_populates="property")
    batch_details: List["BatchDetail"] = Relationship(back_populates="property")
    assay_types: List["AssayType"] = Relationship(link_model=AssayTypeProperty)
    assays: List["Assay"] = Relationship(back_populates="properties", link_model=AssayDetail)
    compounds: List["Compound"] = Relationship(
        back_populates="properties",
        link_model=CompoundDetail,
    )


class AssayResultBase(SQLModel):
    batch_id: int = Field(foreign_key=f"{DB_SCHEMA}.batches.id", nullable=False)
    assay_id: int = Field(foreign_key=f"{DB_SCHEMA}.assays.id", nullable=False)
    property_id: int = Field(foreign_key=f"{DB_SCHEMA}.properties.id", nullable=False)

    value_qualifier: Optional[int] = Field(default=0, nullable=False)  # 0 for "=", 1 for "<", 2 for ">"
    value_num: Optional[float] = Field(default=None)
    value_string: Optional[str] = Field(default=None)
    value_bool: Optional[bool] = Field(default=None)


class AssayResultResponseBase(AssayResultBase):
    id: int = Field(primary_key=True, index=True)


# Extended response model with backward compatibility field
class AssayResultResponse(AssayResultResponseBase):
    # Add a computed field for backward compatibility
    result_value: Optional[Union[float, str, bool]] = None

    @validator("result_value", always=True)
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
    assay: "Assay" = Relationship(back_populates="assay_results")
    property: "Property" = Relationship(back_populates="assay_results")


class AssayTypeDetail(SQLModel, table=True):
    __tablename__ = "assay_type_details"
    __table_args__ = {"schema": DB_SCHEMA}

    assay_type_id: int = Field(foreign_key=f"{DB_SCHEMA}.assay_types.id", primary_key=True)
    property_id: int = Field(foreign_key=f"{DB_SCHEMA}.properties.id", primary_key=True)

    value_datetime: Optional[datetime] = Field(default=None, sa_column=Column(DateTime(timezone=True)))
    value_uuid: Optional[uuid.UUID] = Field(default_factory=uuid.uuid4)
    value_num: Optional[float] = Field(default=None)
    value_string: Optional[str] = Field(default=None)

    assay_type: "AssayType" = Relationship(back_populates="assay_type_details")
    property: "Property" = Relationship()


class BatchAssayResultsCreate(SQLModel):
    assay_id: int
    batch_id: int
    measurements: Dict[str, Union[float, str, bool, Dict[str, Any]]]


class BatchAssayResultsResponse(SQLModel):
    assay_id: int
    batch_id: int
    assay_name: str
    measurements: Dict[str, Union[float, str, bool, Dict[str, Any]]]


class AdditionResponse(SQLModel):
    id: int
    name: str


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

    class Config:
        populate_by_name = True


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


class SynonymTypeBase(SQLModel):
    synonym_level: enums.SynonymLevel = Field(
        sa_column=Column(Enum(enums.SynonymLevel)),
        schema_extra={
            "validation_alias": "level",
        },
    )
    name: str = Field(nullable=False)
    pattern: Optional[str] = None
    description: str = Field(default="")

    class Config:
        populate_by_name = True


class SynonymTypeResponse(SQLModel):
    id: int
    name: str


class SynonymType(SynonymTypeBase, table=True):
    __tablename__ = "synonym_types"
    __table_args__ = (
        CheckConstraint("synonym_level IN ('BATCH', 'COMPOUND')", name="synonym_types_synonym_level_check"),
        {"schema": DB_SCHEMA},
    )

    id: int = Field(primary_key=True, index=True)
    created_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now(), nullable=False))
    updated_at: datetime = Field(
        sa_column=Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False)
    )
    created_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    updated_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)


class CompoundSynonymBase(SQLModel):
    compound_id: int = Field(foreign_key="moltrack.compounds.id", nullable=False)
    synonym_type_id: int = Field(foreign_key="moltrack.synonym_types.id", nullable=False)
    synonym_value: str = Field(nullable=False)


class CompoundSynonym(CompoundSynonymBase, table=True):
    __tablename__ = "compound_synonyms"
    __table_args__ = {"schema": DB_SCHEMA}

    id: int = Field(primary_key=True, index=True)
    created_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now(), nullable=False))
    updated_at: datetime = Field(
        sa_column=Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False)
    )
    created_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    updated_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)

    compound: "Compound" = Relationship(back_populates="compound_synonyms")


class BatchSynonymBase(SQLModel):
    batch_id: int = Field(foreign_key="moltrack.batches.id", nullable=False)
    synonym_type_id: int = Field(foreign_key="moltrack.synonym_types.id", nullable=False)
    synonym_value: str = Field(nullable=False)


class BatchSynonym(BatchSynonymBase, table=True):
    __tablename__ = "batch_synonyms"
    __table_args__ = {"schema": DB_SCHEMA}

    id: int = Field(primary_key=True, index=True)
    created_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now(), nullable=False))
    updated_at: datetime = Field(
        sa_column=Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False)
    )
    created_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)
    updated_by: uuid.UUID = Field(nullable=False, default_factory=uuid.uuid4)

    batch: "Batch" = Relationship(back_populates="batch_synonyms")


class SchemaPayload(SQLModel):
    properties: List["PropertyBase"] = Field(default_factory=list)
    synonym_types: List["SynonymTypeBase"] = Field(default_factory=list)


class AdditionsPayload(SQLModel):
    additions: List["AdditionBase"] = Field(default_factory=list)


class SchemaCompoundResponse(SQLModel):
    properties: List["PropertyBase"] = Field(default_factory=list)
    synonym_types: List["SynonymTypeBase"] = Field(default_factory=list)


class SchemaBatchResponse(SchemaCompoundResponse):
    additions: List["AdditionBase"] = Field(default_factory=list)
