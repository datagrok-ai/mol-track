import os
from datetime import date, datetime
from typing import List, Optional
from uuid import UUID
from sqlmodel import SQLModel, Field, Relationship
from sqlalchemy.sql import func
from sqlalchemy import CheckConstraint

DB_SCHEMA = os.environ.get("DB_SCHEMA", "moltrack")


class CreatedUpdatedDetails(SQLModel, table=False):
    created_at: Optional[datetime] = Field(default=None, sa_column_kwargs={"server_default": func.now()})
    updated_at: Optional[datetime] = Field(
        default=None, sa_column_kwargs={"server_default": func.now(), "onupdate": func.now()}
    )
    created_by: Optional[UUID] = Field(default=None, foreign_key=f"{DB_SCHEMA}.users.id", nullable=True)
    updated_by: Optional[UUID] = Field(default=None, foreign_key=f"{DB_SCHEMA}.users.id", nullable=True)


class Base(CreatedUpdatedDetails, table=False):
    id: int = Field(default=None, primary_key=True)


class SoftDeleteTemplate(SQLModel, table=False):
    is_archived: bool = Field(default=False)
    deleted_at: Optional[datetime] = Field(default=None, sa_column_kwargs={"server_default": func.now()})
    deleted_by: UUID = Field(foreign_key=f"{DB_SCHEMA}.users.id")


class User(CreatedUpdatedDetails, table=True):
    __tablename__ = "users"
    __table_args__ = {"schema": DB_SCHEMA}

    id: UUID = Field(primary_key=True)
    email: str = Field(nullable=False, unique=True)
    first_name: str = Field(nullable=False)
    last_name: str = Field(nullable=False)
    has_password: bool = Field(nullable=False)
    is_active: bool = Field(default=False, nullable=False)
    is_service_account: bool = Field(default=False, nullable=False)


class Settings(SQLModel, table=True):
    __tablename__ = "settings"
    __table_args__ = {"schema": DB_SCHEMA}

    id: Optional[int] = Field(default=None, primary_key=True)
    name: str = Field(nullable=False)
    value: str = Field(nullable=False)
    description: str = Field(nullable=False)


class Compound(Base, SoftDeleteTemplate, table=True):
    __tablename__ = "compounds"
    __table_args__ = {"schema": DB_SCHEMA}

    canonical_smiles: str = Field(nullable=False)
    original_molfile: Optional[str] = None
    molregno: int = Field(nullable=False, unique=True)
    inchi: str = Field(nullable=False)
    inchikey: str = Field(nullable=False, unique=True)
    formula: str = Field(nullable=False)

    hash_mol: UUID = Field(nullable=False, unique=True)
    hash_tautomer: UUID = Field(nullable=False, unique=True)
    hash_canonical_smiles: UUID = Field(nullable=False, unique=True)
    hash_no_stereo_smiles: UUID = Field(nullable=False, unique=True)
    hash_no_stereo_tautomer: UUID = Field(nullable=False, unique=True)

    # Relationships
    # batches: list["Batch"] = Relationship(sa_relationship="Batch", back_populates="compound")


class Batch(Base, table=True):
    __tablename__ = "batches"
    __table_args__ = {"schema": DB_SCHEMA}

    compound_id: int = Field(foreign_key=f"{DB_SCHEMA}.compounds.id", nullable=False)
    batch_number: str = Field(nullable=False)

    amount: Optional[float] = None
    amount_unit: Optional[str] = None
    purity: Optional[float] = None
    notes: Optional[str] = None
    expiry_date: Optional[date] = None

    # Relationships
    compound: "Compound" = Relationship(sa_relationship="Compound", back_populates="batches")
    assay_results: list["AssayResult"] = Relationship(sa_relationship="AssayResult", back_populates="batch")
    batch_details: list["BatchDetail"] = Relationship(sa_relationship="BatchDetail", back_populates="batch")


class Addition(Base, SoftDeleteTemplate, table=True):
    __tablename__ = "additions"
    __table_args__ = {"schema": DB_SCHEMA}

    name: str = Field(nullable=False)
    description: Optional[str] = None
    display_name: Optional[str] = None
    display_order: Optional[int] = None
    is_active: bool = Field(default=True, nullable=False)

    formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    smiles: Optional[str] = None
    molfile: Optional[str] = None
    role: Optional[str] = None


class BatchAddition(Base, table=True):
    __tablename__ = "batch_additions"
    __table_args__ = {"schema": DB_SCHEMA}

    batch_id: int = Field(foreign_key=f"{DB_SCHEMA}.batches.id", primary_key=True)
    addition_id: int = Field(foreign_key=f"{DB_SCHEMA}.additions.id", primary_key=True)
    addition_equivalent: float = Field(default=1.0, nullable=False)


class AssayPropertyLink(SQLModel, table=True):
    __tablename__ = "assay_properties"
    __table_args__ = {"schema": DB_SCHEMA}

    assay_id: int = Field(foreign_key=f"{DB_SCHEMA}.assays.id", primary_key=True)
    property_id: int = Field(foreign_key=f"{DB_SCHEMA}.properties.id", primary_key=True)


class AssayTypePropertyLink(SQLModel, table=True):
    __table_args__ = {"schema": DB_SCHEMA}

    assay_type_id: int = Field(foreign_key=f"{DB_SCHEMA}.assay_types.id", primary_key=True)
    property_id: int = Field(foreign_key=f"{DB_SCHEMA}.properties.id", primary_key=True)


class SemanticType(SQLModel, table=True):
    __tablename__ = "semantic_types"
    __table_args__ = {"schema": DB_SCHEMA}

    id: Optional[int] = Field(default=None, primary_key=True)
    name: str = Field(nullable=False)
    description: Optional[str] = None


class Property(Base, table=True):
    __tablename__ = "properties"
    __table_args__ = (
        CheckConstraint("value_type IN ('int', 'double', 'bool', 'datetime', 'string')", name="check_value_type"),
        CheckConstraint("property_class IN ('CALCULATED', 'MEASURED', 'PREDICTED')", name="check_property_class"),
        {"schema": DB_SCHEMA},
    )

    name: str = Field(nullable=False)
    value_type: Optional[str] = None
    semantic_type_id: Optional[int] = Field(default=None, foreign_key=f"{DB_SCHEMA}.semantic_types.id")
    property_class: Optional[str] = None
    unit: Optional[str] = None

    # Relationships
    assay_types: List["AssayType"] = Relationship(back_populates="properties", link_model=AssayTypePropertyLink)
    assays: List["Assay"] = Relationship(
        sa_relationship="Assay", back_populates="properties", link_model=AssayPropertyLink
    )
    assay_results = Relationship(sa_relationship="AssayResult", back_populates="property")
    batch_details = Relationship(sa_relationship="BatchDetail", back_populates="property")

    class Config:
        orm_mode = True


class BatchDetail(SQLModel, table=True):
    __tablename__ = "batch_details"
    __table_args__ = {"schema": DB_SCHEMA}

    id: Optional[int] = Field(default=None, primary_key=True)
    batch_id: int = Field(foreign_key=f"{DB_SCHEMA}.batches.id", nullable=False)
    property_id: int = Field(foreign_key=f"{DB_SCHEMA}.properties.id", nullable=False)

    value_datetime: Optional[datetime] = Field(default=None, sa_column_kwargs={"server_default": func.now()})
    value_uuid: Optional[UUID] = None
    value_num: Optional[float] = None
    value_string: Optional[str] = None

    # Relationships
    batch: Optional[Batch] = Relationship(sa_relationship="Batch", back_populates="batch_details")
    property: Optional[Property] = Relationship(sa_relationship="Property", back_populates="batch_details")


class SynonymType(Base, table=True):
    __tablename__ = "synonym_types"
    __table_args__ = {"schema": DB_SCHEMA}

    synonym_level: str = Field(nullable=False)
    name: str = Field(nullable=False)
    pattern: str = Field(nullable=False)
    description: str = Field(nullable=False)


class CompoundSynonym(Base, table=True):
    __tablename__ = "compound_synonyms"
    __table_args__ = {"schema": DB_SCHEMA}

    compound_id: int = Field(foreign_key=f"{DB_SCHEMA}.compounds.id", nullable=False)
    synonym_type_id: int = Field(foreign_key=f"{DB_SCHEMA}.synonym_types.id", nullable=False)
    synonym_value: str = Field(nullable=False)


class BatchSynonym(Base, table=True):
    __tablename__ = "batch_synonyms"
    __table_args__ = {"schema": DB_SCHEMA}

    batch_id: int = Field(foreign_key=f"{DB_SCHEMA}.batches.id", nullable=False)
    synonym_type_id: int = Field(foreign_key=f"{DB_SCHEMA}.synonym_types.id", nullable=False)
    synonym_value: str = Field(nullable=False)


class CompoundDetails(CreatedUpdatedDetails, table=True):
    __tablename__ = "compound_details"
    __table_args__ = {"schema": DB_SCHEMA}

    id: Optional[int] = Field(default=None, primary_key=True)
    compound_id: int = Field(foreign_key=f"{DB_SCHEMA}.compounds.id", nullable=False)
    property_id: int = Field(foreign_key=f"{DB_SCHEMA}.properties.id", nullable=False)

    value_datetime: Optional[datetime] = Field(default=None, sa_column_kwargs={"server_default": func.now()})
    value_uuid: Optional[UUID] = None
    value_num: Optional[float] = None
    value_string: Optional[str] = None


class Assay(CreatedUpdatedDetails, table=True):
    __tablename__ = "assays"
    __table_args__ = {"schema": DB_SCHEMA}

    id: Optional[int] = Field(default=None, primary_key=True)
    assay_type_id: int = Field(foreign_key=f"{DB_SCHEMA}.assay_types.id", nullable=False)
    name: str = Field(nullable=False)
    description: Optional[str] = None
    properties: List["Property"] = Relationship(
        sa_relationship="Property", back_populates="assays", link_model=AssayPropertyLink
    )

    assay_type = Relationship(sa_relationship="AssayType", back_populates="assays")
    assay_results = Relationship(sa_relationship="AssayResult", back_populates="assay")
    assay_details = Relationship(sa_relationship="AssayDetail", back_populates="assay")


class AssayType(CreatedUpdatedDetails, table=True):
    __tablename__ = "assay_types"
    __table_args__ = {"schema": DB_SCHEMA}

    id: Optional[int] = Field(default=None, primary_key=True)
    name: str = Field(nullable=False)
    description: Optional[str] = None

    properties: List["Property"] = Relationship(back_populates="assay_types", link_model=AssayTypePropertyLink)

    assays = Relationship(sa_relationship="Assay", back_populates="assay_type")
    assay_type_details = Relationship(sa_relationship="AssayTypeDetails", back_populates="assay_type")
    property_requirements = Relationship(sa_relationship="AssayTypeProperty", back_populates="assay_type")

    class Config:
        orm_mode = True


class AssayTypeDetails(SQLModel, table=True):
    __tablename__ = "assay_type_details"
    __table_args__ = {"schema": DB_SCHEMA}

    assay_type_id: int = Field(foreign_key=f"{DB_SCHEMA}.assay_types.id", primary_key=True)
    property_id: int = Field(foreign_key=f"{DB_SCHEMA}.properties.id", primary_key=True)
    required: bool = Field(default=False, nullable=False)

    value_datetime: Optional[datetime] = Field(default=None, sa_column_kwargs={"server_default": func.now()})
    value_uuid: Optional[UUID] = None
    value_num: Optional[float] = None
    value_string: Optional[str] = None

    assay_type = Relationship(sa_relationship="AssayType", back_populates="assay_type_details")
    property = Relationship(sa_relationship="Property")


class AssayDetail(SQLModel, table=True):
    __tablename__ = "assay_details"
    __table_args__ = {"schema": DB_SCHEMA}

    assay_id: int = Field(foreign_key=f"{DB_SCHEMA}.assays.id", primary_key=True)
    property_id: int = Field(foreign_key=f"{DB_SCHEMA}.properties.id", primary_key=True)

    value_datetime: Optional[datetime] = Field(default=None, sa_column_kwargs={"server_default": func.now()})
    value_uuid: Optional[UUID] = None
    value_num: Optional[float] = None
    value_string: Optional[str] = None

    class Config:
        orm_mode = True

    assay = Relationship(sa_relationship="Assay", back_populates="assay_details")
    property = Relationship(sa_relationship="Property")


class AssayTypeProperty(SQLModel, table=True):
    __tablename__ = "assay_type_properties"
    __table_args__ = {"schema": DB_SCHEMA}

    assay_type_id: int = Field(foreign_key=f"{DB_SCHEMA}.assay_types.id", primary_key=True)
    property_id: int = Field(foreign_key=f"{DB_SCHEMA}.properties.id", primary_key=True)
    required: bool = Field(default=False, nullable=False)

    assay_type = Relationship(sa_relationship="AssayType", back_populates="property_requirements")
    property = Relationship(sa_relationship="Property")


class AssayResult(SQLModel, table=True):
    __tablename__ = "assay_results"
    __table_args__ = {"schema": DB_SCHEMA}

    id: Optional[int] = Field(default=None, primary_key=True)
    batch_id: int = Field(foreign_key=f"{DB_SCHEMA}.batches.id", nullable=False)
    assay_id: int = Field(foreign_key=f"{DB_SCHEMA}.assays.id", nullable=False)
    property_id: int = Field(foreign_key=f"{DB_SCHEMA}.properties.id", nullable=False)

    value_qualifier: int = Field(default=0, nullable=False)
    value_num: Optional[float] = None
    value_string: Optional[str] = None
    value_bool: Optional[bool] = None

    batch = Relationship(sa_relationship="Batch", back_populates="assay_results")
    assay = Relationship(sa_relationship="Assay", back_populates="assay_results")
    property = Relationship(sa_relationship="Property", back_populates="assay_results")
