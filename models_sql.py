from datetime import date, datetime
from typing import Any, Dict, List, Optional
from uuid import UUID
from sqlmodel import SQLModel, Field, Relationship
from pydantic import BaseModel


class CreatedUpdatedDetails(SQLModel, table=False):
    created_at: Optional[datetime] = Field(default_factory=datetime.utcnow)
    updated_at: Optional[datetime] = Field(default_factory=datetime.utcnow)
    created_by: UUID = Field(foreign_key="users.id", nullable=True)
    updated_by: UUID = Field(foreign_key="users.id", nullable=True)


class Base(CreatedUpdatedDetails, table=False):
    id: int = Field(default=None, primary_key=True)


class SoftDeleteTemplate(SQLModel, table=False):
    is_archived: bool = Field(default=False)
    deleted_at: Optional[datetime] = Field(default_factory=datetime.utcnow)
    deleted_by: UUID = Field(foreign_key="users.id")


class User(CreatedUpdatedDetails, table=True):
    __tablename__ = "users"

    id: UUID = Field(primary_key=True)
    email: str = Field(nullable=False, unique=True)
    first_name: str = Field(nullable=False)
    last_name: str = Field(nullable=False)
    has_password: bool = Field(nullable=False)
    is_active: bool = Field(default=False, nullable=False)
    is_service_account: bool = Field(default=False, nullable=False)


class Settings(SQLModel, table=True):
    __tablename__ = "settings"

    id: Optional[int] = Field(default=None, primary_key=True)
    name: str = Field(nullable=False)
    value: str = Field(nullable=False)
    description: str = Field(nullable=False)


class Compound(Base, SoftDeleteTemplate, table=True):
    __tablename__ = "compounds"

    canonical_smiles: str = Field(nullable=False)
    original_molfile: Optional[str] = Field(default=None)
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


class CompoundCreate(Compound, table=False):
    smiles: str
    is_archived: bool = False


class Batch(Base, table=True):
    __tablename__ = "batches"

    compound_id: int = Field(foreign_key="compounds.id", nullable=False)
    batch_number: str = Field(nullable=False)

    amount: Optional[float] = Field(default=None)
    amount_unit: Optional[str] = Field(default=None)
    purity: Optional[float] = Field(default=None)
    notes: Optional[str] = Field(default=None)
    expiry_date: Optional[date] = Field(default=None)

    # Relationships
    compound: "Compound" = Relationship(sa_relationship="Compound", back_populates="batches")
    assay_results: list["AssayResults"] = Relationship(sa_relationship="AssayResult", back_populates="batch")
    batch_details: list["BatchDetails"] = Relationship(sa_relationship="BatchDetail", back_populates="batch")


class Addition(Base, SoftDeleteTemplate, table=True):
    __tablename__ = "additions"

    name: str = Field(nullable=False)
    description: Optional[str] = Field(default=None)
    display_name: Optional[str] = Field(default=None)
    display_order: Optional[int] = Field(default=None)
    is_active: bool = Field(default=True, nullable=False)

    formula: Optional[str] = Field(default=None)
    molecular_weight: Optional[float] = Field(default=None)
    smiles: Optional[str] = Field(default=None)
    molfile: Optional[str] = Field(default=None)

    role: Optional[str] = Field(default=None)


class BatchAddition(CreatedUpdatedDetails, table=True):
    __tablename__ = "batch_additions"

    batch_id: int = Field(foreign_key="batches.id", primary_key=True)
    addition_id: int = Field(foreign_key="additions.id", primary_key=True)
    addition_equivalent: float = Field(default=1.0, nullable=False)

    # Relationships
    # batch: Optional[Batch] = Relationship(back_populates="batch_additions")
    # addition: Optional[Addition] = Relationship(back_populates="batch_additions")


class AssayPropertyLink(SQLModel, table=True):
    __tablename__ = "assay_properties"
    assay_id: int = Field(foreign_key="assays.id", primary_key=True)
    property_id: int = Field(foreign_key="properties.id", primary_key=True)


class AssayTypePropertyLink(SQLModel, table=True):
    assay_type_id: int = Field(foreign_key="assay_types.id", primary_key=True)
    property_id: int = Field(foreign_key="properties.id", primary_key=True)


class SemanticType(SQLModel, table=True):
    __tablename__ = "semantic_types"

    id: Optional[int] = Field(default=None, primary_key=True)
    name: str = Field(nullable=False)
    description: Optional[str] = Field(default=None)


class Property(Base, table=True):
    __tablename__ = "properties"

    name: str = Field(nullable=False)
    value_type: Optional[str] = Field(default=None)
    semantic_type_id: Optional[int] = Field(default=None, foreign_key="semantic_types.id")
    property_class: Optional[str] = Field(default=None)
    unit: Optional[str] = Field(default=None)

    # Relationships
    assay_types: List["AssayType"] = Relationship(back_populates="properties", link_model=AssayTypePropertyLink)
    assays: List["Assay"] = Relationship(
        sa_relationship="Assay", back_populates="properties", link_model=AssayPropertyLink
    )
    assay_results = Relationship(sa_relationship="AssayResult", back_populates="properties")
    batch_details = Relationship(sa_relationship="BatchDetail", back_populates="properties")
    # assay_results = relationship("AssayResult", back_populates="property")
    # batch_details = relationship("BatchDetail", back_populates="property")

    class Config:
        orm_mode = True


class BatchDetails(SQLModel, table=True):
    __tablename__ = "batch_details"

    id: Optional[int] = Field(default=None, primary_key=True)
    batch_id: int = Field(foreign_key="batches.id", nullable=False)
    property_id: int = Field(foreign_key="properties.id", nullable=False)

    value_datetime: Optional[datetime] = Field(default_factory=datetime.utcnow)
    value_uuid: Optional[str] = None
    value_num: Optional[float] = None
    value_string: Optional[str] = None

    # Relationships
    batch: Optional[Batch] = Relationship(sa_relationship="Batch", back_populates="batch_details")
    property: Optional[Property] = Relationship(sa_relationship="Property", back_populates="batch_details")


class SynonymType(Base, table=True):
    __tablename__ = "synonym_types"

    synonym_level: str = Field(nullable=False)
    name: str = Field(nullable=False)
    pattern: str = Field(nullable=False)
    description: str = Field(nullable=False)


class CompoundSynonym(Base, table=True):
    __tablename__ = "compound_synonyms"

    compound_id: int = Field(foreign_key="compounds.id", nullable=False)
    synonym_type_id: int = Field(foreign_key="synonym_types.id", nullable=False)
    synonym_value: str = Field(nullable=False)


class BatchSynonym(Base, table=True):
    __tablename__ = "batch_synonyms"

    batch_id: int = Field(foreign_key="batches.id", nullable=False)
    synonym_type_id: int = Field(foreign_key="synonym_types.id", nullable=False)
    synonym_value: str = Field(nullable=False)


class CompoundDetails(CreatedUpdatedDetails, table=True):
    __tablename__ = "compound_details"

    id: Optional[int] = Field(default=None, primary_key=True)
    compound_id: int = Field(foreign_key="compounds.id", nullable=False)
    property_id: int = Field(foreign_key="properties.id", nullable=False)

    value_datetime: Optional[datetime] = Field(default_factory=datetime.utcnow)
    value_uuid: Optional[str] = None
    value_num: Optional[float] = None
    value_string: Optional[str] = None


class Assay(CreatedUpdatedDetails, table=True):
    __tablename__ = "assays"

    id: Optional[int] = Field(default=None, primary_key=True)
    assay_type_id: int = Field(foreign_key="assay_types.id", nullable=False)
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

    id: Optional[int] = Field(default=None, primary_key=True)
    name: str = Field(nullable=False)
    description: Optional[str] = Field(default=None)

    properties: List["Property"] = Relationship(back_populates="assay_types", link_model=AssayTypePropertyLink)

    assays = Relationship(sa_relationship="Assay", back_populates="assay_type")
    assay_type_details = Relationship(sa_relationship="AssayTypeDetails", back_populates="assay_type")
    property_requirements = Relationship(sa_relationship="AssayTypeProperty", back_populates="assay_type")

    class Config:
        orm_mode = True


class AssayTypeCreate(Base):
    name: str = ""
    description: Optional[str] = None
    property_ids: List[int] = []  # List of property IDs to associate with this assay type
    property_requirements: List[Dict[str, Any]] = []  # List of property requirements
    property_details: List[Dict[str, Any]] = []  # List of property metadata


class AssayTypeDetails(SQLModel, table=True):
    __tablename__ = "assay_type_details"

    assay_type_id: int = Field(foreign_key="assay_types.id", primary_key=True)
    property_id: int = Field(foreign_key="properties.id", primary_key=True)
    required: bool = Field(default=False, nullable=False)

    value_datetime: Optional[datetime] = None
    value_uuid: Optional[str] = None
    value_num: Optional[float] = None
    value_string: Optional[str] = None

    # assay_type = relationship("AssayType", back_populates="assay_type_details")
    # property = relationship("Property")


class AssayDetails(SQLModel, table=True):
    __tablename__ = "assay_details"

    assay_id: int = Field(foreign_key="assays.id", primary_key=True)
    property_id: int = Field(foreign_key="properties.id", primary_key=True)

    value_datetime: Optional[datetime] = Field(default_factory=datetime.utcnow)
    value_uuid: Optional[str] = None
    value_num: Optional[float] = None
    value_string: Optional[str] = None

    class Config:
        orm_mode = True

    # assay = relationship("Assay", back_populates="assay_details")
    # property = relationship("Property")


class AssayTypeProperty(SQLModel, table=True):
    __tablename__ = "assay_type_properties"

    assay_type_id: int = Field(foreign_key="assay_types.id", primary_key=True)
    property_id: int = Field(foreign_key="properties.id", primary_key=True)
    required: bool = Field(default=False, nullable=False)

    # assay_type = relationship("AssayType", back_populates="property_requirements")
    # property = relationship("Property")


class AssayResults(SQLModel, table=True):
    __tablename__ = "assay_results"

    id: Optional[int] = Field(default=None, primary_key=True)
    batch_id: int = Field(foreign_key="batches.id", nullable=False)
    assay_id: int = Field(foreign_key="assays.id", nullable=False)
    property_id: int = Field(foreign_key="properties.id", nullable=False)

    value_qualifier: int = Field(default=0, nullable=False)
    value_num: Optional[float] = None
    value_string: Optional[str] = None
    value_bool: Optional[bool] = None

    # batch = relationship("Batch", back_populates="assay_results")
    # assay = relationship("Assay", back_populates="assay_results")
    # property = relationship("Property", back_populates="assay_results")


class CompoundQueryParams(BaseModel):
    substructure: Optional[str] = None
    skip: int = 0
    limit: int = 100

    class Config:
        from_attributes = True
