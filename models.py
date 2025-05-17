from typing import List, Optional, Union
from sqlalchemy import UUID, Boolean, Column, ForeignKey, Integer, String, Text, Float, DateTime, Date, Enum, CheckConstraint, Table, UniqueConstraint
from sqlmodel import SQLModel, Field, func, Relationship
from sqlalchemy.orm import relationship
from sqlalchemy.sql import func
import enum
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

# This association table is no longer needed since assay_properties isn't in the schema
# Association table for Assay-Property many-to-many relationship
# AssayProperty = Table(
#     "assay_properties",
#     Base.metadata,
#     Column("assay_id", Integer, ForeignKey(f"{DB_SCHEMA}.assays.id"), primary_key=True),
#     Column("property_id", Integer, ForeignKey(f"{DB_SCHEMA}.properties.id"), primary_key=True),
#     schema=DB_SCHEMA
# )

# No longer needed as we've converted to an ORM class
# Association table for AssayType-Property many-to-many relationship
# AssayTypeProperty = Table(

class User(SQLModel, table=True):
    __tablename__ = "users"
    __table_args__ = {"schema": DB_SCHEMA}

    id: uuid.UUID = Field(sa_column = Column(UUID(as_uuid=True), primary_key=True, nullable=False))
    email: str = Field(sa_column=Column(Text, nullable=False))
    first_name: str = Field(sa_column = Column(Text, nullable=False))
    last_name: str = Field(sa_column = Column(Text, nullable=False))
    has_password: bool = Field(sa_column=Column(Boolean, nullable=False))
    is_active: bool = Field(sa_column=Column(Boolean, default=False, nullable=False))
    is_service_account: bool = Field(sa_column=Column(Boolean, default=False, nullable=False))
    created_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now(), nullable=False))
    updated_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False))
    created_by: uuid.UUID = Field(sa_column=Column(UUID(as_uuid=True), nullable=False))
    updated_by: uuid.UUID = Field(sa_column=Column(UUID(as_uuid=True), nullable=False))

    model_config = {
        "arbitrary_types_allowed": True
    }

class CompoundBase(SQLModel):
    canonical_smiles: str = Field(sa_column=Column(Text, nullable=False))
    original_molfile: Optional[str] = Field(sa_column=Column(Text))  # as sketched by the chemist
    inchi: str = Field(sa_column=Column(Text, nullable=False))
    inchikey: str = Field(sa_column=Column(Text, nullable=False, unique=True))
    is_archived: Optional[bool] = Field(sa_column=Column(Boolean, default=False))

    # @validator("inchi", "inchikey", always=True)
    # def set_inchi(cls, v, values):
    #     if v is not None:
    #         return v
    #     # In a real implementation, these would be calculated
    #     # based on canonical_smiles
    #     if "canonical_smiles" in values and values["canonical_smiles"] is not None:
    #         # Placeholder logic
    #         if "inchi" in values:
    #             return "InChI=1S/" + values["canonical_smiles"]
    #         else:
    #             return "INCHIKEY" + values["canonical_smiles"]
    #     return v

class CompoundCreate(CompoundBase):
    smiles: str

class CompoundPublic(CompoundBase):
    molregno: int = Field(sa_column=Column(Integer, nullable=False))
    formula: str = Field(sa_column=Column(Text, nullable=False))
    hash_mol: uuid.UUID  = Field(sa_column=Column(UUID(as_uuid=True), nullable=False))
    hash_tautomer: uuid.UUID = Field(sa_column=Column(UUID(as_uuid=True), nullable=False))
    hash_canonical_smiles: uuid.UUID = Field(sa_column=Column(UUID(as_uuid=True), nullable=False))
    hash_no_stereo_smiles: uuid.UUID = Field(sa_column=Column(UUID(as_uuid=True), nullable=False))
    hash_no_stereo_tautomer: uuid.UUID = Field(sa_column=Column(UUID(as_uuid=True), nullable=False))
    created_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now()))
    updated_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now()))

class CompoundPublicWithBatches(CompoundPublic):
    pass
    # batches: List[Batch] = []

class Compound(CompoundPublic, table=True):
    __tablename__ = "compounds"
    __table_args__ = {"schema": DB_SCHEMA}

    id: int = Field(sa_column=Column(Integer, primary_key=True, index=True))
    created_by: uuid.UUID = Field(sa_column=Column(UUID(as_uuid=True), nullable=False))
    updated_by: uuid.UUID = Field(sa_column=Column(UUID(as_uuid=True), nullable=False))
    deleted_at: Optional[datetime] = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now()))
    deleted_by: Optional[uuid.UUID] = Field(sa_column=Column(UUID(as_uuid=True), nullable=False))
    
    # Relationships
    batches: List["Batch"] = Relationship(sa_relationship="Batch", back_populates="compound")

    class Config:
        from_attributes = True
        orm_mode=True

class BatchBase(SQLModel):
    notes: Optional[str] = Field(sa_column=Column(Text))
    compound_id: int = Field(sa_column=Column(Integer, ForeignKey(f"{DB_SCHEMA}.compounds.id"), nullable=False))

class BatchPublic(BatchBase):
    id: int = Field(sa_column=Column(Integer, primary_key=True, index=True))
    created_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now(), nullable=False))

class BatchPublicWithDetails(BatchPublic):
    pass
    # batch_details: List["BatchDetail"] = []

    # class Config:
    #     from_attributes = True

class Batch(BatchPublic, table=True):
    __tablename__ = "batches"
    __table_args__ = {"schema": DB_SCHEMA}

    updated_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now(), nullable=False))
    created_by: uuid.UUID = Field(sa_column=Column(UUID(as_uuid=True), nullable=False))
    updated_by: uuid.UUID = Field(sa_column=Column(UUID(as_uuid=True), nullable=False))
    batch_regno: int = Field(sa_column=Column(Integer, nullable=False))
    
    # Relationships
    compound = Relationship(sa_relationship="Compound", back_populates="batches")
    assay_results = Relationship(sa_relationship="AssayResult", back_populates="batch")
    batch_details = Relationship(sa_relationship="BatchDetail", back_populates="batch")

    class Config:
        from_attributes = True
        orm_mode=True

class BatchDetailBase(SQLModel):
    batch_id: int = Field(sa_column=Column(Integer, ForeignKey(f"{DB_SCHEMA}.batches.id"), nullable=False))
    property_id: int = Field(sa_column=Column(Integer, ForeignKey(f"{DB_SCHEMA}.properties.id"), nullable=False))
    
    value_qualifier: int = Field(sa_column=Column(Integer, nullable=False, default=0))  # 0 for "=", 1 for "<", 2 for ">"
    value_datetime: Optional[datetime] = Field(sa_column=Column(DateTime(timezone=True)))
    value_uuid: Optional[uuid.UUID] = Field(sa_column=Column(UUID(as_uuid=True)))
    value_num: Optional[float] = Field(sa_column=Column(Float))
    value_string: Optional[str] = Field(sa_column=Column(Text))

class BatchDetailPublic(BatchDetailBase):
    id: int = Field(sa_column=Column(Integer, primary_key=True, index=True))

class BatchDetail(BatchDetailPublic, table=True):
    __tablename__ = "batch_details"
    __table_args__ = {"schema": DB_SCHEMA}

    created_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now(), nullable=False))
    updated_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now(), nullable=False))
    created_by: uuid.UUID = Field(sa_column=Column(UUID(as_uuid=True), nullable=False))
    updated_by: uuid.UUID = Field(sa_column=Column(UUID(as_uuid=True), nullable=False))

    # Relationships
    batch = Relationship(sa_relationship="Batch", back_populates="batch_details")
    property = Relationship(sa_relationship="Property", back_populates="batch_details")


class ValueType(str, enum.Enum):
    int = "int"
    double = "double"
    bool = "bool"
    datetime = "datetime"
    string = "string"


class PropertyClass(str, enum.Enum):
    CALCULATED = "CALCULATED"
    MEASURED = "MEASURED"
    PREDICTED = "PREDICTED"

class ScopeClass(str, enum.Enum):
    BATCH  = "BATCH"
    COMPOUND = "COMPOUND"
    ASSAY = "ASSAY"
    SYSTEM = "SYSTEM" 

class PropertyBase(SQLModel):
    name: str = Field(sa_column=Column(Text, nullable=False))
    value_type: ValueType = Field(sa_column=Column(Enum(ValueType), nullable=False))
    property_class: PropertyClass = Field(sa_column=Column(Enum(PropertyClass), nullable=False))
    unit: Optional[str] = Field(sa_column=Column(Text))
    scope: ScopeClass = Field(sa_column=Column(Enum(ScopeClass), nullable=False)) #define the scope
    semantic_type_id: int = Field(sa_column=Column(Integer, ForeignKey(f"{DB_SCHEMA}.semantic_types.id")))

class PropertyResponse(PropertyBase):
    id: int = Field(sa_column=Column(Integer, primary_key=True, index=True))
    created_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now()))

# Redefining AssayTypeProperty which was previously just a Table
class AssayTypeProperty(SQLModel, table=True):
    __tablename__ = "assay_type_properties"
    __table_args__ = {"schema": DB_SCHEMA}

    assay_type_id: int = Field(sa_column=Column(Integer, ForeignKey(f"{DB_SCHEMA}.assay_types.id"), primary_key=True))
    property_id: int = Field(sa_column=Column(Integer, ForeignKey(f"{DB_SCHEMA}.properties.id"), primary_key=True))
    required: bool = Field(sa_column=Column(Boolean, default=False))
    
    # Relationships
    assay_type = Relationship(sa_relationship="AssayType", back_populates="property_requirements")
    property = Relationship(sa_relationship="Property")

class Property(PropertyResponse, table=True):
    __tablename__ = "properties"
    __table_args__ = (
        CheckConstraint(
            "value_type IN ('int', 'double', 'bool', 'datetime', 'string')",
            name="properties_value_type_check"
        ),
        CheckConstraint(
            "property_class IN ('CALCULATED', 'MEASURED', 'PREDICTED')",
            name="properties_property_class_check"
        ),
        CheckConstraint(
            "scope IN ('BATCH', 'COMPOUND', 'ASSAY', 'SYSTEM')",
            name="properties_scope_check"
        ),
        {"schema": DB_SCHEMA}
    )
    
    updated_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now()))
    created_by: uuid.UUID = Field(sa_column=Column(UUID(as_uuid=True), nullable=False))
    updated_by: uuid.UUID = Field(sa_column=Column(UUID(as_uuid=True), nullable=False))

    # Relationship to assay results
    assay_results = Relationship(sa_relationship="AssayResult", back_populates="property")
    batch_details = Relationship(sa_relationship="BatchDetail", back_populates="property")
    assay_types: List["AssayType"] = Relationship(back_populates="properties", link_model=AssayTypeProperty)


class AssayTypeDetail(SQLModel, table=True):
    __tablename__ = "assay_type_details"
    __table_args__ = {"schema": DB_SCHEMA}

    assay_type_id: int = Field(sa_column=Column(Integer, ForeignKey(f"{DB_SCHEMA}.assay_types.id"), primary_key=True))
    property_id: int = Field(sa_column=Column(Integer, ForeignKey(f"{DB_SCHEMA}.properties.id"), primary_key=True))
    
    value_datetime: Optional[datetime] = Field(sa_column=Column(DateTime(timezone=True)))
    value_uuid: Optional[uuid.UUID] = Field(sa_column=Column(UUID(as_uuid=True)))
    value_num: Optional[float] = Field(sa_column=Column(Float))
    value_string: Optional[str] = Field(sa_column=Column(Text))

    # Relationships
    assay_type = Relationship(sa_relationship="AssayType", back_populates="assay_type_details")
    property = Relationship(sa_relationship="Property")

class AssayDetail(SQLModel, table=True):
    __tablename__ = "assay_details"
    __table_args__ = {"schema": DB_SCHEMA}

    assay_id: int = Field(sa_column=Column(Integer, ForeignKey(f"{DB_SCHEMA}.assays.id"), primary_key=True))
    property_id: int = Field(sa_column=Column(Integer, ForeignKey(f"{DB_SCHEMA}.properties.id"), primary_key=True))
    
    value_datetime: Optional[datetime] = Field(sa_column=Column(DateTime(timezone=True)))
    value_uuid: Optional[uuid.UUID] = Field(sa_column=Column(UUID(as_uuid=True)))
    value_num: Optional[float] = Field(sa_column=Column(Float))
    value_string: Optional[str] = Field(sa_column=Column(Text))

    # Relationships
    assay = Relationship(sa_relationship="Assay", back_populates="assay_details")
    property = Relationship(sa_relationship="Property")

class AssayType(SQLModel, table=True):
    __tablename__ = "assay_types"
    __table_args__ = {"schema": DB_SCHEMA}

    id: int = Field(sa_column=Column(Integer, primary_key=True, index=True))
    name: str = Field(sa_column=Column(Text, nullable=False))
    description: Optional[str] = Field(sa_column=Column(Text))
    created_at: datetime = Field(sa_column=Column(DateTime, server_default=func.now()))
    updated_at: datetime = Field(sa_column=Column(DateTime, server_default=func.now(), onupdate=func.now()))
    created_by: uuid.UUID = Field(sa_column=Column(UUID(as_uuid=True), nullable=False))
    updated_by: uuid.UUID = Field(sa_column=Column(UUID(as_uuid=True), nullable=False))
    
    # Relationships - use secondary join for many-to-many
    properties: List["Property"] = Relationship(
        back_populates="assay_types",
        link_model=AssayTypeProperty
    )

    assays = Relationship(sa_relationship="Assay", back_populates="assay_type")
    assay_type_details = Relationship(sa_relationship="AssayTypeDetail", back_populates="assay_type")
    property_requirements = Relationship(sa_relationship="AssayTypeProperty", back_populates="assay_type")

class AssayBase(SQLModel):
    name: str = Field(sa_column=Column(Text, nullable=False))
    description: Optional[str] = Field(sa_column=Column(Text))
    assay_type_id: int = Field(sa_column=Column(Integer, ForeignKey(f"{DB_SCHEMA}.assay_types.id"), nullable=False))

class AssayPublic(AssayBase):
    id: int = Field(sa_column=Column(Integer, primary_key=True, index=True))
    created_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now()))

class AssayPublicWithProperties(AssayPublic):
    assay_type: Optional["AssayType"] = None
    assay_details: List["AssayDetail"] = []
    properties: List["Property"] = []
    # pass

class AssayCreate(AssayBase):
    property_ids: List[int] = []  # List of property IDs to associate with this assay

class Assay(AssayPublic, table=True):
    __tablename__ = "assays"
    __table_args__ = {"schema": DB_SCHEMA}

    updated_at: datetime = Field(sa_column=Column(DateTime(timezone=True), server_default=func.now()))
    created_by: uuid.UUID = Field(sa_column=Column(UUID(as_uuid=True), nullable=False))
    updated_by: uuid.UUID = Field(sa_column=Column(UUID(as_uuid=True), nullable=False))
    
    # Relationships - use assay_type_properties via assay_type to get list of expected properties
    # No direct properties relationship as assay_properties table no longer exists
    assay_type = Relationship(sa_relationship="AssayType", back_populates="assays")
    assay_results = Relationship(sa_relationship="AssayResult", back_populates="assay")
    assay_details = Relationship(sa_relationship="AssayDetail", back_populates="assay")

# AssayResult schemas
class AssayResultBase(SQLModel):
    batch_id: int = Field(sa_column=Column(Integer, ForeignKey(f"{DB_SCHEMA}.batches.id"), nullable=False))
    assay_id: int = Field(sa_column=Column(Integer, ForeignKey(f"{DB_SCHEMA}.assays.id"), nullable=False))
    property_id: int = Field(sa_column=Column(Integer, ForeignKey(f"{DB_SCHEMA}.properties.id"), nullable=False))
    
    # Value fields
    value_qualifier: int = Field(sa_column=Column(Integer, nullable=False, default=0))  # 0 for "=", 1 for "<", 2 for ">"
    value_num: Optional[float] = Field(sa_column=Column(Float))
    value_string: Optional[str] = Field(sa_column=Column(Text))
    value_bool: Optional[bool] = Field(sa_column=Column(Boolean))

class AssayResultPublic(AssayResultBase):
    id: int = Field(sa_column=Column(Integer, primary_key=True, index=True))

# Extended response model with backward compatibility field
class AssayResultResponse(AssayResultPublic):
    # Add a computed field for backward compatibility
    result_value: Optional[Union[float, str, bool]] = None

    # @validator("result_value", always=True)
    # def compute_result_value(cls, v, values):
    #     """Compute result_value from the appropriate typed value field"""
    #     if "value_num" in values and values["value_num"] is not None:
    #         return values["value_num"]
    #     elif "value_string" in values and values["value_string"] is not None:
    #         return values["value_string"]
    #     elif "value_bool" in values and values["value_bool"] is not None:
    #         return values["value_bool"]
    #     return None

class AssayResult(AssayResultPublic, table=True):
    __tablename__ = "assay_results"
    __table_args__ = {"schema": DB_SCHEMA}
    
    # Relationships
    batch = Relationship(sa_relationship="Batch", back_populates="assay_results")
    assay = Relationship(sa_relationship="Assay", back_populates="assay_results")
    property = Relationship(sa_relationship="Property", back_populates="assay_results")

class SemanticTypeBase(SQLModel):
    name: str = Field(sa_column=Column(Text, nullable=False))
    description: Optional[str] = Field(sa_column=Column(Text, nullable=True))

class SemanticType(SemanticTypeBase, table=True):
    __tablename__ = 'semantic_types'
    __table_args__ = {"schema": DB_SCHEMA}

    id: int = Field(sa_column=Column(Integer, primary_key=True, nullable=False))