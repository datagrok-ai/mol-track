from sqlalchemy import UUID, Boolean, Column, ForeignKey, Integer, String, Text, Float, DateTime, Date, Enum, CheckConstraint, Table, UniqueConstraint
from sqlalchemy.orm import relationship
from sqlalchemy.sql import func
import enum
import os

# Handle both package imports and direct execution
try:
    # When imported as a package (for tests)
    from .database import Base
except ImportError:
    # When run directly
    from database import Base

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

class User(Base):
    __tablename__ = "users"
    __table_args__ = {"schema": DB_SCHEMA}

    id = Column(UUID(as_uuid=True), primary_key=True, nullable=False)
    email = Column(Text, nullable=False)
    first_name = Column(Text, nullable=False)
    last_name = Column(Text, nullable=False)
    has_password = Column(Boolean, nullable=False)
    is_active = Column(Boolean, default=False, nullable=False)
    is_service_account = Column(Boolean, default=False, nullable=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now())
    created_by = Column(UUID(as_uuid=True), nullable=False)
    updated_by = Column(UUID(as_uuid=True), nullable=False)

class Compound(Base):
    __tablename__ = "compounds"
    __table_args__ = {"schema": DB_SCHEMA}

    id = Column(Integer, primary_key=True, index=True)
    canonical_smiles = Column(Text, nullable=False)
    original_molfile = Column(Text, nullable=False)
    inchi = Column(Text, nullable=False)
    inchikey = Column(Text, nullable=False, unique=True)
    molregno = Column(Integer, nullable=False)
    formula = Column(Text, nullable=False)
    hash_mol = Column(UUID(as_uuid=True), nullable=False)
    hash_tautomer = Column(UUID(as_uuid=True), nullable=False)
    hash_canonical_smiles = Column(UUID(as_uuid=True), nullable=False)
    hash_no_stereo_smiles = Column(UUID(as_uuid=True), nullable=False)
    hash_no_stereo_tautomer = Column(UUID(as_uuid=True), nullable=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now())
    created_by = Column(UUID(as_uuid=True), nullable=False)
    updated_by = Column(UUID(as_uuid=True), nullable=False)
    deleted_at = Column(DateTime(timezone=True), server_default=func.now())
    deleted_by = Column(UUID(as_uuid=True), nullable=False)
    is_archived = Column(Boolean, default=False)
    
    # Relationships
    batches = relationship("Batch", back_populates="compound")

class Batch(Base):
    __tablename__ = "batches"
    __table_args__ = {"schema": DB_SCHEMA}

    id = Column(Integer, primary_key=True, index=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), server_default=func.now())
    created_by = Column(UUID(as_uuid=True), nullable=False)
    updated_by = Column(UUID(as_uuid=True), nullable=False)
    compound_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.compounds.id"), nullable=False)
    batch_regno = Column(Integer, nullable=False)
    notes = Column(Text)
    
    # Relationships
    compound = relationship("Compound", back_populates="batches")
    assay_results = relationship("AssayResult", back_populates="batch")
    batch_details = relationship("BatchDetail", back_populates="batch")

class BatchDetail(Base):
    __tablename__ = "batch_details"
    __table_args__ = {"schema": DB_SCHEMA}

    id = Column(Integer, primary_key=True, index=True)
    batch_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.batches.id"), nullable=False)
    property_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.properties.id"), nullable=False)
    
    value_qualifier = Column(Integer, nullable=False, default=0)  # 0 for "=", 1 for "<", 2 for ">"
    value_datetime = Column(DateTime(timezone=True))
    value_uuid = Column(UUID(as_uuid=True))
    value_num = Column(Float)
    value_string = Column(Text)

    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), server_default=func.now())
    created_by = Column(UUID(as_uuid=True), nullable=False)
    updated_by = Column(UUID(as_uuid=True), nullable=False)
    
    # Relationships
    batch = relationship("Batch", back_populates="batch_details")
    property = relationship("Property", back_populates="batch_details")

class ValueType(enum.Enum):
    INT = "int"
    DOUBLE = "double"
    BOOL = "bool"
    DATETIME = "datetime"
    STRING = "string"

class PropertyClass(enum.Enum):
    CALCULATED = "CALCULATED"
    MEASURED = "MEASURED"
    PREDICTED = "PREDICTED"

class Property(Base):
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

    id = Column(Integer, primary_key=True, index=True)
    name = Column(Text, nullable=False)
    value_type = Column(Text, nullable=False)
    property_class = Column(Text, nullable=False)
    unit = Column(Text)
    scope = Column(Text, nullable=False)
    semantic_type_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.semantic_types.id"))
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), server_default=func.now())
    created_by = Column(UUID(as_uuid=True), nullable=False)
    updated_by = Column(UUID(as_uuid=True), nullable=False)

    # Relationship to assay results
    assay_results = relationship("AssayResult", back_populates="property")
    batch_details = relationship("BatchDetail", back_populates="property")

class AssayTypeDetail(Base):
    __tablename__ = "assay_type_details"
    __table_args__ = {"schema": DB_SCHEMA}

    assay_type_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.assay_types.id"), primary_key=True)
    property_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.properties.id"), primary_key=True)
    
    value_datetime = Column(DateTime(timezone=True))
    value_uuid = Column(UUID(as_uuid=True))
    value_num = Column(Float)
    value_string = Column(Text)
    
    # Relationships
    assay_type = relationship("AssayType", back_populates="assay_type_details")
    property = relationship("Property")

class AssayDetail(Base):
    __tablename__ = "assay_details"
    __table_args__ = {"schema": DB_SCHEMA}

    assay_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.assays.id"), primary_key=True)
    property_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.properties.id"), primary_key=True)
    
    value_datetime = Column(DateTime(timezone=True))
    value_uuid = Column(UUID(as_uuid=True))
    value_num = Column(Float)
    value_string = Column(Text)
    
    # Relationships
    assay = relationship("Assay", back_populates="assay_details")
    property = relationship("Property")

# Redefining AssayTypeProperty which was previously just a Table
class AssayTypeProperty(Base):
    __tablename__ = "assay_type_properties"
    __table_args__ = {"schema": DB_SCHEMA}

    assay_type_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.assay_types.id"), primary_key=True)
    property_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.properties.id"), primary_key=True)
    required = Column(Boolean, default=False)
    
    # Relationships
    assay_type = relationship("AssayType", back_populates="property_requirements")
    property = relationship("Property")

class AssayType(Base):
    __tablename__ = "assay_types"
    __table_args__ = {"schema": DB_SCHEMA}

    id = Column(Integer, primary_key=True, index=True)
    name = Column(Text, nullable=False)
    description = Column(Text)
    created_at = Column(DateTime, server_default=func.now())
    updated_at = Column(DateTime, server_default=func.now(), onupdate=func.now())
    created_by = Column(UUID(as_uuid=True), nullable=False)
    updated_by = Column(UUID(as_uuid=True), nullable=False)
    
    # Relationships - use secondary join for many-to-many
    properties = relationship(
        "Property", 
        secondary=f"{DB_SCHEMA}.assay_type_properties",
        primaryjoin=f"AssayType.id==AssayTypeProperty.assay_type_id",
        secondaryjoin=f"AssayTypeProperty.property_id==Property.id",
        backref="assay_types"
    )
    assays = relationship("Assay", back_populates="assay_type")
    assay_type_details = relationship("AssayTypeDetail", back_populates="assay_type")
    property_requirements = relationship("AssayTypeProperty", back_populates="assay_type")

class Assay(Base):
    __tablename__ = "assays"
    __table_args__ = {"schema": DB_SCHEMA}

    id = Column(Integer, primary_key=True, index=True)
    assay_type_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.assay_types.id"), nullable=False)
    name = Column(Text, nullable=False)
    description = Column(Text)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), server_default=func.now())
    created_by = Column(UUID(as_uuid=True), nullable=False)
    updated_by = Column(UUID(as_uuid=True), nullable=False)
    
    # Relationships - use assay_type_properties via assay_type to get list of expected properties
    # No direct properties relationship as assay_properties table no longer exists
    assay_type = relationship("AssayType", back_populates="assays")
    assay_results = relationship("AssayResult", back_populates="assay")
    assay_details = relationship("AssayDetail", back_populates="assay")

class AssayResult(Base):
    __tablename__ = "assay_results"
    __table_args__ = {"schema": DB_SCHEMA}

    id = Column(Integer, primary_key=True, index=True)
    batch_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.batches.id"), nullable=False)
    assay_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.assays.id"), nullable=False)
    property_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.properties.id"), nullable=False)
    
    # Value fields
    value_qualifier = Column(Integer, nullable=False, default=0)  # 0 for "=", 1 for "<", 2 for ">"
    value_num = Column(Float)
    value_string = Column(Text)
    value_bool = Column(Boolean)
    
    # Relationships
    batch = relationship("Batch", back_populates="assay_results")
    assay = relationship("Assay", back_populates="assay_results")
    property = relationship("Property", back_populates="assay_results")

class SemanticType(Base):
    __tablename__ = 'semantic_types'
    __table_args__ = {"schema": DB_SCHEMA}

    id = Column(Integer, primary_key=True, nullable=False)
    name = Column(Text, nullable=False)
    description = Column(Text, nullable=True)