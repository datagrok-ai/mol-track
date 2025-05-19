from sqlalchemy import Boolean, Column, ForeignKey, Integer, String, Text, Float, DateTime, Date, Enum, CheckConstraint, Table, UniqueConstraint
from sqlalchemy.orm import relationship
from sqlalchemy.sql import func
import enum
import os
from datetime import datetime

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

class Compound(Base):
    __tablename__ = "compounds"
    __table_args__ = {"schema": DB_SCHEMA}

    id = Column(Integer, primary_key=True, index=True)
    canonical_smiles = Column(Text, nullable=False)
    original_molfile = Column(Text)  # as sketched by the chemist
    inchi = Column(Text, nullable=False)
    inchikey = Column(Text, nullable=False, unique=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now())
    is_archived = Column(Boolean, default=False)
    
    # Relationships
    batches = relationship("Batch", back_populates="compound")
    compound_synonyms = relationship("CompoundSynonym", back_populates="compound")

class Batch(Base):
    __tablename__ = "batches"
    __table_args__ = {"schema": DB_SCHEMA}

    id = Column(Integer, primary_key=True, index=True)
    compound_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.compounds.id"), nullable=False)
    batch_number = Column(Text, nullable=False)
    amount = Column(Float)
    amount_unit = Column(Text)
    purity = Column(Float)
    notes = Column(Text)
    expiry_date = Column(Date)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    created_by = Column(Integer)
    
    # Relationships
    compound = relationship("Compound", back_populates="batches")
    assay_results = relationship("AssayResult", back_populates="batch")
    batch_details = relationship("BatchDetail", back_populates="batch")
    batch_synonyms = relationship("BatchSynonym", back_populates="batch")

class BatchDetail(Base):
    __tablename__ = "batch_details"
    __table_args__ = {"schema": DB_SCHEMA}

    id = Column(Integer, primary_key=True, index=True)
    batch_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.batches.id"), nullable=False)
    property_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.properties.id"), nullable=False)
    
    value_qualifier = Column(Integer, nullable=False, default=0)  # 0 for "=", 1 for "<", 2 for ">"
    value_datetime = Column(DateTime(timezone=True))
    value_uuid = Column(String(36))
    value_num = Column(Float)
    value_string = Column(Text)
    
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
            name="check_value_type"
        ),
        CheckConstraint(
            "property_class IN ('CALCULATED', 'MEASURED', 'PREDICTED')",
            name="check_property_class"
        ),
        {"schema": DB_SCHEMA}
    )

    id = Column(Integer, primary_key=True, index=True)
    name = Column(Text, nullable=False)
    value_type = Column(Text)
    property_class = Column(Text)
    unit = Column(Text)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    
    # Relationship to assay results
    assay_results = relationship("AssayResult", back_populates="property")
    batch_details = relationship("BatchDetail", back_populates="property")

class AssayTypeDetail(Base):
    __tablename__ = "assay_type_details"
    __table_args__ = {"schema": DB_SCHEMA}

    assay_type_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.assay_types.id"), primary_key=True)
    property_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.properties.id"), primary_key=True)
    
    value_datetime = Column(DateTime(timezone=True))
    value_uuid = Column(String(36))
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
    value_uuid = Column(String(36))
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
    created_on = Column(DateTime, server_default=func.now())
    updated_on = Column(DateTime, server_default=func.now(), onupdate=func.now())
    
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


class SynonymType(Base):
    __tablename__ = "synonym_types"
    __table_args__ = {"schema": DB_SCHEMA}

    id = Column(Integer, primary_key=True, index=True)
    synonym_level = Column(Text)
    name = Column(Text)
    pattern = Column(Text)
    description = Column(Text)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now())
    created_by = Column(UUID(as_uuid=True), nullable=False)
    updated_by = Column(UUID(as_uuid=True), nullable=False)

    # Relationships
    compound_synonyms = relationship("CompoundSynonym", back_populates="synonym_type")
    batch_synonyms = relationship("BatchSynonym", back_populates="synonym_type")

class CompoundSynonym(Base):
    __tablename__ = "compound_synonyms"
    __table_args__ = {"schema": DB_SCHEMA}

    id = Column(Integer, primary_key=True, index=True)
    compound_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.compounds.id"), nullable=False)
    synonym_type_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.synonym_types.id"), nullable=False)
    synonym_value = Column(Text)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now())
    created_by = Column(UUID(as_uuid=True), nullable=False)
    updated_by = Column(UUID(as_uuid=True), nullable=False)

    # Relationships
    compound = relationship("Compound", back_populates="compound_synonyms")
    synonym_type = relationship("SynonymType", back_populates="compound_synonyms")

class BatchSynonym(Base):
    __tablename__ = "batch_synonyms"
    __table_args__ = {"schema": DB_SCHEMA}

    id = Column(Integer, primary_key=True, index=True)
    batch_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.batches.id"), nullable=False)
    synonym_type_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.synonym_types.id"), nullable=False)
    synonym_value = Column(Text)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now())
    created_by = Column(UUID(as_uuid=True), nullable=False)
    updated_by = Column(UUID(as_uuid=True), nullable=False)

    # Relationships
    batch = relationship("Batch", back_populates="batch_synonyms")
    synonym_type = relationship("SynonymType", back_populates="batch_synonyms")
