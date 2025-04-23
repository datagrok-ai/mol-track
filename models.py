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

# Association table for Assay-Property many-to-many relationship
AssayProperty = Table(
    "assay_properties",
    Base.metadata,
    Column("assay_id", Integer, ForeignKey(f"{DB_SCHEMA}.assays.id"), primary_key=True),
    Column("property_id", Integer, ForeignKey(f"{DB_SCHEMA}.properties.id"), primary_key=True),
    schema=DB_SCHEMA
)

# Association table for AssayType-Property many-to-many relationship
AssayTypeProperty = Table(
    "assay_type_properties",
    Base.metadata,
    Column("assay_type_id", Integer, ForeignKey(f"{DB_SCHEMA}.assay_types.id"), primary_key=True),
    Column("property_id", Integer, ForeignKey(f"{DB_SCHEMA}.properties.id"), primary_key=True),
    schema=DB_SCHEMA
)

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

class BatchDetail(Base):
    __tablename__ = "batch_details"
    __table_args__ = {"schema": DB_SCHEMA}

    id = Column(Integer, primary_key=True, index=True)
    batch_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.batches.id"), nullable=False)
    property_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.properties.id"), nullable=False)
    result_value = Column(Float)
    
    # Relationships
    batch = relationship("Batch", back_populates="batch_details")
    property = relationship("Property", back_populates="batch_details")

class ValueType(enum.Enum):
    NUMBER = "NUMBER"
    STRING = "STRING"
    BOOL = "BOOL"

class PropertyClass(enum.Enum):
    CALCULATED = "CALCULATED"
    MEASURED = "MEASURED"
    PREDICTED = "PREDICTED"

class Property(Base):
    __tablename__ = "properties"
    __table_args__ = (
        CheckConstraint(
            "value_type IN ('NUMBER', 'STRING', 'BOOL')",
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
    required = Column(Boolean, default=False)
    value_type = Column(Text)
    property_class = Column(Text)
    unit = Column(Text)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    
    # Relationship to assay results
    assay_results = relationship("AssayResult", back_populates="property")
    batch_details = relationship("BatchDetail", back_populates="property")

class AssayType(Base):
    __tablename__ = "assay_types"
    __table_args__ = {"schema": DB_SCHEMA}

    id = Column(Integer, primary_key=True, index=True)
    name = Column(Text, nullable=False)
    description = Column(Text)
    created_on = Column(DateTime, server_default=func.now())
    updated_on = Column(DateTime, server_default=func.now(), onupdate=func.now())
    
    # Relationships
    properties = relationship("Property", secondary=AssayTypeProperty, backref="assay_types")
    assays = relationship("Assay", back_populates="assay_type")
    
class Assay(Base):
    __tablename__ = "assays"
    __table_args__ = {"schema": DB_SCHEMA}

    id = Column(Integer, primary_key=True, index=True)
    assay_type_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.assay_types.id"), nullable=False)
    name = Column(Text, nullable=False)
    description = Column(Text)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    
    # Relationships
    properties = relationship("Property", secondary=AssayProperty, backref="assays")
    assay_type = relationship("AssayType", back_populates="assays")
    assay_results = relationship("AssayResult", back_populates="assay")

class AssayResult(Base):
    __tablename__ = "assay_results"
    __table_args__ = {"schema": DB_SCHEMA}

    id = Column(Integer, primary_key=True, index=True)
    batch_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.batches.id"), nullable=False)
    assay_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.assays.id"), nullable=False)
    property_id = Column(Integer, ForeignKey(f"{DB_SCHEMA}.properties.id"), nullable=False)
    result_value = Column(Float)
    
    # Relationships
    batch = relationship("Batch", back_populates="assay_results")
    assay = relationship("Assay", back_populates="assay_results")
    property = relationship("Property", back_populates="assay_results")
    