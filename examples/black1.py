# example entry for Black data set of hepatocyte stability
black_assay_setup = {
    "assay.type.name": "Hepatocyte Stability",
    "assay.type.form": "in vitro",
    "assay.biological.system": "cell",
    "assay.cell.species": "Human",
    "assay.properties": [
        { "name": "Reported CLint", "unit": "uL/min/106 cells"},
        { "name": "Mean HTC recovery", "unit": "%"},
        { "name": "SD HTC recovery", "unit": "%"},
        { "name": "Dosed Concentration", "unit": "uM"},
        { "name": "Cell Lot"},
    ],
}

black_compound_setup = {
    "property_types": [],
    "synonym_types": [
        {"name": "batch_corporate_id", "level": "batch", "pattern": ""},
        {"name": "corporate_id", "level": "compound", "pattern": ""},
        {"name": "CAS number", "level": "compound", "pattern": r"\b[1-9]{1}[0-9]{1,6}-\d{2}-\d\b"},
        {"name": "common name", "level": "compound", "pattern": ""},
        {"name": "USAN", "level": "compound", "pattern": ""}
    ]
}
black1 = {
    "batch":
        {
            "batch_corporate_id": "EPA-001-001",
            "compound": {
                "smiles": "OC(=O)c1cc2c(cc1)c(=O)O",
                "common_name": "1,3-benezenedicarboxylic acid",
                "CAS Number": "121-91-5"
            }
        },
    "assay_detail": {
        "assay.type.name": "Hepatocyte Stability",
        "assay.experiment_date": "10-Jan-2024",
        "assay.cell_species": "Human",
        "assay.cell_lot": "H2",
    },
    "assay.results": [
        { 
            "batch_name": "EPA-001-001",
            "Dosed Concentration": 0.1,  # experimental condition
            "Report CLint": None,
            "Mean HTC recovery": 59.51,
            "SD HTC recovery": 84.16
        },
        { 
            "batch_name": "EPA-001-001", 
            "Dosed Concentration": 1,  # experimental condition
            "Report CLint": None,
            "Mean HTC recovery": 69.57,
            "SD HTC recovery": 1.97,
        }
    ]
}

celecoxib = {
    "batch":
        {
            "batch_corporate_id": "EPA-120-001",
            "compound": {
                "smiles": "c1cc(C)ccc1c2cc(C(F)(F)F)nn2c3ccc(cc3)S(=O)(=O)N",
                "USAN": "celecoxib",
                "CAS Number": "169590-42-5"
            }
        },
    "assay_detail": {
        "assay.type.name": "Hepatocyte Stability",
        "assay.experiment_date": "11-Jan-2024",
        "assay.cell_species": "Human",
        "assay.cell_lot": "H2",
    },
    "assay.results": [
        { 
            "batch_name": "EPA-120-001",
            "Dosed Concentration": 0.1,  # experimental condition
            "Report CLint": 24.55,
            "Mean HTC recovery": 105.42,
            "SD HTC recovery": 10.78
        },
        { 
            "batch_name": "EPA-120-001", 
            "Dosed Concentration": 1,  # experimental condition
            "Report CLint": 22.93,
            "Mean HTC recovery": 100.15,
            "SD HTC recovery": 2.49,
        }
    ]
}