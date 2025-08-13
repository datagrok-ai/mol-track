import pytest
import yaml
from rdkit import Chem

from app.utils.chemistry_utils import apply_standardizer_operation, standardize_mol

# Some test molecules where taken from https://github.com/greglandrum/RSC_OpenScience_Standardization_202104/blob/main/MolStandardize%20pieces.ipynb


# Charged and salted mol
@pytest.fixture
def test_molecule():
    smiles = "c1ccccc1C(=O)[O-].[Na+]"
    return Chem.MolFromSmiles(smiles)


# Helper fixture to create a temporary YAML file for operations
@pytest.fixture
def temp_standardizer_config(tmp_path):
    config_data = {
        "operations": [
            {"type": "Cleanup", "enable": True},
            {"type": "FragmentParent", "enable": True},
            {"type": "Uncharger", "enable": True},
        ]
    }
    # Write the YAML data to a temporary file
    config_file = tmp_path / "standardizer_config.yaml"
    with open(config_file, "w") as f:
        yaml.dump(config_data, f)
    return str(config_file)


# Test individual standardization operations
def test_apply_standardizer_operation_cleanup():
    smiles = "C[S+2]([O-])([O-])CC"
    test_mol = Chem.MolFromSmiles(smiles)
    cleaned_mol = apply_standardizer_operation(test_mol, "Cleanup")
    assert Chem.MolToSmiles(cleaned_mol) == "CCS(C)(=O)=O"


def test_apply_standadizer_operation_fragment_parent():
    smiles = "[Na+].CC(=O)O"
    test_mol = Chem.MolFromSmiles(smiles)
    cleaned_mol = apply_standardizer_operation(test_mol, "FragmentParent")
    assert Chem.MolToSmiles(cleaned_mol) == "CC(=O)O"


def test_apply_standardizer_operation_uncharger():
    smiles = "[O-]CC(C(=O)[O-])C[NH+](C)C"
    test_mol = Chem.MolFromSmiles(smiles)
    uncharged_mol = apply_standardizer_operation(test_mol, "Uncharger")
    assert (
        len([a for a in uncharged_mol.GetAtoms() if a.GetFormalCharge() != 0]) == 0
    )  # No formal charges should remain
    assert Chem.MolToSmiles(uncharged_mol) == "CN(C)CC(CO)C(=O)O"


def test_apply_standardizer_operation_invalid():
    smiles = "CC(O)C(=O)O"
    test_mol = Chem.MolFromSmiles(smiles)
    with pytest.raises(ValueError, match="Unknown operation type: InvalidOperation"):
        apply_standardizer_operation(test_mol, "InvalidOperation")


def test_standardize_mol(client, test_db, test_molecule, temp_standardizer_config):
    client.patch(
        "/v1/admin/update-standardization-config",
        files={"file": ("standardize.yaml", open(temp_standardizer_config).read(), "application/x-yaml")},
    )
    standardized_mol = standardize_mol(test_molecule, test_db)
    assert Chem.MolToSmiles(standardized_mol) == "O=C(O)c1ccccc1"


def test_standardize_mol_missing_operation(client, test_db, test_molecule, tmp_path):
    config_data = {"operations": [{"enable": True}]}
    config_file = tmp_path / "invalid_config.yaml"
    with open(config_file, "w") as f:
        yaml.dump(config_data, f)

    client.patch(
        "/v1/admin/update-standardization-config",
        files={"file": ("standardize.yaml", open(config_file).read(), "application/x-yaml")},
    )

    with pytest.raises(ValueError, match="Operation type is missing in the configuration."):
        standardize_mol(test_molecule, test_db)


# Test for disabled operations. Ensure molecule remains unchanged and with salt since all operations are disabled
def test_standardize_mol_disabled_operations(client, test_db, test_molecule, tmp_path):
    config_data = {
        "operations": [
            {"type": "Cleanup", "enable": False},
            {"type": "FragmentParent", "enable": False},
            {"type": "Uncharger", "enable": False},
        ]
    }
    config_file = tmp_path / "disabled_config.yaml"
    with open(config_file, "w") as f:
        yaml.dump(config_data, f)

    client.patch(
        "/v1/admin/update-standardization-config",
        files={"file": ("standardize.yaml", open(config_file).read(), "application/x-yaml")},
    )

    original_smiles = Chem.MolToSmiles(test_molecule)
    standardized_mol = standardize_mol(test_molecule, test_db)
    assert Chem.MolToSmiles(standardized_mol) == original_smiles
