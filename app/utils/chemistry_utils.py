import uuid
import yaml

from rdkit import Chem
from rdkit.Chem import RegistrationHash
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.RegistrationHash import HashLayer

from app.setup.database import SessionLocal
from app import models


class MoleculeStandardizationConfig:
    def __init__(self):
        self._config = None

    @property
    def config(self) -> dict:
        if self._config is None:
            db = SessionLocal()
            try:
                setting = (
                    db.query(models.Settings).filter(models.Settings.name == "Molecule standardization rules").first()
                )
                if not setting:
                    raise Exception("Molecule standardization config not found in settings table.")
                self._config = yaml.safe_load(setting.value)
            finally:
                db.close()
        return self._config


molecule_standardization_config = MoleculeStandardizationConfig()


def standardize_mol(
    mol: Chem.Mol,
) -> Chem.Mol:
    """
    Standardizes a given RDKit molecule using operations defined in the
    molecule standardization settings.

    The operations are dynamically executed in the order defined in the setting, but only if they are enabled.

    Args:
        mol (Chem.Mol): The molecule to standardize.

    Returns:
        Chem.Mol: The standardized molecule after performing all configured operations.
    """
    config = molecule_standardization_config.config

    # Apply only the enabled operations in the order of declaration in the config.
    for operation in config.get("operations", []):
        operation_type = operation.get("type")
        is_enabled = operation.get("enable", True)

        if not is_enabled:
            continue

        if not operation_type:
            raise ValueError("Operation type is missing in the configuration.")

        mol = apply_standardizer_operation(mol, operation_type)

    return mol


def apply_standardizer_operation(mol: Chem.Mol, operation_type: str) -> Chem.Mol:
    """
    Applies a specific operation to the molecule based on the operation type.

    Args:
        mol (Chem.Mol): The molecule to modify.
        operation_type (str): The type of standardization operation to perform.

    Returns:
        Chem.Mol: The transformed molecule.
    """
    operation_map = {
        "Cleanup": rdMolStandardize.Cleanup,
        "FragmentParent": rdMolStandardize.FragmentParent,
        "RemoveHs": Chem.RemoveHs,
        "Uncharger": lambda mol: rdMolStandardize.Uncharger().uncharge(mol),
    }

    if operation_type not in operation_map:
        raise ValueError(f"Unknown operation type: {operation_type}")

    return operation_map[operation_type](mol)


def generate_hash_layers(mol: Chem.Mol) -> dict:
    """
    Generate layers for a given molecule.

    This function calculates the layers using the `RegistrationHash` module.

    Args:
        mol: An RDKit molecule object (`rdkit.Chem.Mol`) for which the layers
                  will be generated.

    Returns:
        dict: A dictionary containing the layers used to compute the MolHash.
    """

    return RegistrationHash.GetMolLayers(mol, enable_tautomer_hash_v2=True)


def generate_uuid_from_string(input_string: str) -> uuid.UUID:
    """
    Generate a UUID hash for a given input string, for hashing different molecule layers.

    Args:
        input_string (str): The input string to hash.

    Returns:
        uuid.UUID: The UUID hash of the input string, ready for PostgreSQL UUID type.
    """
    return uuid.uuid5(uuid.NAMESPACE_DNS, input_string)


def calculate_tautomer_hash(mol: Chem.Mol) -> str:
    """
    Calculate the tautomer hash for a given molecule.
    """
    return generate_uuid_from_string(generate_hash_layers(mol)[HashLayer.TAUTOMER_HASH])


def calculate_no_stereo_smiles_hash(mol: Chem.Mol) -> str:
    """Calculate the no-stereo SMILES hash for a given molecule."""

    return generate_uuid_from_string(generate_hash_layers(mol)[HashLayer.NO_STEREO_SMILES])


def calculate_no_stereo_tautomer_hash(mol: Chem.Mol) -> str:
    """
    Calculate the no-stereo tautomer hash for a given molecule.
    """
    return generate_uuid_from_string(generate_hash_layers(mol)[HashLayer.NO_STEREO_TAUTOMER_HASH])
