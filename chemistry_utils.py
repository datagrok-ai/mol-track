import yaml
from rdkit import Chem
from rdkit.Chem import RegistrationHash
from rdkit.Chem.MolStandardize import rdMolStandardize


def standardize_mol(
    mol: Chem.Mol,
    standardizer_config_file: str = "molecule_standarizer_operations.yaml",
) -> Chem.Mol:
    """
    Standardizes a given RDKit molecule using operations defined in a YAML configuration file.

    The operations are dynamically executed in the order defined in the file, but only if they are enabled.

    Args:
        mol (Chem.Mol): The molecule to standardize.
        standardizer_config_file (str): Path to the YAML configuration file with operation definitions.

    Returns:
        Chem.Mol: The standardized molecule after performing all configured operations.
    """

    with open(standardizer_config_file, "r") as f:
        config = yaml.safe_load(f)

    # Apply only the enabled operations in the order of declaration in the yml file.
    for operation in config.get("operations", []):
        operation_type = operation.get("type")
        is_enabled = operation.get("enable", True)

        if not is_enabled:
            continue

        if not operation_type:
            raise ValueError(
                f"Operation type is missing in the configuration file:{standardizer_config_file}."
            )

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


def generate_molhash(mol: Chem.Mol) -> tuple:
    """
    Generate a molecular hash and corresponding layers for a given molecule.

    This function calculates a molecular hash (MolHash) and its associated layers using
    the `RegistrationHash` module.

    Args:
        mol: An RDKit molecule object (`rdkit.Chem.Mol`) for which the MolHash 
                  and layers will be generated.

    Returns:
        tuple: A tuple containing:
            - mhash (str): The molecular hash string for the input molecule.
            - layers (dict): A dictionary containing the layers used to compute the MolHash.
    """
    # Generate molecular layers
    layers = RegistrationHash.GetMolLayers(molecule, enable_tautomer_hash_v2=True)
    
    # Generate the molecular hash based on the layers
    mhash = RegistrationHash.GetMolHash(layers)
    
    return mhash, layers