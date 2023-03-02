import argparse
from openbabel import openbabel as ob
from rdkit import Chem



def get_absolute_configuration_from_xyz(xyz_file):
    # Convert XYZ to Mol
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("xyz", "mol")
    mol = ob.OBMol()
    conv.ReadFile(mol, xyz_file)

    # Convert Mol to RDKit Mol
    mol_str = conv.WriteString(mol)
    rdkit_mol = Chem.MolFromMolBlock(mol_str, removeHs=False)

    # Get the absolute configuration
    Chem.AssignAtomChiralTagsFromStructure(rdkit_mol)
    absolute_conf = Chem.FindMolChiralCenters(rdkit_mol, includeUnassigned=True)

    return absolute_conf



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get the absolute configuration from XYZ file. Number of atoms start at 0.")
    parser.add_argument("xyz_file", help="Path to the XYZ file")
    args = parser.parse_args()

    absolute_conf = get_absolute_configuration_from_xyz(args.xyz_file)
    print(absolute_conf)
