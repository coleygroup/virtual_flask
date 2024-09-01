# from pysisyphus.optimizers.FIRE import FIRE as Opt
# from pysisyphus.Geometry import Geometry as Geo2
# from pysisyphus.elem_data import INV_ATOMIC_NUMBERS
# from pysisyphus.constants import ANG2BOHR
# import contextlib
# import os
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

# import numpy as np
import subprocess
import re


# def optimize_geometry(geom, calc):
#     """
#     Optimize a geometry with pysisyphus and return the energy

#     Returns:
#         geom: pysisyphus.Geometry
#     """
#     geom.set_calculator(calc)
#     opt = Opt(geom, max_cycles=300)
#     with contextlib.redirect_stdout(open(os.devnull, "w")):
#         opt.run()
#     energy = opt.geometry.energy
#     return energy


# def rdmol2geom(mol, confId, coord_type="cart"):
#     """
#     Convert a rdkit molecule to a pysisyphus geometry

#     Returns:
#         geom: pysisyphus.Geometry
#     """
#     conf = mol.GetConformer(confId)
#     coord = (
#         np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]).flatten()
#         * ANG2BOHR
#     )
#     atoms = [INV_ATOMIC_NUMBERS[a.GetAtomicNum()].lower() for a in mol.GetAtoms()]
#     geom = Geo2(atoms, coord, coord_type=coord_type)
#     return geom


def generate_conformer(smiles):
    """
    Generate a conformer for a given smiles string using RDKit

    Returns:
        mol: rdkit.Chem.rdchem.Mol
    """
    mol1 = Chem.MolFromSmiles(smiles, sanitize=False)
    Chem.SanitizeMol(
        mol1,
        sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
        ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY,
    )
    mol = Chem.AddHs(mol1)
    AllChem.EmbedMolecule(mol)
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except:
        print("Failed to gen conf", smiles)
        return None
    return mol


# def calc2(smiles, calculator):
#     # print(smiles)
#     mol = generate_conformer(smiles)
#     # print(Chem.MolToSmiles(mol))
#     geom = rdmol2geom(mol, 0)
#     energy = optimize_geometry(geom, calculator)
#     return energy


def smiles_to_3d(smiles):
    """Generate 3D coordinates for a molecule from a SMILES string using RDKit."""
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMoleculeConf(mol, ignoreInterfragInteractions=False)
    return mol


def write_xyz(mol, filename):
    """Write the molecule to an XYZ file."""
    xyz = Chem.MolToXYZBlock(mol)
    with open(filename, "w") as file:
        file.write(xyz)


def run_xtb(filename, charge):
    """Run xTB calculation on a given XYZ file."""
    result = subprocess.run(
        ["xtb", filename, "--opt", "--gbsa", "DMSO", "--chrg", charge],
        capture_output=True,
        text=True,
    )

    # result = subprocess.run(['xtb', filename], capture_output=True, text=True)

    output = result.stdout
    # print(output)
    # Regular expression to find the line containing total energy
    energy_line = re.search(r"TOTAL ENERGY\s+([-0-9.]+) Eh", output)
    if energy_line:
        total_energy = energy_line.group(1)
        return total_energy, output

    return 0, output

import os
def calc(smiles, id):
    # print(smiles)
    mol = generate_conformer(smiles)
    if mol is None:
        return -1000000
    xyz_filename = f"{os.getcwd()}/virtual_flask/qm_calcs/{str(id)}.xyz"
    # print(os.getcwd())
    write_xyz(mol, xyz_filename)
    charge = Chem.GetFormalCharge(mol)

    total_energy, xtb_output = run_xtb(xyz_filename, str(charge))
    return float(total_energy)
