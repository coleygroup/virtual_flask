
from pysisyphus.optimizers.FIRE import FIRE as Opt
from pysisyphus.Geometry import Geometry as Geo2
from pysisyphus.elem_data import INV_ATOMIC_NUMBERS
from pysisyphus.constants import ANG2BOHR
import contextlib
import os
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import numpy as np

def optimize_geometry(geom, calc):
    '''
    Optimize a geometry with pysisyphus and return the energy

    Returns:
        geom: pysisyphus.Geometry
    '''
    geom.set_calculator(calc)
    opt = Opt(geom, max_cycles=1)
    with contextlib.redirect_stdout(open(os.devnull, 'w')):
        opt.run()
    energy = opt.geometry.energy
    return energy

def rdmol2geom(mol, confId, coord_type='cart'):
    '''
    Convert a rdkit molecule to a pysisyphus geometry

    Returns:
        geom: pysisyphus.Geometry
    '''
    conf = mol.GetConformer(confId)
    coord = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]).flatten() * ANG2BOHR
    atoms = [INV_ATOMIC_NUMBERS[a.GetAtomicNum()].lower() for a in mol.GetAtoms()]
    geom = Geo2(atoms, coord, coord_type=coord_type)
    return geom

def generate_conformer(smiles):
    '''
    Generate a conformer for a given smiles string using RDKit

    Returns:
        mol: rdkit.Chem.rdchem.Mol
    '''
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol

def calc2(smiles, calculator):
    mol = generate_conformer(smiles)
    geom = rdmol2geom(mol, 0)
    energy = optimize_geometry(geom, calculator)
    return energy