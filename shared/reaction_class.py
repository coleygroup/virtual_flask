from rdkit import Chem
from rdkit.Chem import AllChem
from rdcanon import canon_reaction_smarts
from shared.util import (
    save_img,
    find_convex_hull,
    remove_atom_mapping,
    split_at_period_not_in_parentheses,
    get_networkx_graph_from_state_network,
    get_path_from_init_to_node,
    remove_atom_mapping_mol,
)
import re
import dill
import time
from shared.templates.templates import templates
from multiprocessing import Pool
import json
from psycopg2.extras import execute_batch
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import random
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


comp_to_color = {"A": (1, 0, 0, 0.5), "B": (0, 1, 0, 0.5), "C": (0, 0, 1, 0.5)}

structures = "*1~*~*~*2~*~*~*~*~2~*~1.*1~*~*~*2~*~*~*~*~*~2~*~1.*1~*~*~*2~*(~*~1)~*~*~*1~*~*~*~*~1~2.*1~*~*~*2~*(~*~1)~*~*1~*~*~*~*~*~1~2.*1~*~*~*2~*(~*~1)~*~*~*1~*~*~*~*~*~1~2.*1~*~*~*2~*~*3~*~*~*~*~3~*~*~2~*~1.*1~*~*~*~*~*~1.*1~*~*~*~*~1.*1~*~*~*2~*~*3~*~*~*~*~*~3~*~*~2~*~1.*1~*~*2~*~*~*~*~2~*~1.*1~*~*~*2~*~*~*~*~*~2~*~*~1.*1~*~*~*2~*(~*~1)~*~*~*1~*3~*~*~*~*~3~*~*~*~2~1.*1~*~*~*2~*(~*~1)~*~*1~*~*~*~*~1~2.*1~*~*2~*~*~*3~*~*~*~*~3~*~2~*~1.*1~*~*~*2~*~*~*~*~2~*~*~1.*1~*~*~*2(~*~*~1)~*~*~*~*~2.*1~*~*~*2~*~*3~*~*~*~*~*~3~*~*~*~2~*~1.*1~*~*~*2~*~*3~*(~*~*~2~*~1)~*~*1~*~*~*~*~*~1~3.*1~*~*2~*~*3~*~*~*~*~3~*~*~2~*~1.*1~*~*~*2~*(~*~1)~*~*1~*3~*~*~*~*~*~3~*~*~*~2~1.*1~*~*~*2~*(~*~1)~*~*~*~*1~*~*~*~*~1~2.*1~*~*2~*~*~*~*3~*~*~*(~*~1)~*~2~3.*1~*~*~*2(~*~*~1)~*~*~*~*~*~2.*1~*~*~*~*~*~*~1.*1~*~*~*2~*~*3~*(~*~*~*4~*~*~*~*~*~4~3)~*~*~2~*~1.*1~*~*~*2(~*~*~1)~*~*~*1~*~*~*~*~*~1~2.*1~*~*~*2~*~*3~*~*~*~*~*~3~*~2~*~*~1.*1~*~*2~*~*~*(~*~1)~*~2.*1~*~*~*2(~*~*~1)~*~*~*1~*~*~*~*~*~1~*~2.*1~*~*~*2~*(~*~1)~*~*~*1~*3~*~*~*~*~*~3~*~*~*~2~1.*1~*~*~*2~*(~*~1)~*~*~*1~*~*3~*~*~*~*~*~3~*~1~2.*1~*~*~*2(~*~1)~*~*~*~*~2.*1~*~*2~*~*~*~*3~*~*~*~*(~*~1)~*~2~3.*1~*~*~*2~*(~*~1)~*~*~*~21~*~*~*~*~1.*1~*~*~*2(~*~*~1)~*~*~*~*1~*~*~*~*~*~1~2.*1~*~*2~*~*~*~*(~*~1)~*~2.*1~*~*~*2~*(~*~1)~*~*1~*~*~*3~*~*~*~*~3~*~1~2.*1~*~*2~*~*3~*~*4~*~*~*~*~*~4~*~3~*~*~2~*~1.*1~*~*~*2~*~*3~*~*~*~*~3~*~*~*~2~*~1.*1~*~*~*2(~*~*~1)~*~*~*~2.*1~*~*~*2~*~*~*3~*~*~*~*(~*~1)~*~2~3.*1~*~*~*2~*(~*~1)~*~*~*~*~21~*~*~*~*~1.*1~*~*~*2~*~*3~*~*4~*~*~*~*~*~4~*~*~3~*~*~2~*~1.*1~*~*2~*~*~*~1~*~2.*1~*~*~*2~*~*~*3~*~*~*~*~*~3~*~2~*~*~1.*1~*~*2~*3~*~*~*(~*~3)~*~2~*~1.*1~*~*~*2~*(~*~1)~*~*~*~*1~*~*~*~*~*~1~2.*1~*~*~*2(~*~1)~*~*~*~2"
smarts_rings = [Chem.MolFromSmarts(x) for x in structures.split(".")]


def GetRingSystems(mol, includeSpiro=True):
    ri = mol.GetRingInfo()
    systems = []
    for ring in ri.AtomRings():
        ringAts = set(ring)
        nSystems = []
        for system in systems:
            nInCommon = len(ringAts.intersection(system))
            if nInCommon and (includeSpiro or nInCommon > 1):
                ringAts = ringAts.union(system)
            else:
                nSystems.append(system)
        nSystems.append(ringAts)
        systems = nSystems
    return systems


def mark_if_ring_is_unacceptable(mol):
    Chem.SanitizeMol(
        mol,
        sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE
        ^ Chem.SanitizeFlags.SANITIZE_SYMMRINGS
        ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY,
    )

    ri = GetRingSystems(mol)

    if len(ri) == 0:
        return False

    carbonyl = Chem.MolFromSmarts("C(=O)")
    dicarbonyl = Chem.MolFromSmarts("C(=O)C(=O)")
    sulfoxide = Chem.MolFromSmarts("S(=O)=O")

    for r in ri:
        mol_copy = Chem.RWMol()
        old_idx_to_new = {}
        idx = 0
        for atm in r:
            mol_copy.AddAtom(mol.GetAtomWithIdx(atm))
            old_idx_to_new[atm] = idx
            idx = idx + 1
        for atm in r:
            atom = mol.GetAtomWithIdx(atm)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in r:
                    continue
                mol_copy.AddAtom(nbr)
                old_idx_to_new[nbr.GetIdx()] = idx
                idx = idx + 1
                old_bond = mol.GetBondBetweenAtoms(atm, nbr.GetIdx())
                mol_copy.AddBond(
                    old_idx_to_new[atom.GetIdx()],
                    old_idx_to_new[nbr.GetIdx()],
                    old_bond.GetBondType(),
                )
        for bond in mol.GetBonds():
            if bond.GetBeginAtomIdx() in r and bond.GetEndAtomIdx() in r:
                mol_copy.AddBond(
                    old_idx_to_new[bond.GetBeginAtomIdx()],
                    old_idx_to_new[bond.GetEndAtomIdx()],
                    bond.GetBondType(),
                )

        if len(mol_copy.GetSubstructMatches(carbonyl)) > 2:
            return True
        if len(mol_copy.GetSubstructMatches(dicarbonyl)) > 0:
            return True
        if len(mol_copy.GetSubstructMatches(sulfoxide)) > 1:
            return True

    hits = []
    for i in smarts_rings:
        o = mol.GetSubstructMatches(i)
        for h in o:
            anums = set(h)
            if anums in hits:
                continue
            if anums in ri:
                for anidx, an in enumerate(h):
                    atom_neighbors = mol.GetAtomWithIdx(an).GetNeighbors()
                    neighbors_in_ring = [
                        x.GetIdx() for x in atom_neighbors if x.GetIdx() in anums
                    ]
                    atom_neighbors_smarts = i.GetAtomWithIdx(anidx).GetNeighbors()
                    expected_neighbors_in_ring = [
                        x.GetIdx() for x in atom_neighbors_smarts
                    ]
                    if len(expected_neighbors_in_ring) != len(neighbors_in_ring):
                        return True
                hits.append(set(h))

    if len(hits) == len(ri):
        return False

    return True


class MechanisticReaction:
    def __init__(
        self,
        input_template_smarts,
        reaction_name,
        specified_scope,
        description,
        canonicalize=True,
    ):
        self.input_template = input_template_smarts
        self.name = reaction_name
        self.scope = specified_scope
        self.description = description

        if canonicalize:
            self.canon_template = canon_reaction_smarts(
                input_template_smarts, True, "drugbank", True
            )
        else:
            self.canon_template = input_template_smarts

        self.canon_i_template = self.convert_to_intramolecular(self.canon_template)
        self.canon_r_template = self.convert_to_reverse(self.canon_template)
        self.rxn_rd = AllChem.ReactionFromSmarts(self.canon_template)
        self.rxn_i_rd = AllChem.ReactionFromSmarts(self.canon_i_template)
        self.num_reactants = self.rxn_rd.GetNumReactantTemplates()
        self.num_products = self.rxn_rd.GetNumProductTemplates()
        self.product_templates = [i for i in self.rxn_rd.GetProducts()]
        self.reactant_templates = [i for i in self.rxn_rd.GetReactants()]
        self.product_templates_i = [i for i in self.rxn_i_rd.GetProducts()]
        self.reactant_templates_i = [i for i in self.rxn_i_rd.GetReactants()]

    def serialize(self):
        return {
            "name": self.name,
            "scope": self.scope,
            "description": self.description,
            "template": self.canon_template,
            "input_template": self.input_template,
        }

    def convert_to_intramolecular(self, template):
        reactants = template.split(">>")[0]
        products = template.split(">>")[1]
        if reactants[0] != "(":
            reactants = "(" + reactants + ")"

        if products[0] != "(":
            products = "(" + products + ")"

        i_template = reactants + ">>" + products
        return i_template

    def convert_to_reverse(self, template):
        r_template = ""

        pdts = split_at_period_not_in_parentheses(template.split(">>")[1])
        for p in pdts:
            r_template += p + "."
        r_template = r_template[:-1] + ">>"
        rxts = template.split(">>")[0].split(".")
        for r in rxts:
            r_template += r + "."
        r_template = r_template[:-1]

        return r_template

    def filter_reaction(self, product_sm, ring_filter=False):
        bad_state = False
        output_products = []
        output_product_mols = []
        for pmol in product_sm:
            if ring_filter:
                if mark_if_ring_is_unacceptable(pmol):
                    bad_state = True
                    break

            charges = 0
            for atom in pmol.GetAtoms():
                if atom.GetFormalCharge() != 0:
                    charges += 1

            if charges > 3:
                bad_state = True
                break

            try:
                Chem.SanitizeMol(
                    Chem.Mol(pmol), sanitizeOps=Chem.SanitizeFlags.SANITIZE_KEKULIZE
                )
            except:
                bad_state = True
                break

            output_products.append(Chem.MolToSmiles(pmol))

        if bad_state:
            return None

        return output_products

    def get_highlighed_product_image(self, product_smiles, reacting_atoms, filter={}):
        prod_mol = Chem.MolFromSmiles(product_smiles, sanitize=False)
        atom_cols = {}
        reacting_atoms_color = []
        for atom in prod_mol.GetAtoms():
            a_map = atom.GetAtomMapNum()
            if a_map in reacting_atoms:
                if len(filter) == 0:
                    atom_cols[atom.GetIdx()] = (1, 0, 0, 0.5)
                    reacting_atoms_color.append(atom.GetIdx())
                else:
                    for comp in filter:
                        if a_map in filter[comp]:
                            atom_cols[atom.GetIdx()] = comp_to_color[comp]
                            reacting_atoms_color.append(atom.GetIdx())
                            break

        img = save_img(
            product_smiles,
            highlight_atoms=reacting_atoms_color,
            atom_cols=atom_cols,
        )

        return img

    def _prepare_reaction_i(self, reactants):

        return (
            self.rxn_i_rd,
            reactants,
            self.product_templates_i,
            self.reactant_templates_i,
        )

    def _prepare_reaction(self, reactants, needs_map):
        if len(reactants) != self.num_reactants:
            return None, None, None, None

        rxn_obj = self.rxn_rd
        pt = self.product_templates
        rt = self.reactant_templates

        if type(reactants[0]) == str:
            reactant_objs = [Chem.MolFromSmiles(i, sanitize=False) for i in reactants]
        else:
            reactant_objs = reactants

        if needs_map:
            amap = 1
            for r in reactant_objs:
                for atom in r.GetAtoms():
                    atom.SetAtomMapNum(amap)
                    amap += 1

        return rxn_obj, reactant_objs, pt, rt

    def _get_reactant_map(self, reactant_objs):
        atom_map_to_reactant_atom = {}
        for robj in reactant_objs:
            Chem.SanitizeMol(robj, sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE)
            for a in robj.GetAtoms():
                atom_map_to_reactant_atom[a.GetAtomMapNum()] = a

        return atom_map_to_reactant_atom

    def _compile_products(self, products, rxn_obj, reactant_objs, isomeric, pt, rt):

        atom_map_to_reactant_atom = self._get_reactant_map(reactant_objs)

        output_products = []
        output_reacting_atoms = []
        output_mols = []
        for prod_set_i, p in enumerate(products):
            this_reaction_prods = []
            this_reacting_atoms = []
            this_reaction_mols = []
            for prod_index, prod_mol in enumerate(p):
                Chem.SanitizeMol(
                    prod_mol,
                    sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE,
                )
                atom_map_to_product_atom = {}
                query_map_to_atom_map = {}
                self.adjust_atom_maps(
                    prod_mol,
                    rxn_obj,
                    reactant_objs,
                    atom_map_to_product_atom,
                    query_map_to_atom_map,
                )
                self.fix_protons(
                    atom_map_to_product_atom,
                    query_map_to_atom_map,
                    rt,
                    pt[prod_index],
                    atom_map_to_reactant_atom,
                )
                reacting_atoms = self.get_reacting_atoms(
                    prod_mol, atom_map_to_reactant_atom
                )
                if not reacting_atoms:
                    continue

                smiles = Chem.MolToSmiles(prod_mol, isomericSmiles=isomeric)
                this_reaction_prods.append(smiles)
                this_reacting_atoms.extend(reacting_atoms)
                this_reaction_mols.append(prod_mol)
            output_products.append(this_reaction_prods)
            output_reacting_atoms.append(this_reacting_atoms)
            output_mols.append(this_reaction_mols)

        return output_products, output_reacting_atoms, output_mols

    def run(self, reactants, isomeric=True, needs_map=False, intramolecular=False):

        if intramolecular:
            rxn_obj, reactant_objs, pt, rt = self._prepare_reaction_i(reactants)
        else:
            rxn_obj, reactant_objs, pt, rt = self._prepare_reaction(
                reactants, needs_map
            )

            if not rxn_obj:
                return None, None, None

        products = self.rd_run_reactants(rxn_obj, reactant_objs)

        if not products and not intramolecular:
            if len(reactants) == 1:
                return None, None, None

            reactant_objs.reverse()
            products = self.rd_run_reactants(rxn_obj, reactant_objs)

        if not products:
            return None, None, None

        output_products, output_reacting_atoms, output_mols = self._compile_products(
            products, rxn_obj, reactant_objs, isomeric, pt, rt
        )

        return output_products, output_reacting_atoms, output_mols

    def rd_run_reactants(self, rx_obj, rt_objs):
        return rx_obj.RunReactants(rt_objs)

    def _get_nested_text(self, input_string):
        """
        Splits the input_string by the specified delimiter, but does not split anything within the nested structures.
        """
        parts = []
        current = []
        nested_level = 0
        i = 0

        nested_start = "$("
        nested_end = ")"

        delimiter = [";", "&", ","]

        while i < len(input_string):
            if input_string.startswith(nested_start, i):
                nested_level += 1
                current.append(input_string[i])
            elif input_string.startswith("(", i) and not input_string.startswith(
                nested_start, i - 1
            ):
                nested_level += 1
                current.append(input_string[i])
            elif input_string.startswith(nested_end, i) and nested_level:
                nested_level -= 1
                current.append(input_string[i])
            elif input_string[i] in delimiter and nested_level == 0:
                if current:
                    parts.append("".join(current))
                    current = []
            else:
                if nested_level != 0:
                    current.append(input_string[i])
            i += 1
        if current:
            parts.append("".join(current))

        return parts

    def fix_protons(
        self,
        atom_map_to_product_atom,
        query_map_to_atom_map,
        reactant_templates,
        product_template,
        atom_map_to_reactant_atom,
    ):
        atom_map_to_reactant_smarts = {}
        atom_map_to_product_smarts = {}

        for rt in reactant_templates:
            for a in rt.GetAtoms():
                if a.GetAtomMapNum() in query_map_to_atom_map:
                    atom_map_to_reactant_smarts[
                        query_map_to_atom_map[a.GetAtomMapNum()]
                    ] = a

        for a in product_template.GetAtoms():
            if a.GetAtomMapNum() in query_map_to_atom_map:
                atom_map_to_product_smarts[query_map_to_atom_map[a.GetAtomMapNum()]] = a

        # fix protons
        for match in atom_map_to_product_atom:
            if match not in atom_map_to_reactant_atom:
                continue
            reactant_atom = atom_map_to_reactant_atom[match]
            product_atom = atom_map_to_product_atom[match]
            reactant_template_atom = atom_map_to_reactant_smarts[match]
            product_template_atom = atom_map_to_product_smarts[match]

            r_smarts = reactant_template_atom.GetSmarts()
            p_smarts = product_template_atom.GetSmarts()

            # print(r_smarts, p_smarts)
            recs = self._get_nested_text(r_smarts)
            for rec_patt in recs:
                r_smarts = r_smarts.replace(rec_patt, "")
                p_smarts = p_smarts.replace(rec_patt, "")
            # print(r_smarts, p_smarts)
            # print()

            pattern = r"H\d+(?:,H\d+)*"
            match_before = re.search(pattern, r_smarts)

            match_after = re.search(pattern, p_smarts)

            hydrogens_before = None
            hydrogens_after = None
            if match_before:
                hydrogens_before = match_before[0]
            if match_after:
                hydrogens_after = match_after[0]

            # print(hydrogens_before, hydrogens_after)
            if (
                not hydrogens_before
                and hydrogens_after
                or hydrogens_before
                and not hydrogens_after
            ):
                print("error")
            if not hydrogens_after and not hydrogens_before:
                continue

            first_h_before = int(hydrogens_before.split(",")[0][1])
            first_h_after = int(hydrogens_after.split(",")[0][1])

            if first_h_before == first_h_after:
                product_atom.SetNumExplicitHs(reactant_atom.GetTotalNumHs())
            elif first_h_after - 1 == first_h_before:
                product_atom.SetNumExplicitHs(reactant_atom.GetTotalNumHs() + 1)
            elif first_h_after + 1 == first_h_before:
                product_atom.SetNumExplicitHs(reactant_atom.GetTotalNumHs() - 1)

    def adjust_atom_maps(
        self, i, rxn, rxt_in, atom_map_to_product_atom, query_map_to_atom_map
    ):
        for atom in i.GetAtoms():
            if (
                "molAtomMapNumber" not in atom.GetPropsAsDict()
                and "old_mapno" in atom.GetPropsAsDict()
            ):
                old_atom_idx = atom.GetPropsAsDict()["old_mapno"]
                react_atom_idx = atom.GetPropsAsDict()["react_atom_idx"]
                hit = 0
                for i2, rx in enumerate(rxn.GetReactants()):
                    found = False
                    for i3, j in enumerate(rx.GetAtoms()):
                        if j.GetAtomMapNum() == old_atom_idx:
                            found = True
                            break
                    if found:
                        hit = i2
                        break
                atom_map = rxt_in[hit].GetAtomWithIdx(react_atom_idx).GetAtomMapNum()
                atom.SetAtomMapNum(atom_map)

                atom_map_to_product_atom[atom_map] = atom
                query_map_to_atom_map[old_atom_idx] = atom_map

    def get_reacting_atoms(self, i, atom_map_to_reactant_atom):
        rxt_atoms = []
        for atom in i.GetAtoms():
            # if "old_mapno" in atom.GetPropsAsDict():
            #     old_atom_idx = atom.GetPropsAsDict()["old_mapno"]
            #     react_atom_idx = atom.GetPropsAsDict()["react_atom_idx"]
            #     hit = 0
            #     for i2, rx in enumerate(rxn.GetReactants()):
            #         found = False
            #         for i3, j in enumerate(rx.GetAtoms()):
            #             if j.GetAtomMapNum() == old_atom_idx:
            #                 found = True
            #                 break
            #         if found:
            #             hit = i2
            #             break

            # old_atom = rxt_in[hit].GetAtomWithIdx(react_atom_idx)
            old_atom = atom_map_to_reactant_atom[atom.GetAtomMapNum()]
            old_atom_info = [
                old_atom.GetFormalCharge(),
                old_atom.GetNumExplicitHs(),
                old_atom.GetNumImplicitHs(),
                old_atom.GetDegree(),
                old_atom.GetTotalDegree(),
                old_atom.GetTotalValence(),
                old_atom.GetTotalNumHs(),
                old_atom.GetExplicitValence(),
                old_atom.GetImplicitValence(),
            ]
            # old_atom_neighbors = [x.GetIdx() for x in old_atom.GetNeighbors()]
            # old_atom_neighbors.sort()
            new_atom_info = [
                atom.GetFormalCharge(),
                atom.GetNumExplicitHs(),
                atom.GetNumImplicitHs(),
                atom.GetDegree(),
                atom.GetTotalDegree(),
                atom.GetTotalValence(),
                atom.GetTotalNumHs(),
                atom.GetExplicitValence(),
                atom.GetImplicitValence(),
            ]

            if (
                new_atom_info[0] == old_atom_info[0]
                and new_atom_info[4] != old_atom_info[4]
                and new_atom_info[5] != old_atom_info[5]
                and new_atom_info[6] != old_atom_info[6]
            ):
                old_bonds = [x.GetBondType() for x in old_atom.GetBonds()]
                new_bonds = [x.GetBondType() for x in atom.GetBonds()]

                if old_bonds == new_bonds:
                    # print(Chem.MolToSmiles(rxt_in[0]))
                    # print("yolo")
                    # print(Chem.MolToSmiles(i))
                    # print(atom.GetSymbol())
                    # print(new_atom_info)
                    # print(old_atom_info)
                    # print(old_bonds)
                    # print(new_bonds)
                    # print()
                    # basically in situations where ([N;H1,H2;+0:1]-[c;H0;+0:2].[O;H0;+0:3]=[C;!$([C](=[O])[O]);H0;+0:4])>>([N;+1;H1,H2:1](-[c;H0;+0:2])-[C;!$([C](=[O])[O]);H0;+0:4]-[O;-1;H0:3])
                    # the intramolecular reaction forms a bond that's already there.
                    # rdkit will run the reaction anyways, and mess up the protonation state,
                    # creating a molecule with improper valence
                    return None

            if (
                new_atom_info[3] == old_atom_info[3]
                and new_atom_info[4] == old_atom_info[4]
                and new_atom_info[5] == old_atom_info[5]
                and new_atom_info[6] == old_atom_info[6]
                and new_atom_info[0] != old_atom_info[0]
            ):
                old_bonds = [x.GetBondType() for x in old_atom.GetBonds()]
                new_bonds = [x.GetBondType() for x in atom.GetBonds()]
                if old_bonds == new_bonds:
                    # basically in situations where ([N;H1,H2;+0:1]-[c;H0;+0:2].[O;H0;+0:3]=[C;!$([C](=[O])[O]);H0;+0:4])>>([N;+1;H1,H2:1](-[c;H0;+0:2])-[C;!$([C](=[O])[O]);H0;+0:4]-[O;-1;H0:3])
                    # the intramolecular reaction forms a bond that's already there.
                    # rdkit will run the reaction anyways, and mess up the protonation state,
                    # creating a molecule with improper valence
                    return None

            # new_atom_neighbors = [x.GetAtomMapNum() for x in atom.GetNeighbors()]
            # new_atom_neighbors.sort()
            if old_atom_info != new_atom_info:
                rxt_atoms.append(atom.GetAtomMapNum())

        return rxt_atoms

    def __str__(self):
        return f"{self.name}"

    def __repr__(self):
        return f"{self.name}"

    def __eq__(self, other):
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)


class StateNode:
    def __init__(self, node_data, save_img_flag=False):
        self.propagations = node_data["propagations"]
        self.mapped_smiles = node_data["mapped_smiles"]
        self.unmapped_smiles = node_data["unmapped_smiles"]
        self.reacting_atoms = node_data["reacting_atoms"]
        self.these_reacting_atoms = node_data["these_reacting_atoms"]
        self.count = node_data["count"]
        self.product_in_precalc = node_data["product_in_precalc"]
        self.these_reacting_atoms_path = node_data["these_reacting_atoms_path"]

        self.tcp = None
        self.known_product = None
        self.askcos_hit = None
        self.failed_thermo_state1 = None
        self.hit_sequence = None
        self.known_sequence = None
        self.tcp_connected = None
        self.in_pubchem = None
        self.pubchem_hits = 0
        self.similarity = 0
        self.average_distance_to_db = 0
        self.in_same_mechanism = None
        self.contains_non_participating_transformation = None
        self.structural_failure = None
        self.seen = None
        self.other_data = {}

        if save_img_flag:
            self.img = save_img(self.mapped_smiles)
            self.hull = find_convex_hull(self.img)
        else:
            self.img = ""
            self.hull = []

    def __str__(self):
        return f"Smiles: {self.unmapped_smiles}"

    def __repr__(self):
        return f"Smiles: {self.unmapped_smiles}"

    def __getitem__(self, node_id):
        return self.nodes[node_id]

    def __setitem__(self, attr, val):
        self.nodes[attr] = val

    def serialize(self):
        return {
            "mapped_smiles": self.mapped_smiles,
            "unmapped_smiles": self.unmapped_smiles,
            "reacting_atoms": self.reacting_atoms,
            "these_reacting_atoms": self.these_reacting_atoms,
            "count": self.count,
            "product_in_precalc": self.product_in_precalc,
            "these_reacting_atoms_path": self.these_reacting_atoms_path,
            "tcp": self.tcp,
            "failed_thermo_state1": self.failed_thermo_state1,
            "similarity": self.similarity,
            "average_distance_to_db": self.average_distance_to_db,
            "in_same_mechanism": self.in_same_mechanism,
            "contains_non_participating_transformation": self.contains_non_participating_transformation,
            "structural_failure": self.structural_failure,
        }


def check_if_in_precalc(input_smiles, precalc, sms):

    for sm in input_smiles.split("."):
        if sm in sms:
            continue
        prod_mol = Chem.MolFromSmiles(sm)
        num_charges = Chem.GetFormalCharge(prod_mol)
        if num_charges != 0:
            continue
        prod_inchi = Chem.MolToInchi(prod_mol)
        prod_mol = Chem.MolFromInchi(prod_inchi)
        prod_smiles = Chem.MolToSmiles(prod_mol, isomericSmiles=False)

        if prod_smiles in precalc:
            return True

    return False


class VirtualFlask:
    def __init__(
        self,
        mechanisms,
    ):
        self.nodes = {}
        self.edges = {}
        self.runtime = {}
        self.cum_runtime = {}
        self.mechanisms = mechanisms
        self.mech_count_map = {}
        self.end_state = None
        self.nx = None
        for m in mechanisms:
            self.mech_count_map[m] = 0

    def clear(self):
        self.nodes = {}
        self.edges = {}
        self.runtime = {}
        self.cum_runtime = {}

        self.mech_count_map = {}
        for m in self.mechanisms:
            self.mech_count_map[m] = 0

    def charge(self, input_smiles_list, reagents, remap=True, inp_atom_maps={}):
        if len(reagents) == 0:
            reagents = []
        else:
            reagents = reagents.split(",")
        self.reagents = reagents

        if remap:
            self.mapped_input_smiles_list, self.starting_material_atom_maps = (
                self.map_input_molecules(input_smiles_list, self.reagents)
            )
        else:
            self.mapped_input_smiles_list = input_smiles_list
            self.starting_material_atom_maps = inp_atom_maps

        self.input_mapped_smiles = Chem.MolToSmiles(
            Chem.MolFromSmiles(".".join(self.mapped_input_smiles_list))
        )
        self.input_unmapped_smiles = remove_atom_mapping(self.input_mapped_smiles)

        node_data = {
            "propagations": 0,
            "mapped_smiles": self.input_mapped_smiles,
            "unmapped_smiles": self.input_unmapped_smiles,
            "reacting_atoms": [],
            "these_reacting_atoms": [],
            "count": 0,
            "product_in_precalc": False,
            "these_reacting_atoms_path": [[]],
        }

        self.add_node(self.input_unmapped_smiles, node_data)

    def save_as_dill(self, filename):
        with open(filename, "wb") as f:
            dill.dump(self, f)

    @staticmethod
    def load_from_dill(filename):
        with open(filename, "rb") as f:
            return dill.load(f)

    def push_to_db(self, conn, table):
        cur = conn.cursor()

        cur.execute(
            f"INSERT INTO {table}_networks (mapped_input_smiles_list, reagents, input_mapped_smiles, input_unmapped_smiles, starting_material_atom_maps, runtime, cum_runtime, mech_count_map, end_state, other_data) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s) RETURNING network_id",
            (
                self.mapped_input_smiles_list,
                self.reagents,
                self.input_mapped_smiles,
                self.input_unmapped_smiles,
                json.dumps(self.starting_material_atom_maps),
                json.dumps(self.runtime),
                json.dumps(self.cum_runtime),
                json.dumps(self.mech_count_map),
                self.end_state,
                json.dumps({}),
            ),
        )

        network_id = cur.fetchone()[0]

        # cur.execute("SELECT pg_get_serial_sequence('test_nodes', 'node_id');")
        # sequence_name = cur.fetchone()[0]
        # cur.execute(f"SELECT last_value FROM {sequence_name};")
        # current_index = cur.fetchone()[0]

        data_to_insert = []
        for n in self.nodes:
            if self.nodes[n].these_reacting_atoms_path == [[]]:
                data_to_insert.append(
                    (
                        network_id,
                        self.nodes[n].mapped_smiles,
                        self.nodes[n].unmapped_smiles,
                        None,
                        self.nodes[n].product_in_precalc,
                        json.dumps({}),
                    ),
                )
            else:
                data_to_insert.append(
                    (
                        network_id,
                        self.nodes[n].mapped_smiles,
                        self.nodes[n].unmapped_smiles,
                        json.dumps(self.nodes[n].these_reacting_atoms_path),
                        self.nodes[n].product_in_precalc,
                        json.dumps({}),
                    ),
                )

            # node_to_id_map[self.nodes[n].unmapped_smiles] = current_index
            # current_index += 1
        # for k in data_to_insert:
        # print(k)
        execute_batch(
            cur,
            f"INSERT INTO {table}_nodes (network_id, mapped_smiles, unmapped_smiles, these_reacting_atoms_path, product_in_precalc, other_data) VALUES (%s, %s, %s, %s::jsonb, %s, %s);",
            data_to_insert,
        )

        node_to_id_map = {}
        cur.execute(
            f"SELECT node_id, unmapped_smiles FROM {table}_nodes WHERE network_id = {network_id};",
        )
        for n in cur.fetchall():
            node_to_id_map[n[1]] = n[0]

        data_to_insert = []
        for n in self.nodes:
            source_node_id = node_to_id_map[n]
            if n not in self.edges:
                continue
            for e in self.edges[n]:
                dest_node_id = node_to_id_map[e]
                for edge in self.edges[n][e]:
                    ed = {
                        "propagations": edge["propagations"],
                        "template_obj": edge["template_obj"].serialize(),
                        "template": edge["template"],
                        "reaction_smiles": edge["reaction_smiles"],
                        "reacting_atoms": edge["reacting_atoms"],
                        "these_reacting_atoms": edge["these_reacting_atoms"],
                    }

                    data_to_insert.append(
                        (network_id, source_node_id, dest_node_id, json.dumps(ed))
                    )
                    # cur.execute(
                    #     f"INSERT INTO test_edges (network_id, source_node_id, destination_node_id, other_data) VALUES (%s, %s, %s, %s);",
                    #     (network_id, source_node_id, dest_node_id, json.dumps(ed)),
                    # )

        execute_batch(
            cur,
            f"INSERT INTO {table}_edges (network_id, source_node_id, destination_node_id, other_data) VALUES (%s, %s, %s, %s);",
            data_to_insert,
        )

        conn.commit()
        cur.close()
        return network_id

    def add_node(self, node_id, node_data):
        # print(node_data["these_reacting_atoms_path"])
        self.nodes[node_id] = StateNode(node_data)

    def increment_mech_count(self, mech_name):
        if mech_name in self.mech_count_map:
            self.mech_count_map[mech_name] += 1
        else:
            self.mech_count_map[mech_name] = 1

    def add_edge(self, node1, node2, edge_data):
        if node1 not in self.edges:
            self.edges[node1] = {}
        if node2 not in self.edges[node1]:
            self.edges[node1][node2] = []
        self.edges[node1][node2].append(edge_data)

    def get_node(self, node_id):
        return self.nodes[node_id]

    def __getitem__(self, node_id):
        return self.nodes[node_id]

    def get_edge(self, node1, node2):
        return self.edges[node1][node2]

    def get_incoming_edges(self, node):
        in_edges = []
        for i in self.edges:
            for ii in self.edges[i]:
                if ii == node:
                    in_edges.append((i, self.edges[i][ii]))

        return in_edges

    def update_tree_with_other_tree(self, other_tree):
        for n in other_tree.nodes:
            if n in self.nodes:
                continue
            self.nodes[n] = other_tree.nodes[n]

        for n1 in other_tree.edges:
            for n2 in other_tree.edges[n1]:
                if n1 in self.edges:
                    if n2 in self.edges[n1]:
                        continue
                if n1 not in self.edges:
                    self.edges[n1] = {}
                if n2 not in self.edges[n1]:
                    self.edges[n1][n2] = []
                self.edges[n1][n2] = other_tree.edges[n1][n2]

    def map_input_molecules(self, inp_mols, reagents):
        mapped_inp_mols = []
        atom_map = 1
        for i in inp_mols:
            mol = Chem.MolFromSmiles(i)
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(atom_map)
                atom_map += 1
            mapped_inp_mols.append(Chem.MolToSmiles(mol))

        for i in reagents:
            mol = Chem.MolFromSmiles(i)
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(atom_map)
                atom_map += 1
            mapped_inp_mols.append(Chem.MolToSmiles(mol))

        # lets = ["A", "B", "C", "D", "E", "F", "G", "H"]
        # starting_material_atom_maps = {}
        # for idx, m in enumerate(mapped_inp_mols):
        #     starting_material_atom_maps[lets[idx]] = [
        #         i.GetAtomMapNum() for i in Chem.MolFromSmiles(m).GetAtoms()
        #     ]

        starting_material_atom_maps = {}
        for idx, m in enumerate(mapped_inp_mols):
            starting_material_atom_maps[m] = [
                i.GetAtomMapNum() for i in Chem.MolFromSmiles(m).GetAtoms()
            ]

        return mapped_inp_mols, starting_material_atom_maps

    def propagate_reactions_in_network(
        this,
        state,
        propagations,
        intramolecular=True,
        ring_filter=True,
    ):
        init_reacting_atoms = this.nodes[state].reacting_atoms

        combos = [[this.nodes[state].mapped_smiles]]
        combo_mols = [
            [Chem.MolFromSmiles(this.nodes[state].mapped_smiles, sanitize=False)]
        ]

        reactant_to_reactant_atom_maps = {}

        for rxtnt in this.nodes[state].mapped_smiles.split("."):
            atm_maps = []
            for atom in Chem.MolFromSmiles(rxtnt, sanitize=False).GetAtoms():
                if atom.GetAtomMapNum() > 0:
                    atm_maps.append(atom.GetAtomMapNum())
            reactant_to_reactant_atom_maps[rxtnt] = atm_maps

        new_nodes = []
        new_edges = []

        for r in this.mechanisms:
            rxn_vf = this.mechanisms[r][0]
            for combo_idx, reactants in enumerate(combos):
                # products = list of product lists; each list contains product
                # smiles for each output possibility of the reaction
                # reacting_atoms = list of lists of atom indices; each list
                # contains the atom map numbers of the reactions's reacting atoms

                products, reacting_atoms, p_mols = rxn_vf.run(
                    combo_mols[combo_idx], intramolecular=intramolecular
                )

                if not products:
                    continue

                for p_idx, prods in enumerate(p_mols):
                    # prods = list of product smiles

                    output_products = rxn_vf.filter_reaction(prods, ring_filter)

                    if not output_products:
                        continue

                    new_reacting_atoms = init_reacting_atoms + reacting_atoms[p_idx]
                    new_reacting_atoms = list(set(new_reacting_atoms))

                    product_atom_maps = []
                    for p in prods:
                        for atom in p.GetAtoms():
                            if atom.GetAtomMapNum() > 0:
                                product_atom_maps.append(atom.GetAtomMapNum())

                    to_add = []
                    for input_sm in reactant_to_reactant_atom_maps:
                        has_input = False
                        for atm in product_atom_maps:
                            if atm in reactant_to_reactant_atom_maps[input_sm]:
                                has_input = True
                                break
                        if not has_input:
                            to_add.append(input_sm)

                    mapped_output_state_smiles = ".".join(output_products + to_add)
                    mapped_output_state_mols = Chem.MolFromSmiles(
                        mapped_output_state_smiles, sanitize=False
                    )
                    reaction_smiles = (
                        this.nodes[state].mapped_smiles
                        + ">>"
                        + mapped_output_state_smiles
                    )

                    remove_atom_mapping_mol(mapped_output_state_mols)
                    unmapped_output_state_smiles = Chem.MolToSmiles(
                        mapped_output_state_mols
                    )

                    this.increment_mech_count(r)
                    # print(this.nodes[state].these_reacting_atoms_path)
                    # if this.nodes[state].these_reacting_atoms_path == [[]]:
                    # print("hello.")
                    # if not this.nodes[state].these_reacting_atoms_path:
                    #     new_path = [reacting_atoms[p_idx]]
                    # else:
                    new_path = this.nodes[state].these_reacting_atoms_path.copy()
                    new_path.append(reacting_atoms[p_idx])
                    node_data = {
                        "propagations": propagations,
                        "reaction_description": rxn_vf.description,
                        "template_name": rxn_vf.name + " | " + rxn_vf.scope,
                        "template_obj": rxn_vf,
                        "count": 1,
                        "energy": None,
                        "template": r,
                        "reaction_smiles": reaction_smiles,
                        "mapped_smiles": mapped_output_state_smiles,
                        "unmapped_smiles": unmapped_output_state_smiles,
                        "reacting_atoms": new_reacting_atoms,
                        "these_reacting_atoms": reacting_atoms[p_idx],
                        "these_reacting_atoms_path": new_path,
                        "product_in_precalc": False,
                    }

                    new_nodes.append((unmapped_output_state_smiles, node_data))
                    new_edges.append((state, unmapped_output_state_smiles, node_data))

                    # if unmapped_output_state_smiles not in this.nodes:
                    #     this.add_node(unmapped_output_state_smiles, node_data)
                    # else:
                    #     this.nodes[unmapped_output_state_smiles].count += 1

                    # # if state == "CC(=O)c1ccccc1.Cc1ccc(NC([OH2+])c2ccccc2)cc1":
                    # #     print(propagations, state, "to\n", unmapped_output_state_smiles)
                    # #     print(rxn_vf.description, rxn_vf.name + " | " + rxn_vf.scope)
                    # this.add_edge(state, unmapped_output_state_smiles, node_data)

        return new_nodes, new_edges

    def worker(this, state, propagations, intramolecular, ring_filter):
        # RDLogger.DisableLog('rdApp.*')
        new_nodes, new_edges = this.propagate_reactions_in_network(
            state,
            propagations + 1,
            intramolecular,
            ring_filter,
        )
        return new_nodes, new_edges

    def run_until_done(
        this,
        iters,
        thresh=None,
        intramolecular=True,
        ring_filter=True,
        precalc_prods=[],
        verbose=False,
    ):
        """
        Run the simulation until a stopping condition is met.

        Args:
            iters (int): The maximum number of iterations to run.
            thresh (int): The maximum number of nodes in the state network before stopping.
            map (bool, optional): False if you want to use your own atom mapping. Defaults to True.
            inp_atom_maps (dict, optional): Input atom maps if you want to use your own. Set map to True. Defaults to {}.
            intramolecular (bool, optional): Allow the simulation to apply the templates in an intramolecular fashion. Defaults to True.
            ring_filter (bool, optional): Whether to filter reactions based on ring membership. Defaults to True.
            verbose (bool, optional): Whether to print verbose output. Defaults to False.

        Returns:
            StateNetwork: The state network after running the simulation.
        """

        propagations = 0
        time_00 = time.time()
        early_stop = False
        for i in range(iters):
            nodes_to_run = [
                n
                for n in this.nodes
                if this[n].propagations == propagations
                and this[n].product_in_precalc == False
            ]
            time_0 = time.time()

            args_list = [
                (state, propagations, intramolecular, ring_filter)
                for state in nodes_to_run
            ]
            # with Pool(8) as pool:
            #     updates = pool.starmap(this.worker, args_list)

            # for nodes, edges in updates:
            for arg in args_list:
                nodes, edges = this.worker(*arg)
                for node in nodes:
                    if node[0] not in this.nodes:

                        try:
                            chk = check_if_in_precalc(
                                node[0], precalc_prods, this.input_unmapped_smiles
                            )
                        except:
                            # print(node)
                            chk = False

                        if chk:
                            node[1]["product_in_precalc"] = True
                            # print("hello??", node[0], node[1]["propagations"])
                        else:
                            node[1]["product_in_precalc"] = False

                        this.add_node(node[0], node[1])
                    else:
                        this.nodes[node[0]].count += 1
                for edge in edges:
                    this.add_edge(edge[0], edge[1], edge[2])

                time_n = time.time()
                if time_n - time_00 > 600:
                    this.end_state = "timeout"
                    early_stop = True
                    break

            # for state in nodes_to_run:
            #     this.propagate_reactions_in_network(
            #         state,
            #         propagations + 1,
            #         intramolecular,
            #         ring_filter,
            #     )

            time_1 = time.time()
            this.runtime[propagations + 1] = time_1 - time_0
            this.cum_runtime[propagations + 1] = time_1 - time_00
            # print(this.cum_runtime[propagations + 1], len(this.nodes))
            if verbose:
                print(propagations, len(this.nodes))
            propagations += 1
            if thresh != None:
                if len(this.nodes) > thresh:
                    this.end_state = "threshold"
                    early_stop = True
                    break
            if this.cum_runtime[propagations] > 600:
                this.end_state = "timeout"
                early_stop = True
                break
        if not early_stop:
            this.end_state = "completed"
        this.nx = get_networkx_graph_from_state_network(this)
        print(this.cum_runtime[propagations], len(this.nodes))

    def find_node_with_smiles(self, smiles):
        sm = Chem.CanonSmiles(smiles)
        for n in self.nodes:
            for s in n.split("."):
                if sm == Chem.CanonSmiles(s):
                    return n

    def get_leaf_nodes(self):
        leaf_nodes = []
        for n in self.nodes:
            if n not in self.edges:
                leaf_nodes.append(n)
                continue
            if len(self.edges[n]) == 0:
                leaf_nodes.append(n)
        return leaf_nodes

    def get_path_from_init_to_node(self, targ):
        for h in self.nodes:
            if self.nodes[h].propagations == 0:
                init_node = h
                break

        path = nx.shortest_path(self.nx, init_node, targ)

        return path

    def get_root_node(self):
        for n in self.nodes:
            if self.nodes[n].propagations == 0:
                return n

    def create_subgraph_with_beam_array(self, beam_array):
        root = self.get_root_node()
        current_level_nodes = [root]
        new_graph = nx.DiGraph()
        new_graph.add_node(root)
        for beam_width in beam_array:
            next_level_nodes = []
            for node in current_level_nodes:
                if node in self.edges:
                    children = list(self.edges[node].keys())[:beam_width]
                    for child in children:
                        if child not in new_graph:
                            new_graph.add_node(
                                child, propagations=self.nodes[child].propagations
                            )
                        new_graph.add_edge(node, child)
                        next_level_nodes.append(child)
            current_level_nodes = next_level_nodes
        return new_graph

    def draw_hypergraph_filters(self, name=None):
        G = self.nx

        pos = nx.nx_agraph.graphviz_layout(G, prog="dot")

        propagation_values = [self.nodes[node].propagations for node in self.nodes]
        norm = mcolors.Normalize(vmin=min(propagation_values), vmax=5)
        # print(matplotlib.colormaps)
        cmap = plt.cm._colormaps.get_cmap("viridis")

        node_colors = [cmap(norm(value)) for value in propagation_values]

        node_size = 0.3

        node_colors_final = []
        s = []
        idxes = []
        for i, k in enumerate(self.nodes):
            if (
                self.nodes[k].tcp == False
                or self.nodes[k].failed_thermo_state1 == True
                or self.nodes[k].contains_non_participating_transformation == True
                or self.nodes[k].structural_failure == True
                or self.nodes[k].tcp == None
                or self.nodes[k].failed_thermo_state1 == None
                or self.nodes[k].contains_non_participating_transformation == None
                or self.nodes[k].structural_failure == None
            ):
                node_colors_final.append("lightgrey")
                s.append(node_size)
                continue
            node_colors_final.append(node_colors[i])
            idxes.append(i)
            s.append(node_size)
        print(len(idxes))
        fig, ax = plt.subplots(dpi=300)
        fig.set_size_inches(2, 2)

        nodes_x, nodes_y = zip(*[pos[n] for n in G.nodes()])
        ax.scatter(
            nodes_x,
            nodes_y,
            s=s,
            color=node_colors_final,
            alpha=1,
            zorder=2,
        )

        for i in idxes:
            ax.scatter(
                nodes_x[i],
                nodes_y[i],
                s=s[i],
                color="gold",
                alpha=1,
                zorder=4,
            )

        for edge in G.edges():
            start_x, start_y = pos[edge[0]]
            end_x, end_y = pos[edge[1]]

            # Calculate the direction vector of the edge
            edge_vec = np.array([end_x - start_x, end_y - start_y])
            edge_length = np.linalg.norm(edge_vec)
            edge_direction = edge_vec / edge_length

            # Shorten the edge to stop at the node boundary instead of the center
            start_x, start_y = (
                start_x + edge_direction[0] * node_size,
                start_y + edge_direction[1] * node_size,
            )
            end_x, end_y = (
                end_x - edge_direction[0] * node_size,
                end_y - edge_direction[1] * node_size,
            )

            if edge[0] in idxes and edge[1] in idxes:
                ax.plot(
                    [start_x, end_x],
                    [start_y, end_y],
                    color="gold",
                    alpha=1,
                    zorder=3,
                    linewidth=0.5,
                )
            else:

                ax.plot(
                    [start_x, end_x],
                    [start_y, end_y],
                    color="black",
                    alpha=0.25,
                    zorder=1,
                    linewidth=0.3,
                )

        plt.axis("off")

        if name == None:
            plt.show()
        else:
            plt.savefig(
                f"{name}.png",
                dpi=900,
                transparent=True,
                bbox_inches="tight",
                pad_inches=0.01,
            )

    def draw_hypergraph_sub_no_filter(self, name=None, sub=[], highlight=[]):
        if len(sub) == 0:
            hypx = self.nx
            sub = list(self.nodes.keys())
        else:
            sub = list(set(sub))
            hypx = self.nx.subgraph(sub)

        pos = nx.nx_agraph.graphviz_layout(hypx, prog="dot")

        propagation_values = [self.nodes[node].propagations for node in sub]
        norm = mcolors.Normalize(
            vmin=min(propagation_values) - 1, vmax=max(propagation_values)
        )
        # print(matplotlib.colormaps)
        cmap = plt.cm._colormaps.get_cmap("Blues")

        node_colors = [cmap(norm(value)) for value in propagation_values]

        node_size = 0.3

        node_colors_final = []
        s = []
        idxes = []
        for i, k in enumerate(sub):
            if k in highlight:
                node_colors_final.append("gold")
                s.append(node_size * 2)
                idxes.append(i)
                continue
            node_colors_final.append(node_colors[i])
            s.append(node_size)

        fig, ax = plt.subplots(dpi=300)
        fig.set_size_inches(2, 2)

        G = hypx
        print(len(idxes))
        nodes_x, nodes_y = zip(*[pos[n] for n in G.nodes()])
        print(len(nodes_x), len(nodes_y), len(node_colors_final), len(s))
        ax.scatter(
            nodes_x,
            nodes_y,
            s=s,
            color=node_colors_final,
            alpha=0.6,
            zorder=2,
        )

        # plot highlight overlay
        for i in idxes:
            ax.scatter(
                nodes_x[i],
                nodes_y[i],
                s=s[i],
                color="gold",
                alpha=1,
                zorder=4,
            )

        for edge in G.edges():
            start_x, start_y = pos[edge[0]]
            end_x, end_y = pos[edge[1]]

            # Calculate the direction vector of the edge
            edge_vec = np.array([end_x - start_x, end_y - start_y])
            edge_length = np.linalg.norm(edge_vec)
            edge_direction = edge_vec / edge_length

            # Shorten the edge to stop at the node boundary instead of the center
            start_x, start_y = (
                start_x + edge_direction[0] * node_size,
                start_y + edge_direction[1] * node_size,
            )
            end_x, end_y = (
                end_x - edge_direction[0] * node_size,
                end_y - edge_direction[1] * node_size,
            )

            if edge[0] in highlight and edge[1] in highlight:
                ax.plot(
                    [start_x, end_x],
                    [start_y, end_y],
                    color="gold",
                    alpha=1,
                    zorder=3,
                    linewidth=0.5,
                )
            else:

                ax.plot(
                    [start_x, end_x],
                    [start_y, end_y],
                    color="black",
                    alpha=0.1,
                    zorder=1,
                    linewidth=0.3,
                )

        plt.axis("off")

        if name == None:
            plt.show()
        else:
            plt.savefig(
                f"{name}.png",
                dpi=900,
                transparent=True,
                bbox_inches="tight",
                pad_inches=0.01,
            )

    def draw_hypergraph_single_color(
        self,
        sub,
        name=None,
        highlight=[],
        text=False,
        fig_size=(2, 2),
        node_size=0.3,
        cm="Blues",
        attr="prop",
        highlight_color="lightgrey",
        highlight_scale=1,
        try_rooting=True,
    ):
        # if len(sub) == 0:
        #     G = self.nx
        #     sub = list(self.nodes.keys())
        # else:
        #     sub = list(set(sub))
        #     G = self.nx.subgraph(sub).copy()
        G = sub.copy()

        # G.nodes[self.get_root_node()]["rank"] = "min"
        A = nx.drawing.nx_agraph.to_agraph(G)
        root = self.get_root_node()
        if try_rooting:
            sg = A.add_subgraph(root, name="cluster_root")
            sg.graph_attr["rank"] = "min"
            A.graph_attr["root"] = root
            A.graph_attr["rankdir"] = "TB"
        # pos = nx.nx_agraph.graphviz_layout(A, prog="dot")
        A.layout(prog="dot", args='-Grank="max"')
        pos = {
            node: (
                float(A.get_node(node).attr["pos"].split(",")[0]),
                float(A.get_node(node).attr["pos"].split(",")[1]),
            )
            for node in A.nodes()
        }
        propagation_values = [self.nodes[node].propagations for node in sub]
        # print(propagation_values[0])
        norm = mcolors.Normalize(
            vmin=min(propagation_values) - 1, vmax=max(propagation_values)
        )
        # print(matplotlib.colormaps)
        cmap = plt.cm._colormaps.get_cmap(cm)

        if attr == "prop":
            node_colors = [cmap(norm(value)) for value in propagation_values]
        else:
            rand_vals = [
                random.random() * max(propagation_values)
                for _ in range(len(propagation_values))
            ]
            node_colors = [cmap(norm(value)) for value in rand_vals]
        node_colors_final = []
        s = []
        idxes = []
        for i, k in enumerate(sub):
            if k in highlight:
                node_colors_final.append(highlight_color)
                s.append(node_size * highlight_scale)
                idxes.append(i)
                continue
            node_colors_final.append(node_colors[i])
            s.append(node_size)

        fig, ax = plt.subplots(dpi=300)
        fig.set_size_inches(fig_size[0], fig_size[1])

        nodes_x, nodes_y = zip(*[pos[n] for n in G.nodes()])
        # print(len(nodes_x), len(nodes_y), len(node_colors_final), len(s))
        ax.scatter(
            nodes_x,
            nodes_y,
            s=s,
            color=node_colors_final,
            alpha=1,
            zorder=2,
        )
        if text:
            for i in range(len(nodes_x)):
                ax.text(
                    nodes_x[i],
                    nodes_y[i],
                    str(propagation_values[i]),
                    fontsize=3,
                )

        # plot highlight overlay
        for i in idxes:
            ax.scatter(
                nodes_x[i],
                nodes_y[i],
                s=s[i],
                color=highlight_color,
                alpha=0.5,
                zorder=2,
            )

        for edge in G.edges():
            start_x, start_y = pos[edge[0]]
            end_x, end_y = pos[edge[1]]

            edge_vec = np.array([end_x - start_x, end_y - start_y])
            edge_length = np.linalg.norm(edge_vec)
            edge_direction = edge_vec / edge_length

            start_x, start_y = (
                start_x + edge_direction[0] * node_size,
                start_y + edge_direction[1] * node_size,
            )
            end_x, end_y = (
                end_x - edge_direction[0] * node_size,
                end_y - edge_direction[1] * node_size,
            )

            if edge[0] in highlight and edge[1] in highlight:
                ax.plot(
                    [start_x, end_x],
                    [start_y, end_y],
                    color=highlight_color,
                    alpha=1,
                    zorder=2,
                    linewidth=0.9,
                )
            else:

                ax.plot(
                    [start_x, end_x],
                    [start_y, end_y],
                    color="grey",
                    alpha=0.4,
                    zorder=1,
                    linewidth=0.3,
                )

        plt.axis("off")

        if name == None:
            plt.show()
        else:
            plt.savefig(
                f"{name}.png",
                dpi=900,
                transparent=True,
                bbox_inches="tight",
                pad_inches=0.01,
            )

    def __str__(self):
        return f"Nodes: {self.nodes}\nEdges: {self.edges}"

    def __repr__(self):
        return f"Nodes: {self.nodes}\nEdges: {self.edges}"


def returnReactionTemplates():
    enumerated_templates = {}
    for k in templates:
        for v in templates[k]["scope"]:
            for step in templates[k]["mechanism"]:
                temp = step["template"]
                for variable in v:
                    if variable == "name":
                        continue
                    temp = temp.replace(variable, v[variable])

                # the template, the name of the reaction,
                # the description of the specififed scope,
                # and a description of the transformation
                canon_mech = MechanisticReaction(
                    temp, k, v["name"], step["description"]
                )
                if canon_mech.canon_template not in enumerated_templates:
                    enumerated_templates[canon_mech.canon_template] = []
                enumerated_templates[canon_mech.canon_template].append(canon_mech)

    return enumerated_templates
