from rdkit import Chem
from rdkit import Chem
from rdkit.Chem import AllChem
import networkx as nx
from shared.reaction_class import mark_if_ring_is_unacceptable


def mark_nodes_with_tcp(nodes, network):
    for node in nodes:
        mols_in_state = node.mapped_smiles.split(".")
        for r in mols_in_state:
            mol = Chem.MolFromSmiles(r, sanitize=False)

            atom_maps_in_mol = [atom.GetAtomMapNum() for atom in mol.GetAtoms()]
            found_starting_mols = []
            for inp in network.starting_material_atom_maps:
                for at in network.starting_material_atom_maps[inp]:
                    if at in atom_maps_in_mol:
                        found_starting_mols.append(inp)
                        break
            if len(set(found_starting_mols)) == len(
                network.starting_material_atom_maps
            ):
                node.other_data["tcp"] = True
                node.other_data["target_molecule"] = r
                break
            else:
                node.other_data["tcp"] = False


def mark_charged_states(nodes):
    for node in nodes:
        if not node.other_data["tcp"]:
            node.other_data["failed_thermo_state1"] = None
            continue
        bad_state = False
        mol = Chem.MolFromSmiles(node.other_data["target_molecule"], sanitize=False)
        num_charges = Chem.GetFormalCharge(mol)
        if num_charges != 0:
            bad_state = True
            node.other_data["failed_thermo_state1"] = bad_state
            continue
        for atom in mol.GetAtoms():
            if atom.GetFormalCharge() != 0:
                bad_state = True
                break

        node.other_data["failed_thermo_state1"] = bad_state


def mark_if_too_many_hets_in_simple_ring(smi):
    mol = Chem.MolFromSmiles(smi, sanitize=False)
    Chem.SanitizeMol(
        mol,
        sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE
        ^ Chem.SanitizeFlags.SANITIZE_SYMMRINGS
        ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY,
    )
    ri = mol.GetRingInfo()
    if len(ri.AtomRings()) == 0:
        return False
    for ring in ri.AtomRings():
        het_atoms = 0
        for atom in ring:
            if mol.GetAtomWithIdx(atom).GetIsAromatic():
                break
            if mol.GetAtomWithIdx(atom).GetAtomicNum() != 6:
                het_atoms += 1
        if het_atoms > 2:
            return True

    for ring in ri.AtomRings():
        het_atoms = 0
        for atom in ring:
            if not mol.GetAtomWithIdx(atom).GetIsAromatic():
                break
            if mol.GetAtomWithIdx(atom).GetAtomicNum() != 6:
                het_atoms += 1
        if het_atoms > 3:
            return True

    return False


def mark_non_participating_transformations(nodes):
    for node in nodes:
        if not node.other_data["tcp"]:
            node.other_data["contains_non_participating_transformation"] = None
            continue
        if node.other_data["failed_thermo_state1"] == True:
            node.other_data["contains_non_participating_transformation"] = None
            continue

        atom_network = nx.Graph()
        for a in node.these_reacting_atoms_path:
            for aa in a:
                atom_network.add_node(aa)
        for p in node.these_reacting_atoms_path:
            for a in p:
                for b in p:
                    if a == b:
                        continue
                    atom_network.add_edge(a, b)

        if nx.is_connected(atom_network):
            node.other_data["contains_non_participating_transformation"] = False
        else:
            node.other_data["contains_non_participating_transformation"] = True


patt1 = Chem.MolFromSmarts("[C](-[N,O,S])(-[N,O,S])-[N,O,S]")


def mark_structural_failure(nodes):
    for node in nodes:
        if not node.other_data["tcp"]:
            node.other_data["structural_failure"] = None
            continue
        if node.other_data["failed_thermo_state1"] == True:
            node.other_data["structural_failure"] = None
            continue
        if node.other_data["contains_non_participating_transformation"] == True:
            node.other_data["structural_failure"] = None
            continue

        if mark_if_too_many_hets_in_simple_ring(node.other_data["target_molecule"]):
            node.other_data["structural_failure"] = True
            continue

        mol = Chem.MolFromSmiles(node.other_data["target_molecule"], sanitize=False)
        if mol.HasSubstructMatch(patt1):
            node.other_data["structural_failure"] = True
            continue

        if mark_if_ring_is_unacceptable(mol):
            node.other_data["structural_failure"] = True
            continue
        node.other_data["structural_failure"] = False


def apply_filters_local(network):
    nodes = apply_filters([network[n] for n in network.nodes], network, True)
    for n in nodes:
        network[n.unmapped_smiles].tcp = n.other_data["tcp"]
        network[n.unmapped_smiles].failed_thermo_state1 = n.other_data[
            "failed_thermo_state1"
        ]
        network[n.unmapped_smiles].contains_non_participating_transformation = (
            n.other_data["contains_non_participating_transformation"]
        )
        network[n.unmapped_smiles].structural_failure = n.other_data[
            "structural_failure"
        ]


def apply_filters(nodes, network, verbose=False):

    mark_nodes_with_tcp(nodes, network)
    if verbose:
        print("tcp", len([n for n in nodes if n.other_data["tcp"] == True]))

    mark_charged_states(nodes)
    if verbose:
        print(
            "ions",
            len([n for n in nodes if n.other_data["failed_thermo_state1"] == False]),
        )

    mark_non_participating_transformations(nodes)
    if verbose:
        print(
            "non participating atoms",
            len(
                [
                    n
                    for n in nodes
                    if n.other_data["contains_non_participating_transformation"]
                    == False
                ]
            ),
        )

    mark_structural_failure(nodes)
    if verbose:
        print(
            "structural failures",
            len([n for n in nodes if n.other_data["structural_failure"] == False]),
        )
    return nodes


def print_filter_counts(nodes):
    tcps = len([n for n in nodes if n.other_data["tcp"] == False])
    charged = len([n for n in nodes if n.other_data["failed_thermo_state1"] == True])
    non_part = len(
        [
            n
            for n in nodes
            if n.other_data["contains_non_participating_transformation"] == True
        ]
    )
    struct = len([n for n in nodes if n.other_data["structural_failure"] == True])

    print("initial nodes", len(nodes))
    print("failed tcp:", tcps)
    print("failed ion:", charged)
    print("failed non participating:", non_part)
    print("failed str:", struct)
    print(
        "passing:",
        len([n for n in nodes if n.other_data["structural_failure"] == False]),
    )
    for n in nodes:
        if n.other_data["structural_failure"] == False:
            print(n.product_in_precalc)
            print(n.other_data["target_molecule"])
            print(n)


def get_drugbank_fps(loc="./data/drugbank_all_structures_20231226.sdf"):
    drugbank = Chem.SDMolSupplier(loc)

    drugbank_fps = []
    for mol in drugbank:
        if mol == None:
            continue
        drugbank_fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, 1024))
    return drugbank_fps
