import flask
import webkit
from webkit.api.graph_util import generateSMARTSfromDrawer, getDrawing
from shared.templates.templates import templates
from shared.reaction_class import MechanisticReaction
import copy
from rdkit import Chem
from rdkit.Chem import AllChem
import re
import pandas as pd

random_smiles_dataset = pd.read_excel(
    webkit.app.config["DATA_FOLDER"] + "/random_25000_smiles.xlsx"
)


@webkit.app.route("/api/get_all_template_sets", methods=["GET"])
def fetch_all_name_mechanisms2():
    context = {}

    context["mechanisms"] = templates

    return flask.jsonify(**context)


@webkit.app.route("/api/run_test_input", methods=["POST"])
def run_test_mechanism():
    context = {}

    json_data = flask.request.get_json()

    if json_data["input"] == "":
        context["output"] = "no input entered"
        return flask.jsonify(**context)

    inputs = json_data["input"].split(".")

    temp = json_data["mechanism"]["template"]

    if json_data["scope"] == None:
        context["output"] = "no scope selected"
        return flask.jsonify(**context)
    else:
        for rr in json_data["scope"]:
            if rr == "name":
                continue
            temp = temp.replace(rr, json_data["scope"][rr])

        canon_mech = MechanisticReaction(
            temp, json_data["mechanism"]["description"], json_data["scope"]["name"], ""
        )

        reactants = []
        atom_map_label = 1
        for inp in inputs:
            mol = Chem.MolFromSmiles(inp, sanitize=False)
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE)
            if not mol:
                context["output"] = "invalid input"
                return flask.jsonify(**context)

            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(atom_map_label)
                atom_map_label += 1

            reactants.append(Chem.MolToSmiles(mol))

        reactant_objs = [Chem.MolFromSmiles(".".join(reactants), sanitize=False)]
        amap = 1
        for r in reactant_objs:
            for atom in r.GetAtoms():
                atom.SetAtomMapNum(amap)
                amap += 1

        products, reacting_atoms, _ = canon_mech.run(
            reactant_objs, intramolecular=json_data["intra"]
        )

        if not products:
            context["output"] = "number of inputs do not match number of templates"
            return flask.jsonify(**context)

        all_outs = []
        for prod in products:
            these_outs = []
            for p in prod:
                these_outs.append(p)
            all_outs.append(".".join(these_outs))
        context["output"] = ".".join(all_outs)

    if context["output"] == "":
        context["output"] = "no output"

    return flask.jsonify(**context)


def standardize_smiles_and_tautomer(input_smiles):

    try:
        prod_mol = Chem.MolFromSmiles(input_smiles, sanitize=False)

        for atom in prod_mol.GetAtoms():
            atom.SetAtomMapNum(0)

        prod_smiles = Chem.MolToSmiles(prod_mol, isomericSmiles=False)
    except:
        return None
    return prod_smiles


@webkit.app.route("/api/search_mechanisms", methods=["POST"])
def search_mechanisms():
    context = {}

    json_data = flask.request.get_json()

    if json_data["search"] == "":
        context["output"] = "no input entered"
        return flask.jsonify(**context)

    reactants = json_data["search"].split(">>")[0]

    reactants2 = []
    atom_map_label = 1
    for inp in reactants.split("."):
        mol = Chem.MolFromSmiles(inp, sanitize=False)
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE)
        if not mol:
            context["output"] = "invalid input"
            return flask.jsonify(**context)

        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(atom_map_label)
            atom_map_label += 1

        reactants2.append(Chem.MolToSmiles(mol))

    reactant_objs = [Chem.MolFromSmiles(".".join(reactants2), sanitize=False)]

    products = json_data["search"].split(">>")[1]

    prod_search = standardize_smiles_and_tautomer(products)

    mechs = json_data["mechanisms"]["mechanisms"]
    hits = []

    for temp in mechs:
        t = mechs[temp]

        for step in mechs[temp]["mechanism"]:
            for rr in mechs[temp]["scope"]:
                t_scoped = step["template"]
                for rr2 in rr:
                    if rr2 == "name":
                        continue
                    t_scoped = t_scoped.replace(rr2, rr[rr2])

                canon_mech = MechanisticReaction(
                    t_scoped, step["description"], rr["name"], ""
                )

                products, reacting_atoms, _ = canon_mech.run(
                    reactant_objs, intramolecular=True
                )
                if not products:
                    continue

                for p in products:
                    for pp in p:
                        for ppp in pp.split("."):
                            this_prod = standardize_smiles_and_tautomer(ppp)
                            if this_prod == prod_search:
                                print(f"Match found: {temp}")
                                dat = {
                                    "template": t_scoped,
                                    "description": temp + " " + step["description"],
                                    "scope": rr,
                                }
                                hits.append(dat)

    context["hits"] = hits
    return flask.jsonify(**context)


@webkit.app.route("/api/query_pattern", methods=["POST"])
def query_pattern():
    context = {}

    json_data = flask.request.get_json()
    temp = json_data["mechanism"]["template"]

    if json_data["scope"] == None:
        context["output"] = "no scope selected"
        return flask.jsonify(**context)
    else:
        for rr in json_data["scope"]:
            if rr == "name":
                continue
            temp = temp.replace(rr, json_data["scope"][rr])

        rxn = AllChem.ReactionFromSmarts(temp)

        reactants = []
        for inp in rxn.GetReactants():
            for i, rand_sm in random_smiles_dataset.iterrows():
                if "." in rand_sm["smiles"]:
                    continue
                rand_mol = Chem.MolFromSmiles(rand_sm["smiles"])
                if not rand_mol:
                    continue

                mol_wt = Chem.rdMolDescriptors.CalcExactMolWt(rand_mol)
                if mol_wt > 400:
                    continue
                if rand_mol.HasSubstructMatch(inp):
                    reactants.append(rand_sm["smiles"])
                    break

        context["output"] = ".".join(reactants)

    if context["output"] == "":
        context["output"] = "no hits found"

    return flask.jsonify(**context)


@webkit.app.route("/api/get_drawing_from_smiles", methods=["POST"])
def get_drawing_from_smiles():
    json_data = flask.request.get_json()
    context = {}

    try:
        out = getDrawing(json_data)
    except:
        context["error"] = True
        context["output"] = {"nodes": [], "edges": []}
        return flask.jsonify(**context)

    context["output"] = out
    context["error"] = False

    return flask.jsonify(**context)


@webkit.app.route("/api/get_smarts_from_drawing", methods=["POST"])
def get_smarts_from_drawing():
    json_data = flask.request.get_json()
    context = {}

    out = generateSMARTSfromDrawer(json_data)

    context["smarts"] = out

    return flask.jsonify(**context)


def interpolate_old_range_to_new_range(old_range, new_range, value):
    old_min, old_max = old_range
    new_min, new_max = new_range

    old_range = old_max - old_min
    new_range = new_max - new_min

    if old_range == 0:
        new_value = new_min
    else:
        new_value = (((value - old_min) * new_range) / old_range) + new_min

    return new_value


@webkit.app.route("/api/get_drawer_nodes_and_edges", methods=["POST"])
def get_mechanism_drawings():
    json_data = flask.request.get_json()
    context = {}

    context["nodes"] = []
    context["edges"] = []
    context["arrows"] = []
    for r in json_data["mechanism"]:

        if json_data["scope"] != None:
            for rr in json_data["scope"]:
                if rr == "name":
                    continue
                r["template"] = r["template"].replace(rr, json_data["scope"][rr])

        reacts = r["template"].split(">>")[0]
        if reacts[0] == "(":
            reacts = reacts[1:-1]

        reactants = reacts.split(".")

        all_nodes = {"nodes": [], "edges": []}
        max_x = 50
        new_id = 0
        node_to_id_shift = {}
        new_edges = []
        edge_id = 0
        for react in reactants:
            out = getDrawing({"smiles": react})

            mol = Chem.MolFromSmarts(react)
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE)
            for i, rr in enumerate(out["nodes"]):
                rr["x"] = rr["x"] * 30 + max_x
                rr["y"] = rr["y"] * 30 + 125

                atom = mol.GetAtomWithIdx(i)

                rr["label"] = atom.GetAtomMapNum()

                rr["atom"] = atom.GetSymbol()
                pattern = r"H-?\d+"
                matches = re.findall(pattern, atom.GetSmarts())
                rr["smarts"] = atom.GetSmarts()
                if len(matches) == 0:
                    rr["hydrogens"] = 0
                else:
                    rr["hydrogens"] = int(matches[0][1:])

                if "-:" in atom.GetSmarts():
                    rr["charge"] = -1
                elif "+:" in atom.GetSmarts():
                    rr["charge"] = 1
                else:
                    rr["charge"] = atom.GetFormalCharge()
                node_to_id_shift[rr["id"]] = new_id
                rr["id"] = new_id
                new_id += 1

                all_nodes["nodes"].append(rr)

            for rr in out["nodes"]:
                if rr["x"] > max_x:
                    max_x = rr["x"]

            max_x = max_x + 50
            all_nodes["edges"].extend(out["edges"])

            for rr in all_nodes["edges"]:
                _new_edge = copy.deepcopy(rr)
                _new_edge["id"] = edge_id
                edge_id += 1
                _new_edge["source"] = node_to_id_shift[rr["source"]]
                _new_edge["target"] = node_to_id_shift[rr["target"]]
                new_edges.append(_new_edge)

            all_nodes["edges"] = []
            node_to_id_shift = {}

        context["arrows"].append([{"x": max_x, "y": 125, "direction": "right"}])
        max_x += 150

        prods = r["template"].split(">>")[1]
        if prods[0] == "(":
            prods = prods[1:-1]
        products = prods.split(".")
        for prod in products:
            out = getDrawing({"smiles": prod})

            mol = Chem.MolFromSmarts(prod)
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE)
            for i, rr in enumerate(out["nodes"]):
                rr["x"] = rr["x"] * 30 + max_x
                rr["y"] = rr["y"] * 30 + 125
                atom = mol.GetAtomWithIdx(i)
                rr["atom"] = atom.GetSymbol()
                rr["label"] = atom.GetAtomMapNum()
                rr["smarts"] = atom.GetSmarts()

                pattern = r"H-?\d+"
                matches = re.findall(pattern, atom.GetSmarts())
                if len(matches) == 0:
                    rr["hydrogens"] = 0
                else:
                    rr["hydrogens"] = int(matches[0][1:])

                if "-:" in atom.GetSmarts():
                    rr["charge"] = -1
                elif "+:" in atom.GetSmarts():
                    rr["charge"] = 1
                else:
                    rr["charge"] = atom.GetFormalCharge()

                node_to_id_shift[rr["id"]] = new_id

                rr["id"] = new_id
                new_id += 1

                all_nodes["nodes"].append(rr)

            all_nodes["edges"].extend(out["edges"])
            for rr in out["nodes"]:
                if rr["x"] > max_x:
                    max_x = rr["x"]

            max_x += 50

            for rr in all_nodes["edges"]:
                _new_edge = copy.deepcopy(rr)
                _new_edge["id"] = edge_id
                edge_id += 1
                _new_edge["source"] = node_to_id_shift[rr["source"]]
                _new_edge["target"] = node_to_id_shift[rr["target"]]
                new_edges.append(_new_edge)

            all_nodes["edges"] = []
            node_to_id_shift = {}

        context["nodes"].append(all_nodes["nodes"])
        context["edges"].append(new_edges)

    return flask.jsonify(**context)
