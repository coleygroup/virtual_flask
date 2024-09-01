from shared.reaction_class import returnReactionTemplates, VirtualFlask
from shared.util import get_path_from_init_to_node
import dill
import sys
from rdkit import RDLogger, Chem
import os
from shared.filters3 import (
    apply_filters,
    print_filter_counts,
    get_drugbank_fps,
)
from shared.hpc.energy_util import calc
import psycopg2
from psycopg2 import sql
from psycopg2.extras import Json
import time
import requests
import itertools
import json
from rdkit.Chem import rdFingerprintGenerator, DataStructs, AllChem
from psycopg2.extras import execute_batch
import re
import networkx as nx
import random

RDLogger.DisableLog("rdApp.*")

table_name_global = "test_table5"


class Node:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __repr__(self):
        return f"Node({', '.join(f'{k}={v}' for k, v in self.__dict__.items())})"

    def __getattr__(self, key):
        return self.__dict__.get(key, None)

    def __setattr__(self, key, value):
        self.__dict__[key] = value


class Network:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __repr__(self):
        return f"Network({', '.join(f'{k}={v}' for k, v in self.__dict__.items())})"

    def __getattr__(self, key):
        return self.__dict__.get(key, None)

    def __setattr__(self, key, value):
        self.__dict__[key] = value

    def load_data(self, nodes, edges):
        self.nodes = nodes
        self.edges = edges

        G = nx.DiGraph()

        nx_nodes = []
        for k in nodes:
            nx_nodes.append((k.node_id, k))

        nx_edges = []
        for k in edges:
            nx_edges.append((k.source_node_id, k.destination_node_id, k.serialize()))

        G.add_nodes_from(nx_nodes)
        G.add_edges_from(nx_edges)

        self.nx_graph = G

        self.node_map = {}
        for i in nodes:
            self.node_map[i.node_id] = i

    def get_shortest_path(self, node1, node2):
        return nx.shortest_path(self.nx_graph, node1, node2)

    def get_edge(self, node1, node2):
        # print("n", node1, node2)
        return self.nx_graph.get_edge_data(node1, node2)

    def get_edges_of_path(self, node_path):
        edges = []
        for i in range(len(node_path) - 1):
            edge = self.get_edge(node_path[i], node_path[i + 1])
            # print(edge)
            edges.append(edge)

        return edges

    def get_origin_node(self):
        for node in self.nodes:
            if node.these_reacting_atoms_path == None:
                return node.node_id

        return None

    def draw_graph(self):
        nx.draw(self.nx_graph, with_labels=False)


class Edge:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __repr__(self):
        return f"Edge({', '.join(f'{k}={v}' for k, v in self.__dict__.items())})"

    def __getattr__(self, key):
        return self.__dict__.get(key, None)

    def __setattr__(self, key, value):
        self.__dict__[key] = value

    def serialize(self):
        return {
            "edge_id": self.edge_id,
            "network_id": self.network_id,
            "source_node_id": self.source_node_id,
            "destination_node_id": self.destination_node_id,
            "other_data": self.other_data,
        }


def standardize_smiles_and_tautomer(input_smiles):

    try:
        prod_mol = Chem.MolFromSmiles(input_smiles)

        prod_inchi = Chem.MolToInchi(prod_mol)
        prod_mol = Chem.MolFromInchi(prod_inchi)
        prod_smiles = Chem.MolToSmiles(prod_mol, isomericSmiles=False)
    except:
        return None
    return prod_smiles


def get_top_n_products_from_askcos(inputs, n=30):
    json_data = {
        "model_name": "pistachio_23Q3",
        "smiles": inputs,
        "reagents": ["O=C(O)C1=CC=CC=C1"],
    }

    headers = {}
    response = requests.post(  # 18.224.138.75
        "http://3.138.199.28:9510/predictions/pistachio_23Q3",
        headers=headers,
        json=json_data,
    )

    if "code" in response.json():
        if response.json()["code"] == 500 or response.json()["code"] == 503:
            print(response.json())
            print(json_data)
            return None, None

    out_products = []
    out_scores = []
    for idx, i in enumerate(response.json()):
        try:
            products, scores = [], []
            for idx2, x in enumerate(i["products"]):
                o = standardize_smiles_and_tautomer(x)
                if not o:
                    continue
                scores.append(i["scores"][idx2])
                products.append(o)
            out_products.append(products)
            out_scores.append(scores)
        except:
            out_products.append([])
            out_scores.append([])

    return out_products[0], out_scores[0]


def calculate_morgan_fp(smiles, radius=2, n_bits=1024):
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    Chem.SanitizeMol(
        mol,
        sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE
        ^ Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
    )

    if mol is None:
        return None
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    fp = gen.GetFingerprint(mol)
    # Convert the fingerprint to a bit string for storing in the database
    return fp, "".join([str(b) for b in fp])


def connect_to_rds(
    host="localhost",
    port=6432,
    dbname="testdb",
    user="myuser",
    password="mypassword",
):
    """Create a connection to the PostgreSQL database"""
    try:
        conn = psycopg2.connect(
            host=host, port=port, dbname=dbname, user=user, password=password
        )
        print("Connection established")
        return conn
    except Exception as e:
        print(f"Unable to connect to the database: {e}")
        print("No connection to RDS")
        sys.exit(1)


def upload_data(conn, table_name, data):
    """Upload data to the specified table"""
    try:
        with conn.cursor() as cur:
            columns = [column for column in data]

            # data["path_data"] = Json(data["path_data"])
            # print(data["path_data"])
            data["path_data"] = json.dumps(data["path_data"])
            # print(json.dumps(data["path_data"]))
            # print(data["path_data"])
            values = [data[column] for column in columns]
            placeholders = sql.SQL(", ").join([sql.Placeholder() for _ in values])
            query = sql.SQL("INSERT INTO {table} ({fields}) VALUES ({values})").format(
                table=sql.Identifier(table_name),
                fields=sql.SQL(", ").join(map(sql.Identifier, columns)),
                values=placeholders,
            )

            # print(query)
            # print(values)
            cur.execute(query, values)
            conn.commit()
            print("Data uploaded successfully")
    except Exception as e:
        print(f"Failed to upload data: {e}")
        conn.rollback()


def update_database(id_list, product_list, hits_3cm):
    conn = connect_to_rds()
    cursor = conn.cursor()
    for row_id, prod in zip(id_list, product_list):
        cursor.execute(
            "UPDATE %s SET found_in_askcos_forward = %s WHERE id = %s;",
            (table_name_global, hits_3cm[prod], row_id),
        )

    conn.commit()
    cursor.close()
    conn.close()
    print("Connection closed")


def map_1(dir, data, index_value):

    mechs = returnReactionTemplates()
    hyps = []
    for i, k in enumerate(data):
        print("running", k)

        conn = connect_to_rds()

        cursor = conn.cursor()
        cursor.execute(
            "SELECT askcos_product_array FROM inputs_to_askcos_products WHERE inputs = %s;",
            (k,),
        )

        found_hit = cursor.fetchone()

        if not found_hit:
            print("skipping (no precalc)", k)
            continue

        precalc_prods = found_hit[0]
        print(precalc_prods)
        print(precalc_prods[0])
        print(len(precalc_prods))

        cursor.close()
        conn.close()
        # precalc_prods, precalc_scores = precalculate_novelty_askcos(k.split("."))

        state_network = VirtualFlask(mechs)
        state_network.charge(k.split("."), [])
        state_network.run_until_done(
            iters=5, thresh=50000, ring_filter=False, precalc_prods=precalc_prods
        )

        conn = connect_to_rds()
        hyp_id = state_network.push_to_db(conn)
        conn.close()
        print("conn closed, hyp_id:", hyp_id)

        hyps.append(str(hyp_id))

        print("finished", k)
        print()

    with open(f"{dir}/map_2_input_data/input_{index_value}.txt", "w") as f:
        f.write("\n".join(hyps))

    # return state_network


def get_network(cursor, network_id):
    cursor.execute(
        "SELECT * FROM test_networks WHERE network_id = %s;",
        (network_id,),
    )

    network = cursor.fetchone()

    colnames = [desc[0] for desc in cursor.description]
    row_dict = dict(zip(colnames, network))
    network = Network(**row_dict)
    return network


def get_nodes(cursor, network_id, query=None, values="*"):
    if query == None:
        cursor.execute(
            "SELECT " + values + " FROM test_nodes WHERE network_id = %s;",
            (network_id,),
        )
    else:
        cursor.execute(
            "SELECT "
            + values
            + " FROM test_nodes WHERE network_id = %s AND "
            + query
            + ";",
            (network_id,),
        )

    onodes = cursor.fetchall()

    colnames = [desc[0] for desc in cursor.description]
    nodes = []
    for row in onodes:
        row_dict = dict(zip(colnames, row))
        test_node = Node(**row_dict)
        nodes.append(test_node)

    return nodes


def get_edges(cursor, network_id):
    cursor.execute(
        "SELECT * FROM test_edges WHERE network_id = %s;",
        (network_id,),
    )

    oedges = cursor.fetchall()

    colnames = [desc[0] for desc in cursor.description]
    edges = []
    for row in oedges:
        row_dict = dict(zip(colnames, row))
        test_edge = Edge(**row_dict)
        edges.append(test_edge)

    return edges


def map_2(dir, data, index_value):
    files = []
    for i in data:
        conn = connect_to_rds()
        cursor = conn.cursor()

        nodes = get_nodes(cursor, i)
        network = get_network(cursor, i)

        cursor.close()
        conn.close()

        apply_filters(nodes, network)
        print_filter_counts(nodes)

        updates = []
        for node in nodes:
            updates.append((json.dumps(node.other_data), node.node_id))

        conn = connect_to_rds()
        cursor = conn.cursor()
        execute_batch(
            cursor,
            "UPDATE test_nodes SET other_data = %s WHERE node_id = %s;",
            updates,
        )

        conn.commit()
        cursor.close()
        conn.close()

        files.append(i)
        print("finished", i)
        print()

    with open(f"{dir}/thermo_input_data/input_{index_value}.txt", "w") as f:
        f.write("\n".join(files))


def assign_energy(sm_ar, idx, precalc):
    node_energy = 0
    for idr, i in enumerate(sm_ar):
        if i in precalc:
            ene = precalc[i]
        else:
            ene = calc(i, str(idx) + "-" + str(idr))
            precalc[i] = ene
        node_energy += ene

    return node_energy


def thermo(dir, data, index_value):
    files = []
    precalc = {}

    conn = connect_to_rds()
    cursor = conn.cursor()

    cursor.execute("SELECT smiles, score FROM xtb_calc;")
    orig = []
    for row in cursor.fetchall():
        orig.append(row[0])
        precalc[row[0]] = row[1]

    cursor.close()
    conn.close()

    for i in data:
        conn = connect_to_rds()
        cursor = conn.cursor()

        nodes = get_nodes(cursor, i)
        edges = get_edges(cursor, i)
        network = get_network(cursor, i)

        cursor.close()
        conn.close()

        network.load_data(nodes, edges)

        origin_node_id = network.get_origin_node()
        print("network", i, "origin", origin_node_id)

        input_smiles = ".".join(network.mapped_input_smiles_list)
        input_mol = Chem.MolFromSmiles(input_smiles)

        total_energy_input = assign_energy(
            network.input_unmapped_smiles.split("."), origin_node_id, precalc
        )

        mol_formula = Chem.rdMolDescriptors.CalcMolFormula(input_mol)
        # pattern = r"H\d+(?:,H\d+)*"
        # hydrogens = re.search(pattern, mol_formula)

        # num_input_hs = int(hydrogens[0][1:])
        print("inp", total_energy_input, mol_formula, input_smiles)

        for n in nodes:
            if n.other_data["structural_failure"] == False:

                path = network.get_shortest_path(origin_node_id, n.node_id)

                energy_path = []

                path_failed = False
                for p in path:
                    nn = network.node_map[p]
                    if "energy" in nn.other_data:
                        if nn.other_data["energy"] != None:
                            print(
                                "prec",
                                p,
                                nn.other_data["energy"],
                                network.node_map[p].unmapped_smiles,
                            )
                            energy_path.append(
                                (nn.other_data["energy"], nn.unmapped_smiles)
                            )
                            continue
                    nn_sms = nn.unmapped_smiles.split(".")
                    nn_node_energy = assign_energy(nn_sms, nn.node_id, precalc)
                    if nn_node_energy < -90000:
                        path_failed = True
                    nn.other_data["energy"] = nn_node_energy
                    energy_path.append((nn_node_energy, nn.unmapped_smiles))
                    print(p, nn_node_energy, network.node_map[p].unmapped_smiles)

                # mol = Chem.MolFromSmiles(n.unmapped_smiles)
                # mol_formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
                # print(mol_formula, n.unmapped_smiles)
                # prod_hydrogens = re.search(pattern, mol_formula)

                # num_prod_hs = int(prod_hydrogens[0][1:])

                # if num_prod_hs > num_input_hs:
                #     print("wtf?")
                # elif num_prod_hs < num_input_hs:
                #     for i in range(num_input_hs - num_prod_hs):
                #         prod_sms.append("[H]")
                # prod_sms = n.unmapped_smiles.split(".")

                # node_energy = assign_energy(prod_sms, n.node_id, precalc)

                # if node_energy < -90000:
                # path_failed = True

                # n.other_data["energy"] = node_energy
                # energy_path.append((node_energy, n.unmapped_smiles))
                n.other_data["energy_path"] = energy_path
                n.other_data["energy_failed"] = path_failed
                # print(node_energy, energy_path)

            else:
                n.other_data["energy"] = None
                n.other_data["energy_path"] = None
                n.other_data["energy_failed"] = None

        updates = []
        for node in nodes:
            updates.append((json.dumps(node.other_data), node.node_id))

        conn = connect_to_rds()
        cursor = conn.cursor()

        execute_batch(
            cursor,
            "UPDATE test_nodes SET other_data = %s WHERE node_id = %s;",
            updates,
        )

        calc_updates = []
        for k in precalc:
            if k not in orig:
                calc_updates.append((k, precalc[k]))
                orig.append(k)

        execute_batch(
            cursor,
            "INSERT INTO xtb_calc (smiles, score) VALUES (%s, %s);",
            calc_updates,
        )

        conn.commit()
        cursor.close()
        conn.close()

        files.append(i)
        print("finished xtb", i)
        # print()

    with open(f"{dir}/dock_input_data/input_{index_value}.txt", "w") as f:
        f.write("\n".join(files))


def extract_docking_score(out_file):
    best_docking_score = None
    try:
        with open(out_file, "r") as f:
            for line in f:
                if line.strip().startswith("REMARK VINA RESULT:"):
                    parts = line.split()
                    if len(parts) > 3:
                        try:
                            docking_score = float(parts[3])
                            if (
                                best_docking_score is None
                                or docking_score < best_docking_score
                            ):
                                best_docking_score = docking_score
                        except ValueError:
                            pass
    except FileNotFoundError:
        print(f"Out file {out_file} not found.")
    return best_docking_score


def dock(dir, data, index_value):
    files = []
    for i in data:
        conn = connect_to_rds()
        cursor = conn.cursor()

        nodes = get_nodes(cursor, i)
        edges = get_edges(cursor, i)
        network = get_network(cursor, i)

        cursor.close()
        conn.close()

        network.load_data(nodes, edges)

        origin_node_id = network.get_origin_node()
        print("network", i, "origin", origin_node_id)

        for n in nodes:
            if n.other_data["energy_failed"] == None:
                continue
            if n.other_data["energy_failed"] == False:
                ligand_filename = f"target_ligand_{str(n.node_id)}"
                sdf_writer = Chem.SDWriter("sdf/" + ligand_filename + ".sdf")

                mol = Chem.MolFromSmiles(n.other_data["target_molecule"])
                print("docking", n.other_data["target_molecule"])
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol)
                AllChem.MMFFOptimizeMolecule(mol)

                sdf_writer.write(mol)
                sdf_writer.close()

                os.system(
                    f"obabel -isdf sdf/{ligand_filename}.sdf -osdf -O sdf/{ligand_filename}_protonated.sdf -p 7"
                )
                os.system(
                    f"obabel -isdf sdf/{ligand_filename}_protonated.sdf -opdbqt -O pdbqt/{ligand_filename}.pdbqt -m"
                )
                os.system(
                    f"/home/bmahjour/autodock_vina_1_1_2_linux_x86/bin/vina --config data/file.conf --ligand pdbqt/{ligand_filename}1.pdbqt --out dock_out/{ligand_filename}.out --log dock_out/{ligand_filename}.log --cpu 4"
                )

                docking_score = extract_docking_score(f"dock_out/{ligand_filename}.out")
                print(docking_score)
                n.other_data["docking_score"] = docking_score
            else:
                n.other_data["docking_score"] = None

        updates = []
        for node in nodes:
            updates.append((json.dumps(node.other_data), node.node_id))

        conn = connect_to_rds()
        cursor = conn.cursor()

        execute_batch(
            cursor,
            "UPDATE test_nodes SET other_data = %s WHERE node_id = %s;",
            updates,
        )

        conn.commit()
        cursor.close()
        conn.close()

        files.append(i)
        print("finished docking", i)
        print()


def average_difference(values):
    if len(values) < 2:
        raise ValueError(
            "At least two elements are required to calculate average difference."
        )

    differences = [values[i + 1] - values[i] for i in range(len(values) - 1)]
    avg_diff = sum(differences) / len(differences)

    return avg_diff


def reduce(dir, data, index_value):
    conn = connect_to_rds()
    cur = conn.cursor()

    for i, k in enumerate(data):
        cur.execute(f"SELECT * FROM test_nodes WHERE node_id = {k};")
        k = int(k)
        colnames = [desc[0] for desc in cur.description]
        out = cur.fetchone()
        row_dict = dict(zip(colnames, out))
        test_node = Node(**row_dict)

        mol = Chem.MolFromSmiles(test_node.other_data["target_molecule"])
        mol_weight = Chem.Descriptors.ExactMolWt(mol)
        if mol_weight < 90:
            continue

        network_id = test_node.network_id

        net = get_network(cur, network_id)
        print("net", network_id, i, len(data))
        edges = get_edges(cur, network_id)
        print("e", len(edges))
        nodes = get_nodes(cur, network_id)
        print("n", len(nodes))
        net.load_data(nodes, edges)

        origin_node_id = net.get_origin_node()

        node_ids = [n.node_id for n in nodes]
        # print(node_ids)
        print(origin_node_id, origin_node_id in node_ids)
        print(k, k in node_ids)

        print()
        path = net.get_shortest_path(origin_node_id, k)

        ep = []
        energies = []
        for kk in path:

            cur.execute(f"SELECT * FROM test_nodes WHERE node_id = {kk};")

            colnames = [desc[0] for desc in cur.description]
            out = cur.fetchone()
            row_dict = dict(zip(colnames, out))
            node_kk = Node(**row_dict)

            # print(node_kk.other_data)
            ep.append((node_kk.other_data["energy"], node_kk.mapped_smiles))
            energies.append(node_kk.other_data["energy"])

        try:
            avg_diff = average_difference(energies)
        except:
            print(k, "Failed", energies)
            continue

        edge_path = net.get_edges_of_path(path)
        path_data = []
        mech_path = []
        for ed in edge_path:
            path_data.append(ed["other_data"]["template_obj"])
            mech_path.append(ed["other_data"]["template"])

        dat = {
            "unmapped_input_smiles": net.input_unmapped_smiles,
            "mapped_input_smiles": net.input_mapped_smiles,
            "propagations": len(ep),
            "mech_path": mech_path,
            "path_data": path_data,
            "unmapped_product_smiles": test_node.unmapped_smiles,
            "target_molecule": test_node.other_data["target_molecule"],
            "overall_reaction": test_node.unmapped_smiles
            + ">>"
            + test_node.other_data["target_molecule"],
            "mapped_product_smiles": test_node.mapped_smiles,
            "network_size": len(nodes),
            "average_delta_g": avg_diff,
            "docking_score": test_node.other_data["docking_score"],
            "energies": energies,
            "mechanistic_path": [e[1] for e in ep],
        }

        upload_data(conn, table_name_global, dat)


def reduce_1(dir, data):

    conn = connect_to_rds()
    db_fps = get_drugbank_fps()

    if not conn:
        print("No connection to RDS")
        sys.exit(1)

    for i in data:
        filename = f"{dir}/map_1_output_data/{i}"
        hyp = dill.load(open(filename, "rb"))

        sm = hyp.mapped_input_smiles_list
        nn = 0
        for i in hyp.nodes:
            if not hyp[i].tcp:
                continue
            if hyp[i].failed_thermo_state1:
                continue
            if hyp[i].known_product:
                continue
            if hyp[i].in_same_mechanism:
                continue
            if hyp[i].contains_non_participating_transformation:
                continue
            if hyp[i].structural_failure:
                continue

            path = get_path_from_init_to_node(hyp, i)
            this_path = []
            this_mech_path = []
            this_path_data = []
            for ip, p in enumerate(path):
                if ip + 1 == len(path):
                    break
                edge = hyp.get_edge(path[ip], path[ip + 1])[0]
                this_path.append(edge["template"])
                this_path_data.append(edge["template_obj"].serialize())
                this_mech_path.append(p)

            out = "~".join(this_path)

            unmapped_smiles = hyp[i].unmapped_smiles

            fpr, fp = calculate_morgan_fp(hyp[i].target_molecule)
            dists = []
            for d in db_fps:
                dist = 1 - DataStructs.cDataStructs.TanimotoSimilarity(fpr, d)
                dists.append(dist)
            min_dist = min(dists)

            reduced_data = {
                "unmapped_input_smiles": hyp.input_unmapped_smiles,
                "mapped_input_smiles": sm,
                "propagations": hyp[i].propagations,
                "path": out,
                "mech_path": this_mech_path,
                "path_data": this_path_data,
                "unmapped_product_smiles": unmapped_smiles,
                "target_molecule": hyp[i].target_molecule,
                "overall_reaction": ".".join(sm) + ">>" + hyp[i].target_molecule,
                "mapped_product_smiles": hyp[i].mapped_smiles,
                "network_size": len(hyp.nodes),
                "min_dist_to_db": min_dist,
                "morgan_fp": fp,
                # "mol": Chem.MolFromSmiles(hyp[i].target_molecule),
            }

            upload_data(conn, table_name_global, reduced_data)
            nn = nn + 1


def precalculate_novelty_askcos(inputs):
    time_00 = time.time()
    all_2mers = list(itertools.combinations(inputs, 2))

    layer_2_third_component = []
    for pair in all_2mers:
        remaining_component = [x for x in inputs if x not in pair]
        layer_2_third_component.append(remaining_component[0])

    prods, scores = get_top_n_products_from_askcos([".".join(inputs)], n=5)
    print(prods)
    # print("finished mcr", len(prods), time.time() - time_00)

    for idx, i in enumerate(all_2mers):
        prods_layer_1, scores_l1 = get_top_n_products_from_askcos([".".join(i)], n=5)
        for idx2, j in enumerate(prods_layer_1):
            prods_layer_2, scores_l2 = get_top_n_products_from_askcos(
                [j + "." + layer_2_third_component[idx]], n=5
            )

            scores_out = [
                (scores_l1[idx2], scores_l2[idx3]) for idx3 in range(len(scores_l2))
            ]
            # print(
            #     "finished layered", idx, idx2, len(prods_layer_2), time.time() - time_00
            # )
            prods.extend(prods_layer_2)
            scores.extend(scores_out)

    print("finished precalc askcos all", len(prods), time.time() - time_00)
    return prods, scores
