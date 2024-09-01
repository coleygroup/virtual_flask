import flask
import webkit
import dill
import pandas as pd
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, Descriptors, Draw
import numpy as np
import re
import copy
from rdkit.Chem.Draw import rdMolDraw2D
import base64
import io
from PIL import Image
import psycopg2
import json
from shared.hpc.commands import connect_to_rds


def getReactionImage(reaction_smarts):
    if reaction_smarts is None:
        return None

    # Convert reaction SMARTS to a reaction object
    # print(reaction_smarts)
    rxn = AllChem.ReactionFromSmarts(reaction_smarts, useSmiles=True)

    if not rxn:
        return None

    # Create a depiction of the reaction
    for r in rxn.GetReactants():
        Chem.SanitizeMol(r, sanitizeOps=Chem.SanitizeFlags.SANITIZE_CLEANUP)
        for a in r.GetAtoms():
            a.SetAtomMapNum(0)
        # Draw.rdMolDraw2D.PrepareMolForDrawing(r)
    for p in rxn.GetProducts():
        Chem.SanitizeMol(p, sanitizeOps=Chem.SanitizeFlags.SANITIZE_CLEANUP)
        for a in p.GetAtoms():
            a.SetAtomMapNum(0)
        # Draw.rdMolDraw2D.PrepareMolForDrawing(p)

    drawer = rdMolDraw2D.MolDraw2DCairo(-1, -1)
    dopts = drawer.drawOptions()
    dopts.bondLineWidth = 1.5  # default is 2.0
    drawer.DrawReaction(rxn)
    drawer.FinishDrawing()

    img = drawer.GetDrawingText()
    base64_str = base64.b64encode(img).decode("utf-8")

    return base64_str


def getImage(row, removeAtomMap = False):
    if row == None:
        return None

    mol = Chem.MolFromSmiles(row, sanitize=False)

    Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE)

    if not mol:
        return None

    if removeAtomMap:
        for a in mol.GetAtoms():
            a.SetAtomMapNum(0)

    mc = Chem.Mol(mol.ToBinary())
    AllChem.Compute2DCoords(mc)

    coords = mc.GetConformer(-1).GetPositions()

    # if the margin +/-1 is removed linear molecules are broken
    min_p = Geometry.Point2D(*coords.min(0)[:2] - 1)
    max_p = Geometry.Point2D(*coords.max(0)[:2] + 1)

    dpa = 20  # dots per angstrom
    # try to catch 0's with max()
    w = int(dpa * (max_p.x - min_p.x)) + 1
    h = int(dpa * (max_p.y - min_p.y)) + 1

    # drawer = rdMolDraw2D.MolDraw2DSVG(max(w, dpa), max(h, dpa))
    drawer = rdMolDraw2D.MolDraw2DCairo(max(w, dpa), max(h, dpa))
    dopts = drawer.drawOptions()
    dopts.bondLineWidth = 3

    drawer.DrawMolecule(mc)

    drawer.FinishDrawing()
    img = drawer.GetDrawingText()
    png_bytes_io = io.BytesIO(img)
    img = Image.open(png_bytes_io)
    img = img.convert("RGBA")
    datas = img.getdata()

    # Replace white pixels with transparent
    newData = []
    for item in datas:
        if item[0] == 255 and item[1] == 255 and item[2] == 255:  # Finding white pixels
            newData.append((255, 255, 255, 0))  # Changing white pixels to transparent
        else:
            newData.append(item)
    # Update image data with transparency
    img.putdata(newData)

    output_buffer = io.BytesIO()
    img.save(
        output_buffer, format="PNG", dpi=(300, 300)
    )  # You can change 'PNG' to 'JPEG'
    byte_data = output_buffer.getvalue()
    base64_str = base64.b64encode(byte_data).decode("utf-8")
    return base64_str

    # except:
    #     return None


# @webkit.app.route("/api/update_list", methods=["POST"])
# def update_list():
#     new_list = []
#     json_data = flask.request.get_json()

#     smiles = json_data["smiles"]

#     if smiles == "" or len(smiles) == 0:
#         for idx, h in hits.iterrows():
#             for k in h["unmapped_smiles"].split("."):
#                 if k not in new_list:
#                     new_list.append(k)

#         return flask.jsonify(list(set(new_list)))

#     for i in smiles:
#         new_list.append(i)

#     for idx, h in hits.iterrows():
#         all_smiles_found = True
#         for s in smiles:
#             if s not in h["unmapped_smiles"].split("."):
#                 all_smiles_found = False
#                 break
#         if all_smiles_found:
#             for k in h["unmapped_smiles"].split("."):
#                 if k not in new_list:
#                     new_list.append(k)

#     # print(len(new_list))

#     return flask.jsonify(list(set(new_list)))


@webkit.app.route("/api/get_product_masses", methods=["POST"])
def get_product_masses():
    json_data = flask.request.get_json()
    context = {}

    masses = []
    for r in json_data["smiles"]:
        mol = Chem.MolFromSmiles(r)
        masses.append(Chem.Descriptors.ExactMolWt(mol))

    context["masses"] = masses
    return flask.jsonify(**context)


@webkit.app.route("/api/get_buffer", methods=["POST"])
def get_buffer():
    json_data = flask.request.get_json()
    context = {}

    if type(json_data["buffer"][0]) == int:
        buffer = json_data["buffer"]
    else:
        buffer = [x[0] for x in json_data["buffer"]]

    # print(buffer)

    conn = None
    try:
        # Connect to your PostgreSQL database
        conn = connect_to_rds()
        cur = conn.cursor()

        # SQL query to fetch details based on bucket1 and bucket2
        query = """
        SELECT *
        FROM test_table5
        WHERE id = ANY(%s);
        """
        # Execute the query with parameters
        cur.execute(query, (buffer,))

        # Fetch all results
        records = cur.fetchall()

        column_headers = [desc[0] for desc in cur.description]
        column_headers.append("reaction_img")

        col_idx = column_headers.index("overall_reaction")
        # ask_idx = column_headers.index("found_in_askcos_forward")
        for idx, r in enumerate(records):
            records[idx] = list(r)

            overall_img = r[1] + ">>" + r[7]
            records[idx].append(getReactionImage(overall_img))



        # print(column_headers)
        context["headers"] = column_headers
        context["hits"] = records
        cur.close()
    except psycopg2.DatabaseError as e:
        print(f"Database error: {e}")
    finally:
        if conn is not None:
            conn.close()

    return flask.jsonify(**context)


@webkit.app.route("/api/get_bucket_data", methods=["POST"])
def get_bucket_data():

    context = {}
    bucket1 = flask.request.get_json()["bucket1"]
    bucket2 = flask.request.get_json()["bucket2"]
    ids = flask.request.get_json()["ids"]
    context["hit_ids"] = None
    """Fetches details from the histogram_details table for specific buckets."""
    conn = None
    try:
        # Connect to your PostgreSQL database
        conn = connect_to_rds()
        cur = conn.cursor()

        # SQL query to fetch details based on bucket1 and bucket2
        if len(ids) > 0:
            query = """
            SELECT id
            FROM histogram_details
            WHERE bucket1 = %s AND bucket2 = %s AND id = ANY(%s);
            """
            # Execute the query with parameters
            cur.execute(query, (bucket1, bucket2, ids))

        else:
            query = """
            SELECT id
            FROM histogram_details
            WHERE bucket1 = %s AND bucket2 = %s;
            """
            # Execute the query with parameters
            cur.execute(query, (bucket1, bucket2))

        # Fetch all results
        records = cur.fetchall()
        context["hit_ids"] = records
        # print(f"Data retrieved for Bucket1: {bucket1}, Bucket2: {bucket2}")
        # for record in records:
        # print(record)

        cur.close()
    except psycopg2.DatabaseError as e:
        print(f"Database error: {e}")
    finally:
        if conn is not None:
            conn.close()

    return flask.jsonify(**context)


@webkit.app.route("/api/get_hist", methods=["POST"])
def get_hist():

    d = flask.request.get_json()["selectedDisplay"]
    if d == "min_delta_lambda_max":
        query = "SELECT * FROM histogram_data3;"
    elif d == "network_size":
        query = "SELECT * FROM histogram_data2;"
    conn = connect_to_rds()

    # query = "SELECT * FROM histogram_data3;"

    out_data = []
    max_freq = 0

    with conn.cursor() as cursor:
        cursor.execute(query)
        histogram_data = cursor.fetchall()
        # print(cursor.description)
        for bucket1, bucket2, frequency, start1, end1, start2, end2, ids in histogram_data:
            if frequency > max_freq:
                max_freq = frequency
            out_data.append(
                {
                    "bucket1": bucket1,
                    "bucket2": bucket2,
                    "frequency": frequency,
                    "start1": round(start1, 0),
                    "end1": round(end1, 0),
                    "start2": round(start2, 3),
                    "end2": round(end2, 3),
                    "ids": ids,
                }
            )

            # print(f"Bucket {bucket1} ({start1:.2f} to {end1:.2f}): {frequency}")
            # print(f"Bucket {bucket2} ({start2:.2f} to {end2:.2f}): {frequency}")
            # print()

        # for bucket, frequency, start, end in histogram_data:
        #     out_data.append(
        #         {"bucket": bucket, "frequency": frequency, "start": start, "end": end}
        #     )

        # for bucket, frequency, start, end in histogram_data:
        #     print(f"Bucket {bucket} ({start:.2f} to {end:.2f}): {frequency}")

    context = {"hist": out_data, "max_freq": max_freq}

    return flask.jsonify(**context)


@webkit.app.route("/api/update_hist", methods=["POST"])
def update_hist():

    conn = connect_to_rds()

    ids = flask.request.get_json()["ids"]

    cursor = conn.cursor()
    out_data = []

    num_buckets = 100
    query = f"""
    SELECT
        width_bucket(min_dist_to_db, rp.min_val1, rp.max_val1, {num_buckets}) AS bucket1,
        width_bucket(network_size, rp.min_val2, rp.max_val2, {num_buckets}) AS bucket2,
        COUNT(*) AS frequency,
        rp.min_val1 + ((width_bucket(min_dist_to_db, rp.min_val1, rp.max_val1, {num_buckets}) - 1) * rp.bucket_width1) AS range_start1,
        rp.min_val1 + (width_bucket(min_dist_to_db, rp.min_val1, rp.max_val1, {num_buckets}) * rp.bucket_width1) AS range_end1,
        rp.min_val2 + ((width_bucket(network_size, rp.min_val2, rp.max_val2, {num_buckets}) - 1) * rp.bucket_width2) AS range_start2,
        rp.min_val2 + (width_bucket(network_size, rp.min_val2, rp.max_val2, {num_buckets}) * rp.bucket_width2) AS range_end2
    FROM
        test_table3,
        range_params_view rp
    WHERE
        id = ANY(%s)
    GROUP BY
        bucket1, bucket2, rp.min_val1, rp.max_val1, rp.bucket_width1, rp.min_val2, rp.max_val2, rp.bucket_width2
    ORDER BY
        bucket1, bucket2;
    """
    cursor.execute(query, (ids,))
    histogram_data = cursor.fetchall()
    # for k in histogram_data:
        # print(k)
        # break
    max_freq = 0
    for bucket1, bucket2, frequency, start1, end1, start2, end2 in histogram_data:
        if frequency > max_freq:
            max_freq = frequency
        out_data.append(
            {
                "bucket1": bucket1,
                "bucket2": bucket2,
                "frequency": frequency,
                "start1": round(start1, 2),
                "end1": round(end1, 2),
                "start2": round(start2, 0),
                "end2": round(end2, 0),
            }
        )

    cursor.close()
    conn.close()

    context = {"hist": out_data, "max_freq": max_freq}

    return flask.jsonify(**context)

# hits = pd.read_excel(webkit.app.config["DATA_FOLDER2"] + "/input_smiles_index.xlsx")


# @webkit.app.route("/api/get_hits", methods=["GET"])
# def get_hits():
#     context = {}

#     context["hits"] = hits.to_dict(orient="records")

#     ums = []
#     for k in context["hits"]:
#         ums.extend(k["unmapped_smiles"].split("."))

#     # print(len(set(ums)))

#     out_inputs = []

#     for i in set(ums):
#         out_inputs.append({"smiles": i, "img": getImage(i)})
#     # print(len(out_inputs), "hello")
#     context["inputs"] = out_inputs

#     return flask.jsonify(**context)


def getDrawing(data):
    mol = Chem.MolFromSmarts(data["smiles"])

    Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE)
    # optimize structure in 2D
    # AllChem.EmbedMolecule(mol)
    # AllChem.MMFFOptimizeMolecule(mol)
    AllChem.Compute2DCoords(mol)

    nodes = []
    edges = []
    conf = mol.GetConformer()
    min_x = 0
    for atom in mol.GetAtoms():
        x_cord = conf.GetAtomPosition(atom.GetIdx()).x
        y_cord = conf.GetAtomPosition(atom.GetIdx()).y
        if x_cord < min_x:
            min_x = x_cord
        nodes.append(
            {
                "id": atom.GetIdx(),
                "atomSmarts": atom.GetSmarts(),
                "x": x_cord,
                "y": y_cord,
            }
        )

    for n in nodes:
        n["x"] = n["x"] + np.abs(min_x)
        # print(n["x"])

    for bond in mol.GetBonds():
        if bond.GetSmarts() == "":
            bond2 = "-"
        else:
            bond2 = bond.GetSmarts()
        edges.append(
            {
                "id": bond.GetIdx(),
                "source": bond.GetBeginAtomIdx(),
                "target": bond.GetEndAtomIdx(),
                "bondSmarts": bond2,
            }
        )

    return {"nodes": nodes, "edges": edges}


@webkit.app.route("/api/run_query", methods=["POST"])
def run_query():
    json_data = flask.request.get_json()
    context = {}

    conn = connect_to_rds()
    cur = conn.cursor()
    # print(json_data)

    query = json_data["query"]
    d = json_data["selectedDisplay"]

    # buffer = [x[0] for x in selected]

    # cur.execute(query, (buffer,))
    cur.execute(query)

    records = cur.fetchall()

    # print length
    # print(len(records))

    hit_ids = []
    for r in records:
        hit_ids.append(r[0])



    conn = connect_to_rds()

    if d == "min_delta_lambda_max":
        query = "SELECT * FROM histogram_data3;"
    elif d == "network_size":
        query = "SELECT * FROM histogram_data2;"

    # query = "SELECT * FROM histogram_data3;"

    out_data = []
    max_freq = 0
    removed = 0
    with conn.cursor() as cursor:
        cursor.execute(query)
        histogram_data = cursor.fetchall()
        # print(cursor.description)
        for bucket1, bucket2, frequency, start1, end1, start2, end2, ids in histogram_data:

            # filter out ids that are not in hit_ids
            orig = len(ids)
            ids = [x for x in ids if x in hit_ids]
            removed = removed + (orig - len(ids))
            # update frequency
            frequency = len(ids)

            if frequency > max_freq:
                max_freq = frequency
            out_data.append(
                {
                    "bucket1": bucket1,
                    "bucket2": bucket2,
                    "frequency": frequency,
                    "start1": round(start1, 0),
                    "end1": round(end1, 0),
                    "start2": round(start2, 3),
                    "end2": round(end2, 3),
                    "ids": ids,
                }
            )


    context["hist"] = out_data
    context["max_freq"] = max_freq
    context["removed"] = removed

    cur.close()
    conn.close()

    return flask.jsonify(**context)


@webkit.app.route("/api/get_drawer_nodes_and_edges", methods=["POST"])
def get_mechanism_drawings():
    json_data = flask.request.get_json()
    context = {}


    # for k in json_data:
        # print(k)

    context["intermediates"] = []
    for r in json_data["intermediates"]:
        out = getImage(r, True)
        # print(r)
        context["intermediates"].append(out)

    context["intermediates"].append(getImage(json_data["unmapped_product_smiles"], True))


    # print(json_data["unmapped_product_smiles"])
    # print(context["intermediates"][-1])

    context["nodes"] = []
    context["edges"] = []
    context["arrows"] = []
    for r in json_data["mechanism"]:
        reacts = r.split(">>")[0]
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
            # print(react)
            out = getDrawing({"smiles": react})

            mol = Chem.MolFromSmarts(react)
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE)
            for i, rr in enumerate(out["nodes"]):
                rr["x"] = rr["x"] * 30 + max_x
                rr["y"] = rr["y"] * 30 + 125

                # print(rr["id"], new_id)
                atom = mol.GetAtomWithIdx(i)

                rr["label"] = atom.GetAtomMapNum()

                rr["atom"] = atom.GetSymbol()
                pattern = r"H-?\d+"
                # pattern = r"H\d+(?:,H\d+)*"
                matches = re.findall(pattern, atom.GetSmarts())
                # print(matches, atom.GetSmarts())
                rr["smarts"] = atom.GetSmarts()
                if len(matches) == 0:
                    rr["hydrogens"] = 0
                else:
                    rr["hydrogens"] = int(matches[0][1:])

                # print(atom.GetSmarts())
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

        prods = r.split(">>")[1]
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

        # print(len(all_nodes["edges"]), len(new_edges))
        # print(new_id, len(all_nodes["nodes"]))
        context["nodes"].append(all_nodes["nodes"])
        context["edges"].append(new_edges)

    return flask.jsonify(**context)
