from rdkit import Chem, Geometry
from PIL import Image
import base64
import cv2
from rdkit.Chem.Draw import rdMolDraw2D
import numpy as np
from rdkit.Chem import AllChem
import io
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib


def get_networkx_graph_from_state_network(hyp):
    # convert into networkx graph
    hypx = nx.DiGraph()
    for h in hyp.nodes:
        hypx.add_node(h)

    for h in hyp.nodes:
        if h in hyp.edges:
            for n in hyp.edges[h]:
                hypx.add_edge(h, n)

    return hypx


def get_path_from_init_to_node(hyp, targ):
    for h in hyp.nodes:
        if hyp[h].propagations == 0:
            init_node = h
            break

    path = nx.shortest_path(hyp.nx, init_node, targ)

    return path


def remove_atom_mapping(in_smiles):
    mol = Chem.MolFromSmiles(in_smiles, sanitize=False)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
        atom.SetIsotope(0)
    return Chem.MolToSmiles(mol)


def remove_atom_mapping_mol(in_mol):
    for atom in in_mol.GetAtoms():
        atom.SetAtomMapNum(0)
        atom.SetIsotope(0)


def split_at_period_not_in_parentheses(s):
    parts = []
    current_part = []
    depth = 0  # Track the depth of parentheses nesting

    for char in s:
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        elif char == "." and depth == 0:
            # If we're not inside parentheses, treat the period as a delimiter
            parts.append("".join(current_part))
            current_part = []
            continue

        # Add the current character to the part we're building
        current_part.append(char)

    # Add the last part, if there is one
    if current_part:
        parts.append("".join(current_part))

    return parts


def find_convex_hull(base64_image_str):

    image_data = base64.b64decode(base64_image_str)
    image_array = np.frombuffer(image_data, dtype=np.uint8)
    image = cv2.imdecode(image_array, cv2.IMREAD_UNCHANGED)

    if image.shape[2] == 4:
        # Separate the alpha channel
        _, _, _, alpha = cv2.split(image)
        mask = cv2.threshold(alpha, 1, 255, cv2.THRESH_BINARY)[1]
        bgr_image = cv2.cvtColor(image, cv2.COLOR_BGRA2BGR)
        image_with_background = cv2.bitwise_and(bgr_image, bgr_image, mask=mask)
    else:
        image_with_background = cv2.cvtColor(image, cv2.COLOR_BGRA2BGR)

    gray = cv2.cvtColor(image_with_background, cv2.COLOR_BGR2GRAY)
    _, thresh = cv2.threshold(gray, 1, 255, cv2.THRESH_BINARY)
    contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    all_contours_points = np.vstack([contours[i] for i in range(len(contours))])
    hull = cv2.convexHull(all_contours_points)
    epsilon = 0.01 * cv2.arcLength(hull, True)
    simplified_hull = cv2.approxPolyDP(hull, epsilon, True)
    hull_list = simplified_hull[:, 0, :].tolist()

    return hull_list


def save_img(mol, img_path="", highlight_atoms=[], atom_cols={}):
    if type(mol) == str:
        smiles = mol
        mol = Chem.MolFromSmiles(smiles, sanitize=False)

        if not mol:
            print("mol not formed for: ", smiles)
            mol = Chem.MolFromSmarts(smiles)

    highlight_bonds = []
    bond_cols = {}
    highlight_atom_radii = {idx: 0.5 for idx in highlight_atoms}

    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)

    if highlight_atoms:
        for bond in mol.GetBonds():
            beg = bond.GetBeginAtomIdx()
            end = bond.GetEndAtomIdx()
            if beg in highlight_atoms and end in highlight_atoms:
                if atom_cols[beg] == atom_cols[end]:
                    highlight_bonds.append(bond.GetIdx())
                    bond_cols[bond.GetIdx()] = atom_cols[end]

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

    if highlight_atoms:
        rdMolDraw2D.PrepareAndDrawMolecule(
            drawer,
            mc,
            highlightAtoms=highlight_atoms,
            highlightAtomColors=atom_cols,
            highlightAtomRadii=highlight_atom_radii,
            highlightBonds=highlight_bonds,
            highlightBondColors=bond_cols,
        )
    else:
        drawer.DrawMolecule(mc)

    drawer.FinishDrawing()
    img = drawer.GetDrawingText()
    png_bytes_io = io.BytesIO(img)
    img = Image.open(png_bytes_io)
    img = img.convert("RGBA")
    datas = img.getdata()

    newData = []
    for item in datas:
        if item[0] == 255 and item[1] == 255 and item[2] == 255:  # Finding white pixels
            newData.append((255, 255, 255, 0))  # Changing white pixels to transparent
        else:
            newData.append(item)
    img.putdata(newData)

    output_buffer = io.BytesIO()
    img.save(
        output_buffer, format="PNG", dpi=(600, 600)
    )  # You can change 'PNG' to 'JPEG'
    byte_data = output_buffer.getvalue()
    base64_str = base64.b64encode(byte_data).decode("utf-8")
    return base64_str


def get_state_img(mol, highlight_atoms=[], atom_cols={}):
    if type(mol) == str:
        smiles = mol
        mol = Chem.MolFromSmiles(smiles, sanitize=True)

        if not mol:
            print("mol not formed for: ", smiles)
            mol = Chem.MolFromSmarts(smiles)

    atom_map_to_index = {}
    for atom in mol.GetAtoms():
        atom_map_to_index[atom.GetAtomMapNum()] = atom.GetIdx()
        atom.SetAtomMapNum(0)

    highlight_atoms_r = []
    highlight_bonds = []
    bond_cols = {}
    atom_cols_r = {}
    highlight_atom_radii = {idx: 0.5 for idx in highlight_atoms}
    if highlight_atoms:

        for atom in highlight_atoms:
            highlight_atoms_r.append(atom_map_to_index[atom])
            atom_cols_r[atom_map_to_index[atom]] = atom_cols[atom]

        for bond in mol.GetBonds():
            beg = bond.GetBeginAtomIdx()
            end = bond.GetEndAtomIdx()
            if beg in highlight_atoms_r and end in highlight_atoms_r:
                if atom_cols_r[beg] == atom_cols_r[end]:
                    highlight_bonds.append(bond.GetIdx())
                    bond_cols[bond.GetIdx()] = atom_cols_r[end]

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

    if highlight_atoms:
        rdMolDraw2D.PrepareAndDrawMolecule(
            drawer,
            mc,
            highlightAtoms=highlight_atoms_r,
            highlightAtomColors=atom_cols_r,
            highlightAtomRadii=highlight_atom_radii,
            highlightBonds=highlight_bonds,
            highlightBondColors=bond_cols,
        )
    else:
        drawer.DrawMolecule(mc)

    drawer.FinishDrawing()
    img = drawer.GetDrawingText()
    png_bytes_io = io.BytesIO(img)
    img = Image.open(png_bytes_io)
    img = img.convert("RGBA")
    datas = img.getdata()

    newData = []
    for item in datas:
        if item[0] == 255 and item[1] == 255 and item[2] == 255:  # Finding white pixels
            newData.append((255, 255, 255, 0))  # Changing white pixels to transparent
        else:
            newData.append(item)
    img.putdata(newData)

    output_buffer = io.BytesIO()
    img.save(
        output_buffer, format="PNG", dpi=(600, 600)
    )  # You can change 'PNG' to 'JPEG'
    byte_data = output_buffer.getvalue()
    base64_str = base64.b64encode(byte_data).decode("utf-8")
    return base64_str


def draw_hypergraph(name, hyp, highlight):
    hypx = get_networkx_graph_from_state_network(hyp)
    pos = nx.nx_agraph.graphviz_layout(hypx, prog="dot")

    # propagation_values = [hyp.nodes[node]["propagations"] for node in hyp.nodes]
    # norm = mcolors.Normalize(vmin=min(propagation_values), vmax=max(propagation_values))
    # cmap = plt.cm.inferno

    node_colors_final = []
    for i, k in enumerate(hyp.nodes):
        # if k in highlight:
        #     node_colors_final.append("red")
        #     continue
        if not hyp[k].tcp or hyp[k].contains_non_participating_transformation:
            node_colors_final.append("grey")
            continue
        if hyp[k].failed_thermo_state1 or hyp[k].structural_failure:
            node_colors_final.append("#ffcb3f")
            continue
        if hyp[k].known_product or hyp[k].in_same_mechanism:
            node_colors_final.append("#3b558d")
            continue
        # print(hyp[k].average_distance_to_db)
        if hyp[k].average_distance_to_db < 0.9 or hyp[k].average_distance_to_db == None:
            node_colors_final.append("#bb5ca3")
            continue
        if hyp[k].askcos_hit:
            node_colors_final.append("#3b558d")
            continue
        node_colors_final.append("green")

    fig, ax = plt.subplots()
    fig.set_size_inches(6, 4)

    node_size = 0.005
    G = hypx

    nodes_x, nodes_y = zip(*[pos[n] for n in G.nodes()])

    ax.scatter(
        nodes_x, nodes_y, s=node_size * 1e3, color=node_colors_final, alpha=1, zorder=2
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
        ax.plot(
            [start_x, end_x],
            [start_y, end_y],
            color="black",
            alpha=0.5,
            zorder=1,
            linewidth=0.3,
        )

    plt.axis("off")

    plt.savefig(
        f"{name}.png",
        dpi=900,
        transparent=True,
        bbox_inches="tight",
        pad_inches=0.01,
    )


def draw_hypergraph_no_filter(name, hyp, highlight=[]):
    hypx = get_networkx_graph_from_state_network(hyp)
    pos = nx.nx_agraph.graphviz_layout(hypx, prog="dot")

    propagation_values = [hyp[node].propagations for node in hyp.nodes]
    norm = mcolors.Normalize(vmin=min(propagation_values) - 1, vmax=10)
    # print(matplotlib.colormaps)
    cmap = plt.cm._colormaps.get_cmap("BrBG")

    node_colors = [cmap(norm(value)) for value in propagation_values]

    node_colors_final = []
    for i, k in enumerate(hyp.nodes):
        if k in highlight:
            node_colors_final.append("black")
            continue
        node_colors_final.append(node_colors[i])

    fig, ax = plt.subplots(dpi=300)
    fig.set_size_inches(2, 2)

    node_size = 0.3
    G = hypx

    nodes_x, nodes_y = zip(*[pos[n] for n in G.nodes()])
    ax.scatter(
        nodes_x, nodes_y, s=node_size * 1, color=node_colors_final, alpha=1, zorder=2
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


def draw_mechanism_sequence(G, target):

    path = get_path_from_init_to_node(G, target)

    num_images = len(path)

    fig, axes = plt.subplots(num_images, 1, figsize=(5, num_images * 5))

    # Ensure axes is iterable even if there's only one image
    if num_images == 1:
        axes = [axes]

    for ax, p in zip(axes, path):
        img = G[p].img
        image_bytes = base64.b64decode(img)
        # Convert bytes to a PIL image
        pil_image = Image.open(io.BytesIO(image_bytes))
        # Display the image on the respective subplot
        ax.imshow(pil_image)
        ax.axis("off")  # Hide the axis

    plt.show()


def show_image(img):
    fig, ax = plt.subplots(figsize=(2, 2), dpi=300)
    image_bytes = base64.b64decode(img)
    # Convert bytes to a PIL image
    image = Image.open(io.BytesIO(image_bytes))

    ax.imshow(image)
    plt.axis("off")
    plt.show()


def getImage(row):
    if row == None:
        return None

    mol = Chem.MolFromSmiles(row)

    if not mol:
        return None

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
