from rdkit import Chem, Geometry
from PIL import Image
import base64
import cv2
from rdkit.Chem.Draw import rdMolDraw2D
import numpy as np
from rdkit.Chem import AllChem, rdchem
import io


def remove_atom_mapping(in_smiles):
    mol = Chem.MolFromSmiles(in_smiles, sanitize=False)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
        atom.SetIsotope(0)
    return Chem.MolToSmiles(mol)


def find_convex_hull(base64_image_str):
    # Load the image with alpha channel
    # image = cv2.imread(image_path, cv2.IMREAD_UNCHANGED)

    # print(base64_image_str)
    image_data = base64.b64decode(base64_image_str)

    image_array = np.frombuffer(image_data, dtype=np.uint8)

    # Decode the image from the array.
    image = cv2.imdecode(image_array, cv2.IMREAD_UNCHANGED)

    # Check if the image has an alpha channel
    if image.shape[2] == 4:
        # Separate the alpha channel
        _, _, _, alpha = cv2.split(image)

        # Use the alpha channel as the mask to identify non-transparent areas
        mask = cv2.threshold(alpha, 1, 255, cv2.THRESH_BINARY)[1]

        # Apply the mask to the image (assuming the background is white)
        # You might need to adjust this if your background is different
        bgr_image = cv2.cvtColor(image, cv2.COLOR_BGRA2BGR)
        image_with_background = cv2.bitwise_and(bgr_image, bgr_image, mask=mask)
    else:
        # If there's no alpha channel, convert to BGR (3 channels)
        image_with_background = cv2.cvtColor(image, cv2.COLOR_BGRA2BGR)

    # Convert image to grayscale
    gray = cv2.cvtColor(image_with_background, cv2.COLOR_BGR2GRAY)

    # Threshold the image to create a binary image
    _, thresh = cv2.threshold(gray, 1, 255, cv2.THRESH_BINARY)

    # Find contours in the binary image
    # contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    # simplified_hull_points = []
    contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    all_contours_points = np.vstack([contours[i] for i in range(len(contours))])

    # Calculate the convex hull of the merged contour points
    hull = cv2.convexHull(all_contours_points)

    # Optional: Simplify the hull
    epsilon = 0.01 * cv2.arcLength(hull, True)
    simplified_hull = cv2.approxPolyDP(hull, epsilon, True)

    # Convert hull points to a simple list format
    hull_list = simplified_hull[:, 0, :].tolist()
    # if contours:
    #     # Find the largest contour (optional)
    #     largest_contour = max(contours, key=cv2.contourArea)
    #     hull = cv2.convexHull(largest_contour)
    #     # print(hull)
    #     simplified_hull_points.append(cv2.approxPolyDP(
    #         hull, 1e-10 * cv2.arcLength(hull, True), True
    #     ))

    #     # second biggest contour
    #     if len(contours) > 1:
    #         second_largest_contour = sorted(contours, key=cv2.contourArea)[-2]
    #         second_hull = cv2.convexHull(second_largest_contour)
    #         simplified_hull_points.append(cv2.approxPolyDP(second_hull, 1e-10 * cv2.arcLength(second_hull, True), True)
    #         )

    # cv2.drawContours(image, [hull], -1, (0, 255, 0), 2)

    # Find the convex hull object for each contour
    # hull_points = [cv2.convexHull(cnt) for cnt in contours]

    # Simplify hull points (optional, to reduce the number of vertices)
    # simplified_hull_points = [
    #     cv2.approxPolyDP(hull, 0.00001 * cv2.arcLength(hull, True), True)
    #     for hull in hull_points
    # ]

    # print(len(simplified_hull_points), simplified_hull_points[0])

    # Convert hull points to a simple list format if necessary
    # hull_list = [point[0].tolist() for hull in simplified_hull_points for point in hull]
    # hull_list = [point[0].tolist() for point in simplified_hull_points]

    return hull_list


def save_img(mol, img_path="", highlight_atoms=[], atom_cols={}):
    if type(mol) == str:
        smiles = mol
        mol = Chem.MolFromSmiles(smiles)

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
        output_buffer, format="PNG", dpi=(600, 600)
    )  # You can change 'PNG' to 'JPEG'
    byte_data = output_buffer.getvalue()
    base64_str = base64.b64encode(byte_data).decode("utf-8")
    return base64_str


def generateSMARTSfromDrawer(drawing_data):


    print(drawing_data)

    mol = Chem.RWMol()

    atom_id_to_idx = {}
    n = 0
    for atom_data in drawing_data["drawerNodes"]:
        try:
            atom = Chem.AtomFromSmarts(atom_data["atomSmarts"])
            mol.AddAtom(atom)
        except:
            return "invalid atom"
        atom_id_to_idx[atom_data["id"]] = n
        n += 1

    n2 = 0
    for bond_data in drawing_data["drawerEdges"]:
        begin_atom = atom_id_to_idx[bond_data["source"]]
        end_atom = atom_id_to_idx[bond_data["target"]]
        bond_type = bond_data["bondSmarts"]

        if begin_atom == end_atom:
            continue

        try:
            bond = Chem.BondFromSmarts(bond_type)
            mol.AddBond(begin_atom, end_atom, rdchem.BondType.SINGLE)
            mol.ReplaceBond(n2, bond)
        except:
            return "invalid bond"
        n2 = n2 + 1

    out = Chem.MolToSmarts(mol)
    return out


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
                "bondSmarts": bond2
            }
        )

    return {"nodes": nodes, "edges": edges}
