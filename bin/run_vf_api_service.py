#!/usr/bin/env python3
"""Small API service for running VirtualFlask networks.

Request body (POST /api/v1/virtual-flask/run):
{
  "inputs": "SMI1.SMI2" or ["SMI1", "SMI2"],
  "iterations": 10,
  "threshold": 5000,
  "reagents": "",                    # optional, comma-separated if used
  "ring_filter": false,               # optional
  "intramolecular": true,             # optional
  "precalc_prods": ["SMI"]           # optional
}
"""

from __future__ import annotations

import base64
import io
import os
import sys
from pathlib import Path
from typing import Any

from flask import Flask, jsonify, request
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Draw

# Ensure imports work when script is run directly from any location.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from shared.reaction_class import VirtualFlask, returnReactionTemplates  # noqa: E402


RDLogger.DisableLog("rdApp.*")

app = Flask(__name__)
MECHS = returnReactionTemplates()


def _normalize_inputs(inputs: Any) -> list[str]:
    if isinstance(inputs, str):
        return [s for s in inputs.split(".") if s]
    if isinstance(inputs, list):
        return [str(s).strip() for s in inputs if str(s).strip()]
    return []


def _to_molecule(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return mol

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None
    Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE)
    return mol


def _molecule_png_base64(mol) -> str:
    image = Draw.MolToImage(mol, size=(320, 320))
    buf = io.BytesIO()
    image.save(buf, format="PNG")
    return base64.b64encode(buf.getvalue()).decode("utf-8")


@app.get("/api/health")
def health():
    return jsonify({"ok": True})


@app.post("/api/v1/virtual-flask/run")
def run_virtual_flask():
    payload = request.get_json(silent=True) or {}

    inputs = _normalize_inputs(payload.get("inputs"))
    if not inputs:
        return jsonify({"error": "'inputs' is required. Use dot-separated SMILES or a list."}), 400

    iterations = int(payload.get("iterations", 10))
    threshold = payload.get("threshold", 5000)
    threshold = None if threshold is None else int(threshold)

    reagents = str(payload.get("reagents", ""))
    ring_filter = bool(payload.get("ring_filter", False))
    intramolecular = bool(payload.get("intramolecular", True))
    precalc_prods = payload.get("precalc_prods", []) or []

    network = VirtualFlask(MECHS)
    network.charge(inputs, reagents)
    network.run_until_done(
        iters=iterations,
        thresh=threshold,
        ring_filter=ring_filter,
        intramolecular=intramolecular,
        precalc_prods=precalc_prods,
    )

    compounds = []
    seen = set()

    for node_smiles in network.nodes:
        for frag in node_smiles.split("."):
            mol = _to_molecule(frag)
            if mol is None:
                continue

            canonical = Chem.MolToSmiles(mol)
            if canonical in seen:
                continue
            seen.add(canonical)

            compounds.append(
                {
                    "smiles": canonical,
                    "exact_mass": float(Descriptors.ExactMolWt(mol)),
                    "image_png_base64": _molecule_png_base64(mol),
                }
            )

    response = {
        "inputs": inputs,
        "iterations": iterations,
        "threshold": threshold,
        "end_state": network.end_state,
        "node_count": len(network.nodes),
        "compound_count": len(compounds),
        "results": compounds,
    }
    return jsonify(response)


def main():
    host = os.getenv("VF_API_HOST", "127.0.0.1")
    port = int(os.getenv("VF_API_PORT", "8010"))
    debug = os.getenv("VF_API_DEBUG", "1") == "1"
    app.run(host=host, port=port, debug=debug)


if __name__ == "__main__":
    main()
