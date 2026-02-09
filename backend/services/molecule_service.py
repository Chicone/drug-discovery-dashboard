from fastapi.responses import Response
from rdkit import Chem
from rdkit.Chem import AllChem
from fastapi import UploadFile, File
import os
from pathlib import Path
import json
import requests


def generate_3d_from_smiles(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return {"error": "Invalid SMILES"}
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
    mol_block = Chem.MolToMolBlock(mol)
    return Response(mol_block, media_type="text/plain")

async def save_uploaded_pdb(file: UploadFile):
    path = f"uploaded_structures/{file.filename}"
    os.makedirs("uploaded_structures", exist_ok=True)
    with open(path, "wb") as f:
        f.write(await file.read())
    return {"filename": file.filename}

def fetch_pdb_from_rcsb(pdb_id: str):
    url = f"https://files.rcsb.org/view/{pdb_id}.pdb"
    r = requests.get(url)
    if r.status_code != 200:
        return {"error": "PDB not found"}
    return {"pdb": r.text}

from fastapi import UploadFile, File
def load_sidecar_box(ligand_path: Path):
    sidecar = ligand_path.with_suffix(".box.json")
    if not sidecar.exists():
        return None

    data = json.loads(sidecar.read_text())

    center = data.get("box_center")
    size = data.get("box_size")

    if not (isinstance(center, list) and len(center) == 3):
        return None
    if not (isinstance(size, list) and len(size) == 3):
        return None

    return {
        "center": [float(x) for x in center],
        "size": [float(x) for x in size],
    }

async def load_ligand(ligand: UploadFile = File(...)):
    smiles_dir = Path("smiles")
    smiles_dir.mkdir(exist_ok=True)

    ligand_path = smiles_dir / ligand.filename

    with open(ligand_path, "wb") as f:
        f.write(await ligand.read())

    box = load_sidecar_box(ligand_path)

    resp = {
        "ligand": ligand.filename,
        "has_box": box is not None,
    }

    if box:
        resp["box"] = box
        resp["box_source"] = "sidecar_json"

    return resp


