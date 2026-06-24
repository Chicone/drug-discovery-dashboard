from fastapi import APIRouter, Query
from fastapi import APIRouter, Query, UploadFile, File

from backend.services.molecule_service import (
    generate_3d_from_smiles,
    save_uploaded_pdb,
    fetch_pdb_from_rcsb,
    load_ligand,
)

router = APIRouter()

@router.get("/mol3d")
def mol3d(smiles: str = Query(...)):
    return generate_3d_from_smiles(smiles)

@router.post("/upload_pdb")
async def upload_pdb_route(file: UploadFile = File(...)):
    return await save_uploaded_pdb(file)

@router.get("/fetch_pdb")
def fetch_pdb_route(pdb_id: str):
    return fetch_pdb_from_rcsb(pdb_id)

@router.post("/load_ligand")
async def load_ligand_route(file: UploadFile = File(...)):
    return await load_ligand(file)
