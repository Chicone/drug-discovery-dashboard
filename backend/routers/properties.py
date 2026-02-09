# backend/routers/properties.py

from fastapi import APIRouter, Query
from backend.services.properties_service import compute_properties

router = APIRouter()


@router.get("")
def get_properties(smiles: str = Query(..., description="SMILES string")):
    return compute_properties(smiles)
