from fastapi import APIRouter
from backend.services.pdb_service import (
    fetch_pdbe_uniprot_mapping,
    build_pdb_author_to_uniprot,
)

router = APIRouter()

@router.get("/{pdb_id}/pdb_to_uniprot")
def pdb_to_uniprot(pdb_id: str, accession: str = "P29274"):
    pdbe_json = fetch_pdbe_uniprot_mapping(pdb_id)
    mapping = build_pdb_author_to_uniprot(pdbe_json, pdb_id, accession)

    return {
        "pdb_id": pdb_id.lower(),
        "accession": accession,
        "pdb_to_uniprot": {str(k): v for k, v in mapping.items()},
    }
