from functools import lru_cache
import requests

PDBe_UNIPROT_MAPPING_URL = (
    "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"
)

@lru_cache(maxsize=32)
def fetch_pdbe_uniprot_mapping(pdb_id: str) -> dict:
    url = PDBe_UNIPROT_MAPPING_URL.format(pdb_id=pdb_id.lower())
    r = requests.get(url, timeout=20.0)
    r.raise_for_status()
    return r.json()


def build_pdb_author_to_uniprot(pdbe_json: dict, pdb_id: str, accession: str) -> dict[int, int]:
    pdb_id = pdb_id.lower()
    root = pdbe_json.get(pdb_id, {})
    uniprot = root.get("UniProt", {})
    acc_block = uniprot.get(accession, {})
    mappings = acc_block.get("mappings", [])

    out = {}

    for m in mappings:
        unp_start = int(m["unp_start"])
        unp_end = int(m["unp_end"])

        start_author = m["start"].get("author_residue_number")
        end_author = m["end"].get("author_residue_number")

        start_internal = int(m["start"]["residue_number"])
        end_internal = int(m["end"]["residue_number"])

        if start_author is not None:
            start_author = int(start_author)
        else:
            start_author = start_internal

        if end_author is not None:
            end_author = int(end_author)
        else:
            offset = start_author - start_internal
            end_author = end_internal + offset

        for i in range(unp_end - unp_start + 1):
            out[start_author + i] = unp_start + i

    return out
