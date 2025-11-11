# -------------------------------------------------------------
#  main.py  —  Backend for the Drug Discovery Dashboard
#  ------------------------------------------------------------
#  Responsibilities:
#   • Serve the static frontend (HTML + JS) under /app
#   • Provide API endpoints for chemical operations
#     (e.g. validating and canonicalizing SMILES strings)
# -------------------------------------------------------------

from fastapi import FastAPI, Query                   # FastAPI core framework + query param helper
from fastapi.responses import FileResponse, JSONResponse  # Used to send JSON or files back to client
from fastapi.staticfiles import StaticFiles          # Lets FastAPI serve static frontend files
from rdkit import Chem                               # RDKit chemistry toolkit for molecule parsing
import os                                            # To build OS-independent file paths
from fastapi.responses import PlainTextResponse
from rdkit.Chem import AllChem

# -------------------------------------------------------------
#  Create the FastAPI application
# -------------------------------------------------------------
app = FastAPI(title="Drug Discovery Dashboard API")

# -------------------------------------------------------------
#  Serve the frontend (static files)
# -------------------------------------------------------------
#  We mount the /frontend directory under the URL path /app
#  so that visiting http://127.0.0.1:8001/app loads index.html.
#  html=True tells FastAPI to serve index.html automatically
#  when someone accesses the folder root.
# -------------------------------------------------------------
frontend_dir = os.path.join(os.path.dirname(__file__), "..", "frontend")
app.mount("/app", StaticFiles(directory=frontend_dir, html=True), name="frontend")


# -------------------------------------------------------------
#  Root route (simple status endpoint)
# -------------------------------------------------------------
#  When you visit http://127.0.0.1:8001/ directly,
#  this just confirms the backend is running.
# -------------------------------------------------------------
@app.get("/")
def root():
    """Simple root message"""
    return {"message": "Backend is running. Visit /app for the dashboard."}


# -------------------------------------------------------------
#  /validate-smiles endpoint
# -------------------------------------------------------------
#  Purpose:
#    • Receive a SMILES string as query parameter.
#    • Use RDKit to check if it’s valid and canonicalize it.
#    • Return {valid: True/False, canonical: "..."} to the frontend.
# -------------------------------------------------------------
@app.get("/validate-smiles")
def validate_smiles(smiles: str = Query(...)):
    # Log the raw input in the backend console for debugging
    print("🧪 Raw SMILES received:", smiles)

    # Try to create a molecule object from the SMILES
    mol = Chem.MolFromSmiles(smiles)

    # If RDKit fails to parse, return an error JSON
    if mol is None:
        return JSONResponse({"valid": False, "error": "Invalid SMILES"})

    # Generate canonical (standardized) SMILES representation
    canonical = Chem.MolToSmiles(mol)

    # Log for backend visibility
    print("✅ Canonical SMILES:", canonical)

    # Return success response
    return {"valid": True, "canonical": canonical}


@app.get("/generate-3d", response_class=PlainTextResponse)
def generate_3d(smiles: str = Query(...)):
    """
    Generate a 3D conformer for the given SMILES and return it as PDB text.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "ERROR: Invalid SMILES"

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)
    pdb_block = Chem.MolToPDBBlock(mol)
    return pdb_block