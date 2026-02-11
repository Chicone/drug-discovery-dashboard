from fastapi import FastAPI, Query, UploadFile, File, Form
from backend.routers.md import router as md_router
from backend.routers.docking import router as docking_router
from backend.routers.pdb import router as pdb_router
from backend.routers.properties import router as properties_router
from backend.routers.molecules import router as molecules_router
from backend.routers.analysis import router as analysis_router


app = FastAPI()

app.include_router(md_router, prefix="/api/md")
app.include_router(docking_router)
app.include_router(pdb_router, prefix="/api/pdb", tags=["PDB Mapping"])
app.include_router(properties_router, prefix="/api/properties", tags=["Properties"])
app.include_router(molecules_router, prefix="/api/molecules", tags=["Molecules"])
app.include_router(analysis_router, prefix="/api/analysis", tags=["Analysis"])

print(">>> RUNTIME: runtime.py LOADED", __file__)

