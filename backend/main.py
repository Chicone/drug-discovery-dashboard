from fastapi import FastAPI, Query
from fastapi.responses import Response
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
from rdkit.Chem import QED

app = FastAPI()


@app.get("/api/properties")
def get_properties(smiles: str = Query(..., description="SMILES string")):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return {"error": "Invalid SMILES"}

    # --- base property computation ---
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hba = Lipinski.NumHAcceptors(mol)
    hbd = Lipinski.NumHDonors(mol)
    rb = Lipinski.NumRotatableBonds(mol)
    qed_score = round(QED.qed(mol), 3)

    props = {
        "MolecularWeight": mw,
        "LogP": logp,
        "TPSA": tpsa,
        "NumHAcceptors": hba,
        "NumHDonors": hbd,
        "RotatableBonds": rb,
        "QED": qed_score,
    }

    insights = []

    # ================================================================
    # 1️⃣  SIZE AND COMPLEXITY CHECKS
    # ================================================================
    heavy_atoms = mol.GetNumHeavyAtoms()

    if mw < 150 or heavy_atoms < 10:
        insights.append(
            "Molecule is very small (<150 Da or <10 heavy atoms). "
            "Such compounds (e.g., ethanol, acetone) can be biologically active "
            "but usually lack the size and complexity required for selective, drug-like behaviour."
        )

        if logp < 0:
            insights.append("LogP is low/negative, indicating a polar compound.")
        elif logp > 5:
            insights.append("LogP is unusually high for such a small molecule.")
        props["insights"] = insights
        return props

    # ================================================================
    # 2️⃣  HIGH-MW REGION
    # ================================================================
    if mw > 700:
        insights.append(
            "Molecular weight is high (>700 Da), moving it beyond classical small-molecule space; "
            "oral absorption may require special design strategies."
        )

    # ================================================================
    # 3️⃣  LIPINSKI SUMMARY (apply only for 150–700 Da)
    # ================================================================
    violations = 0
    if mw > 500:
        violations += 1
    if logp > 5:
        violations += 1
    if hba > 10:
        violations += 1
    if hbd > 5:
        violations += 1

    if 150 <= mw <= 700:
        if violations == 0:
            insights.append(
                "Physicochemical properties are consistent with Lipinski's oral drug-like region."
            )
        elif violations == 1:
            insights.append(
                "One Lipinski parameter is outside the usual range; may still be acceptable depending on the compound series."
            )
        else:
            insights.append(
                f"{violations} Lipinski parameters are outside the usual range; "
                "consider optimization for oral bioavailability."
            )
    else:
        insights.append(
            "Compound lies outside the size range where Lipinski's rules are typically informative."
        )

    # ================================================================
    # 4️⃣  POLARITY / PERMEABILITY
    # ================================================================
    if tpsa > 140:
        insights.append("TPSA is high (>140 Å²), associated with poor passive permeability.")
    elif 60 <= tpsa <= 140:
        insights.append("TPSA is moderate (60–140 Å²), often acceptable for oral drugs.")
    else:
        insights.append("TPSA is low (<60 Å²), generally favourable for permeability.")

    # ================================================================
    # 5️⃣  LIPOPHILICITY (LogP)
    # ================================================================
    if logp < 0:
        insights.append("LogP is low/negative; compound is polar and may need transporters.")
    elif 0 <= logp <= 5:
        insights.append("LogP is within the broadly acceptable 0–5 range.")
    else:
        insights.append("LogP is high (>5), which may impair solubility and increase off-target risk.")

    # ================================================================
    # 6️⃣  HYDROGEN BONDING
    # ================================================================
    if hba > 10 or hbd > 5:
        insights.append(
            "Hydrogen bond donors/acceptors exceed typical Lipinski limits; "
            "this may reduce membrane permeability."
        )
    else:
        insights.append("Hydrogen bonding parameters are within typical drug-like limits.")

    # ================================================================
    # 7️⃣  FLEXIBILITY
    # ================================================================
    if rb > 10:
        insights.append(
            "Rotatable bond count is high (>10), which can reduce bioavailability and binding selectivity."
        )
    else:
        insights.append("Molecular flexibility is within a commonly acceptable range.")

    props["insights"] = insights
    return props


# ================================================================
# 3D COORDINATE ENDPOINT
# ================================================================
@app.get("/api/mol3d")
def mol3d(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return {"error": "Invalid SMILES"}
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
    mol_block = Chem.MolToMolBlock(mol)
    return Response(mol_block, media_type="text/plain")


from fastapi import UploadFile, File


@app.post("/api/upload_pdb")
async def upload_pdb(file: UploadFile = File(...)):
    path = f"uploaded_structures/{file.filename}"
    os.makedirs("uploaded_structures", exist_ok=True)
    with open(path, "wb") as f:
        f.write(await file.read())
    return {"filename": file.filename}


@app.get("/api/fetch_pdb")
def fetch_pdb(pdb_id: str):
    import requests
    url = f"https://files.rcsb.org/view/{pdb_id}.pdb"
    r = requests.get(url)
    if r.status_code != 200:
        return {"error": "PDB not found"}
    return {"pdb": r.text}


# ================================================================
# DOCKING ENDPOINT (AutoDock Vina)
# ================================================================
import os
import shutil
import tempfile

import numpy as np
from fastapi import UploadFile
from fastapi.responses import FileResponse, JSONResponse

from rdkit import Chem
from rdkit.Chem import AllChem

from meeko import MoleculePreparation, PDBQTMolecule, RDKitMolCreate
from vina import Vina


@app.post("/api/dock_vina")
async def dock_vina(receptor: UploadFile, ligand: UploadFile):
    """
    Dock a ligand into a receptor using Meeko + Vina.
    - Receptor: PDB upload (protein).
    - Ligand: PDB or SMILES (.smi / .txt).
    - Output: merged receptor + docked ligand PDB.
    """
    try:
        tmp = tempfile.mkdtemp(prefix="vina_meeko_")
        print(f"[DEBUG] Working directory: {tmp}")

        # ---------------------------------------------------------
        # Helper: save UploadFile to disk safely
        # ---------------------------------------------------------
        def save_upload(upload: UploadFile, path: str) -> None:
            upload.file.seek(0)
            with open(path, "wb") as f:
                shutil.copyfileobj(upload.file, f)

        # ---------------------------------------------------------
        # Helper: fix receptor PDB so RDKit can parse it
        # (only rewrite ATOM/HETATM lines, keep headers)
        # ---------------------------------------------------------
        def fix_pdb_format(input_pdb: str, output_pdb: str) -> None:
            with open(input_pdb) as fin, open(output_pdb, "w") as fout:
                for line in fin:
                    if line.startswith(("ATOM", "HETATM")):
                        atom = line[:6]
                        serial = line[6:11]
                        name = line[12:16]
                        altLoc = line[16]
                        resName = line[17:20]
                        chainID = line[21]
                        resSeq = line[22:26]
                        iCode = line[26]
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        element = line[76:78].strip()

                        fout.write(
                            f"{atom}{serial} {name}{altLoc}{resName} "
                            f"{chainID}{resSeq}{iCode}   "
                            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00"
                            f"          {element:>2}\n"
                        )
                    else:
                        fout.write(line)

        # ---------------------------------------------------------
        # 1) Save uploaded receptor and ligand
        # ---------------------------------------------------------
        rec_pdb = os.path.join(tmp, "receptor.pdb")
        lig_path = os.path.join(tmp, ligand.filename)

        save_upload(receptor, rec_pdb)
        save_upload(ligand, lig_path)

        print("[DEBUG] receptor size:", os.path.getsize(rec_pdb))
        print("[DEBUG] ligand size:", os.path.getsize(lig_path))

        # ---------------------------------------------------------
        # 2) Convert receptor to PDBQT for Vina (Open Babel)
        # ---------------------------------------------------------
        rec_pdbqt = os.path.join(tmp, "receptor.pdbqt")
        cmd = (
            f"obabel -ipdb {rec_pdb} -xr -xp -d "
            f"-O {rec_pdbqt}"
        )
        ret = os.system(cmd)
        if ret != 0 or not os.path.exists(rec_pdbqt):
            raise RuntimeError("Open Babel failed to make receptor.pdbqt")

        # ---------------------------------------------------------
        # 3) Load ligand as RDKit molecule
        # ---------------------------------------------------------
        if ligand.filename.endswith((".smi", ".txt")):
            with open(lig_path) as f:
                smiles = f.read().strip().split()[0]
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("Could not parse SMILES from ligand file")
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)
        else:
            mol = Chem.MolFromPDBFile(lig_path, removeHs=False)

        if mol is None:
            raise ValueError("Could not load ligand as RDKit molecule")

        # ---------------------------------------------------------
        # 4) Prepare ligand for Vina with Meeko
        # ---------------------------------------------------------
        prep = MoleculePreparation()
        result = prep.prepare(mol)
        if isinstance(result, list):
            result = result[0]

        # This works with your current Meeko version
        prep.molsetup = result
        pdbqt_string = prep.write_pdbqt_string()

        # ---------------------------------------------------------
        # 5) Set up Vina and compute grid center from receptor PDB
        # ---------------------------------------------------------
        v = Vina(sf_name="vina")
        v.set_receptor(rec_pdbqt)
        v.set_ligand_from_string(pdbqt_string)

        coords = []
        with open(rec_pdb) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))

        if not coords:
            raise ValueError("No ATOM/HETATM lines in receptor PDB")

        center = np.mean(coords, axis=0)
        print(f"[DEBUG] Auto-center: {center}")

        v.compute_vina_maps(center=center.tolist(), box_size=[20, 20, 20])

        # ---------------------------------------------------------
        # 6) Run docking
        # ---------------------------------------------------------
        v.dock(exhaustiveness=8, n_poses=5)
        score = v.score()
        if isinstance(score, (list, np.ndarray)):
            best_score = score[0]
        else:
            best_score = score
        print(f"[DEBUG] Best score: {best_score:.3f}")

        # ---------------------------------------------------------
        # 7) Get best pose and convert back to RDKit with Meeko
        # ---------------------------------------------------------
        vina_output = v.poses()
        if isinstance(vina_output, (list, tuple)):
            vina_pdbqt_str = vina_output[0]
        else:
            vina_pdbqt_str = vina_output

        pdbqt_mol = PDBQTMolecule(vina_pdbqt_str, skip_typing=True)
        rdkit_mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)

        if not rdkit_mols or rdkit_mols[0] is None:
            raise ValueError("Failed to create RDKit mol from Vina PDBQT pose")

        mol_docked = rdkit_mols[0]
        docked_pdb = os.path.join(tmp, "docked_meeko.pdb")
        Chem.MolToPDBFile(mol_docked, docked_pdb)

        # ---------------------------------------------------------
        # 8) Fix receptor PDB formatting and load both with RDKit
        # ---------------------------------------------------------
        fixed_rec_pdb = os.path.join(tmp, "receptor_fixed.pdb")
        fix_pdb_format(rec_pdb, fixed_rec_pdb)

        rec_mol = Chem.MolFromPDBFile(fixed_rec_pdb, removeHs=False)
        lig_mol = Chem.MolFromPDBFile(docked_pdb, removeHs=False)

        if rec_mol is None:
            raise ValueError(
                "Receptor PDB could not be parsed by RDKit "
                "even after cleaning"
            )
        if lig_mol is None:
            raise ValueError(
                "Docked ligand PDB could not be parsed by RDKit"
            )

        # ensure they have conformers (3D coords)
        try:
            rec_mol.GetConformer()
            lig_mol.GetConformer()
        except Exception:
            raise ValueError(
                "Receptor or ligand has no conformer (3D coordinates missing)"
            )

        # ---------------------------------------------------------
        # 9) Merge receptor + ligand and write final complex
        # ---------------------------------------------------------
        merged = Chem.CombineMols(rec_mol, lig_mol)
        merged_pdb = os.path.join(tmp, "complex_meeko.pdb")
        Chem.MolToPDBFile(merged, merged_pdb)

        print(f"[DEBUG] RDKit merged complex written to {merged_pdb}")

        return FileResponse(
            merged_pdb,
            media_type="chemical/x-pdb",
            filename="docked_complex.pdb",
        )

    except Exception as e:
        print(f"[ERROR] {e}")
        return JSONResponse(
            status_code=500,
            content={"error": str(e)},
        )
