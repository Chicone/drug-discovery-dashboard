from fastapi import FastAPI, Query, UploadFile, File, Form
from fastapi.responses import Response, JSONResponse
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
from rdkit.Chem import QED
from typing import Optional
from pathlib import Path
import json
from typing import Any
import numpy as np
from fastapi.responses import JSONResponse

def _iter_pdbqt_atoms(pdbqt_text: str) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
    drop_prefixes = ("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")

    for line in pdbqt_text.splitlines():
        s = line.strip()
        if not s or s.startswith(drop_prefixes):
            continue
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue

        atom_name = line[12:16].strip()
        res_name = line[17:20].strip()
        chain_id = (line[21:22].strip() or "?")
        res_seq = line[22:26].strip()

        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except Exception:
            continue

        # Element: PDBQT usually has it in atom name, but be robust.
        # Example: " C1 ", " OA ", " NA ", etc.
        elem = atom_name[:2].strip().capitalize()
        if len(elem) == 2 and elem[1].isdigit():
            elem = elem[0]
        if elem not in {"C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "H"}:
            elem = atom_name[:1].strip().capitalize()

        atoms.append({
            "record": "HETATM" if line.startswith("HETATM") else "ATOM",
            "atom_name": atom_name,
            "res_name": res_name,
            "chain_id": chain_id,
            "res_seq": res_seq,
            "element": elem,
            "xyz": (x, y, z),
        })

    return atoms


def _split_pdbqt_models(poses_pdbqt_text: str) -> list[str]:
    # Keep full "MODEL ... ENDMDL" blocks.
    chunks = poses_pdbqt_text.split("MODEL")
    models: list[str] = []
    for c in chunks[1:]:
        models.append("MODEL" + c)
    return models


def analyze_pose_contacts_from_files(
    receptor_pdbqt_text: str,
    ligand_pose_pdbqt_text: str,
    cutoff_res: float = 4.0,
    cutoff_hbond: float = 3.5,
) -> dict[str, Any]:
    rec_atoms = _iter_pdbqt_atoms(receptor_pdbqt_text)
    lig_atoms = _iter_pdbqt_atoms(ligand_pose_pdbqt_text)

    # Heuristic: receptor is ATOM, ligand is HETATM. If ligand got marked
    # as ATOM by conversion, fall back to "everything in pose block".
    rec = [a for a in rec_atoms if a["record"] == "ATOM"]
    lig = [a for a in lig_atoms if a["record"] == "HETATM"] or lig_atoms

    if not rec or not lig:
        return {
            "residues": [],
            "counts": {"hbond_candidates": 0, "hydrophobic": 0, "polar": 0},
        }

    rec_xyz = np.array([a["xyz"] for a in rec], dtype=float)
    lig_xyz = np.array([a["xyz"] for a in lig], dtype=float)

    # Pairwise squared distances (L x R)
    d2 = np.sum((lig_xyz[:, None, :] - rec_xyz[None, :, :]) ** 2, axis=2)

    # (1) residues within cutoff_res
    cutoff2 = cutoff_res * cutoff_res
    close = d2 <= cutoff2
    rec_idx = np.where(close)[1]

    residues_set: set[tuple[str, str, str]] = set()
    for j in rec_idx.tolist():
        a = rec[j]
        residues_set.add((a["chain_id"], a["res_seq"], a["res_name"]))

    residues = [
        {"chain": c, "res_seq": rs, "res_name": rn}
        for (c, rs, rn) in sorted(residues_set)
    ]

    # (2) rough contact counts
    rec_el = np.array([a["element"] for a in rec], dtype=object)
    lig_el = np.array([a["element"] for a in lig], dtype=object)

    is_rec_polar = np.isin(rec_el, ["N", "O", "S"])
    is_lig_polar = np.isin(lig_el, ["N", "O", "S"])
    is_rec_c = (rec_el == "C")
    is_lig_c = (lig_el == "C")

    # H-bond candidates: polar-polar within 3.5 Å (no angle check)
    hb2 = cutoff_hbond * cutoff_hbond
    hb_pairs = (d2 <= hb2) & (is_lig_polar[:, None] & is_rec_polar[None, :])
    hbond_candidates = int(np.count_nonzero(hb_pairs))

    # Hydrophobic contacts: C-C within cutoff_res
    hyd_pairs = (d2 <= cutoff2) & (is_lig_c[:, None] & is_rec_c[None, :])
    hydrophobic = int(np.count_nonzero(hyd_pairs))

    # Polar contacts: any pair with at least one polar atom within cutoff_res
    pol_pairs = (d2 <= cutoff2) & (is_lig_polar[:, None] | is_rec_polar[None, :])
    polar = int(np.count_nonzero(pol_pairs))

    return {
        "residues": residues,
        "counts": {
            "hbond_candidates": hbond_candidates,
            "hydrophobic": hydrophobic,
            "polar": polar,
        },
    }


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
import uuid
from datetime import datetime
from zoneinfo import ZoneInfo

import numpy as np
from fastapi import UploadFile
from fastapi.responses import FileResponse, JSONResponse

from rdkit import Chem
from rdkit.Chem import AllChem

from meeko import MoleculePreparation, PDBQTMolecule, RDKitMolCreate
from vina import Vina

RUNS_DIR = Path(__file__).resolve().parent / "data" / "docking_runs"

@app.post("/api/dock_vina")
async def dock_vina(
    receptor: UploadFile = File(...),
    ligand: UploadFile = File(...),

    # Box (optional: default to your current behavior)
    center_x: Optional[float] = Form(None),
    center_y: Optional[float] = Form(None),
    center_z: Optional[float] = Form(None),
    size_x: float = Form(22.0),
    size_y: float = Form(22.0),
    size_z: float = Form(22.0),

    # Vina params
    exhaustiveness: int = Form(8),
    num_modes: int = Form(9),
    seed: Optional[int] = Form(None),
):

    """
    Dock a ligand into a receptor using Meeko + Vina.
    - Receptor: PDB upload (protein).
    - Ligand: PDB or SMILES (.smi / .txt).
    - Output: merged receptor + docked ligand PDB.
    """
    try:
        tmp = tempfile.mkdtemp(prefix="vina_meeko_")
        print(f"[DEBUG] Working directory: {tmp}")

        def split_pdbqt_models(pdbqt_text: str) -> list[str]:
            """
            Split a multi-model PDBQT string into per-model strings.
            Works for Vina outputs containing MODEL/ENDMDL blocks.
            """
            if "MODEL" not in pdbqt_text:
                return [pdbqt_text]

            models = []
            cur = []
            in_model = False
            for line in pdbqt_text.splitlines(keepends=True):
                if line.startswith("MODEL"):
                    in_model = True
                    cur = [line]
                elif line.startswith("ENDMDL") and in_model:
                    cur.append(line)
                    models.append("".join(cur))
                    in_model = False
                elif in_model:
                    cur.append(line)

            return models if models else [pdbqt_text]

        def parse_vina_score(pdbqt_model: str) -> Optional[float]:
            """
            Extract score from a Vina pose block.
            Looks for: 'REMARK VINA RESULT:   -7.5 ...'
            """
            for line in pdbqt_model.splitlines():
                if "VINA RESULT:" in line:
                    # Example: REMARK VINA RESULT: -7.5  0.0  0.0
                    parts = line.split()
                    for p in parts:
                        try:
                            return float(p)
                        except ValueError:
                            continue
            return None

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

        # Choose center
        def is_zero_center(x, y, z, eps=1e-6) -> bool:
            return (x is not None and y is not None and z is not None
                    and abs(float(x)) < eps and abs(float(y)) < eps and abs(float(z)) < eps)

        if (
                center_x is None
                or center_y is None
                or center_z is None
                or is_zero_center(center_x, center_y, center_z)
        ):
            # auto-center on receptor (your current behavior)
            center = np.mean(coords, axis=0).tolist()
            print(f"[DEBUG] Auto-center (receptor centroid): {center}")
        else:
            center = [float(center_x), float(center_y), float(center_z)]
            print(f"[DEBUG] User center: {center}")

        box_size = [float(size_x), float(size_y), float(size_z)]

        v.compute_vina_maps(center=center, box_size=box_size)

        # ---------------------------------------------------------
        # 6) Run docking
        # ---------------------------------------------------------
        dock_kwargs = {
            "exhaustiveness": int(exhaustiveness),
            "n_poses": int(num_modes),
        }
        if seed is not None:
            dock_kwargs["seed"] = int(seed)

        v.dock(**dock_kwargs)


        # ---------------------------------------------------------
        # 7) Prepare receptor RDKit mol once (we'll reuse it)
        # ---------------------------------------------------------
        fixed_rec_pdb = os.path.join(tmp, "receptor_fixed.pdb")
        fix_pdb_format(rec_pdb, fixed_rec_pdb)

        rec_mol = Chem.MolFromPDBFile(fixed_rec_pdb, removeHs=False)
        if rec_mol is None:
            raise ValueError(
                "Receptor PDB could not be parsed by RDKit even after cleaning"
            )
        try:
            rec_mol.GetConformer()
        except Exception:
            raise ValueError("Receptor has no conformer (3D coordinates missing)")

        # ---------------------------------------------------------
        # 8) Get ALL poses from Vina as multi-model PDBQT text
        # ---------------------------------------------------------
        try:
            vina_pdbqt_all = v.poses(n_poses=int(num_modes))
        except TypeError:
            # some vina versions don't accept n_poses here
            vina_pdbqt_all = v.poses()

        if isinstance(vina_pdbqt_all, (list, tuple)):
            # Some builds return list of pose strings already
            pose_pdbqt_models = list(vina_pdbqt_all)
        else:
            # Most builds return one multi-model PDBQT string
            pose_pdbqt_models = split_pdbqt_models(vina_pdbqt_all)

        # ---------------------------------------------------------
        # 9) Convert each pose to RDKit + merge with receptor
        # ---------------------------------------------------------
        poses = []
        for i, pose_pdbqt in enumerate(pose_pdbqt_models):
            score_i = parse_vina_score(pose_pdbqt)

            pdbqt_mol = PDBQTMolecule(pose_pdbqt, skip_typing=True)
            rdkit_mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)

            if not rdkit_mols or rdkit_mols[0] is None:
                continue

            lig_mol = rdkit_mols[0]
            try:
                lig_mol.GetConformer()
            except Exception:
                continue

            merged = Chem.CombineMols(rec_mol, lig_mol)
            merged_pdb_text = Chem.MolToPDBBlock(merged)

            poses.append(
                {
                    "mode": i + 1,
                    "score": score_i,
                    "pdb": merged_pdb_text,
                }
            )

        if not poses:
            raise ValueError("No valid poses produced by Vina/Meeko conversion")

        best_score = min(
            [p["score"] for p in poses if p["score"] is not None],
            default=None,
        )
        print(f"[DEBUG] Best score from parsed poses: {best_score:.3f}")

        # ---------------------------------------------------------
        # Persist docking run (Run History v1)
        # ---------------------------------------------------------
        RUNS_DIR.mkdir(parents=True, exist_ok=True)

        run_id = str(uuid.uuid4())
        run_dir = RUNS_DIR / run_id
        run_dir.mkdir(exist_ok=False)

        created_at = datetime.now(
            ZoneInfo("Europe/Paris")
        ).isoformat(timespec="seconds")

        # Save poses (raw Vina multi-model output)
        poses_pdbqt_path = run_dir / "poses.pdbqt"

        with open(poses_pdbqt_path, "w", encoding="utf-8") as f:
            f.write("\n".join(pose_pdbqt_models).strip() + "\n")

        # Save receptor used for docking (needed to reload full complex)
        receptor_pdbqt_out = run_dir / "receptor.pdbqt"
        shutil.copyfile(rec_pdbqt, receptor_pdbqt_out)

        scores = [p.get("score") for p in poses]

        if best_score is None:
            best_pose_index = None
        else:
            best_pose_index = min(
                range(len(scores)),
                key=lambda i: scores[i] if scores[i] is not None else float("inf"),
            )

        run_record = {
            "run_id": run_id,
            "created_at": created_at,
            "ligand": Path(ligand.filename).stem,
            "receptor": receptor.filename,
            "box": {
                "center": [float(x) for x in center],
                "size": [float(x) for x in box_size],
            },
            "vina": {
                "exhaustiveness": int(exhaustiveness),
                "num_modes": int(num_modes),
            },
            "results": {
                "scores": scores,
                "best_pose_index": best_pose_index,
            },
        }

        (run_dir / "run.json").write_text(
            json.dumps(run_record, indent=2),
            encoding="utf-8",
        )

        return JSONResponse(
            {
                "run_id": run_id,
                "created_at": created_at,
                "best_score": best_score,
                "center": center,
                "box_size": box_size,
                "exhaustiveness": int(exhaustiveness),
                "num_modes": int(num_modes),
                "seed": int(seed) if seed is not None else None,
                "poses": poses,
            }
        )


    except Exception as e:
        print(f"[ERROR] {e}")
        return JSONResponse(
            status_code=500,
            content={"error": str(e)},
        )

@app.get("/api/docking/runs")
def list_docking_runs():
    runs_dir = RUNS_DIR
    if not runs_dir.exists():
        return []

    runs = []
    for run_dir in runs_dir.iterdir():
        if not run_dir.is_dir():
            continue
        run_json = run_dir / "run.json"
        if not run_json.exists():
            continue

        try:
            data = json.loads(run_json.read_text(encoding="utf-8"))
            scores = (
                data.get("results", {})
                    .get("scores", [])
            )
            best = None
            numeric_scores = [s for s in scores if isinstance(s, (int, float))]
            if numeric_scores:
                best = min(numeric_scores)

            runs.append({
                "run_id": data.get("run_id"),
                "created_at": data.get("created_at"),
                "ligand": data.get("ligand"),
                "receptor": data.get("receptor"),
                "best_score": best,
            })
        except Exception:
            continue

    runs.sort(key=lambda r: r.get("created_at") or "", reverse=True)
    return runs

def pdbqt_to_pdb_text(pdbqt_text: str) -> str:
    keep_prefixes = ("ATOM", "HETATM", "MODEL", "ENDMDL", "TER", "END", "REMARK")
    drop_prefixes = ("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")

    out_lines = []
    for line in pdbqt_text.splitlines():
        s = line.strip()
        if not s:
            continue

        if s.startswith(drop_prefixes):
            # AutoDock torsion-tree directives: not valid PDB, breaks many viewers
            continue

        if s.startswith(("ATOM", "HETATM")):
            # Strip PDBQT extras (charge + AutoDock atom type) by truncating
            out_lines.append(line[:78])
            continue

        if s.startswith(keep_prefixes):
            out_lines.append(line)
            continue

        # Drop any other unknown lines
        continue

    return "\n".join(out_lines) + "\n"



@app.get("/api/docking/runs/{run_id}")
def get_docking_run(run_id: str):
    run_dir = RUNS_DIR / run_id
    run_json = run_dir / "run.json"
    poses_file = run_dir / "poses.pdbqt"

    if not run_json.exists() or not poses_file.exists():
        return JSONResponse(
            status_code=404,
            content={"error": "Run not found"},
        )

    data = json.loads(run_json.read_text(encoding="utf-8"))

    # Rebuild poses in the SAME format as /api/dock_vina returns
    receptor_pdbqt_file = run_dir / "receptor.pdbqt"
    receptor_pdb = ""
    if receptor_pdbqt_file.exists():
        receptor_pdbqt_text = receptor_pdbqt_file.read_text(encoding="utf-8")
        receptor_pdb = pdbqt_to_pdb_text(receptor_pdbqt_text).strip() + "\n"

    poses = []
    raw_models = poses_file.read_text(encoding="utf-8").split("MODEL")
    scores = data["results"]["scores"]

    for i, model in enumerate(raw_models[1:]):
        ligand_pdbqt_pose = "MODEL" + model
        ligand_pdb_pose = pdbqt_to_pdb_text(ligand_pdbqt_pose).strip() + "\n"

        # Return a "complex" PDB-like text: receptor + this ligand pose
        complex_pdb = receptor_pdb + ligand_pdb_pose

        poses.append({
            "mode": i + 1,
            "score": scores[i] if i < len(scores) else None,
            "pdb": complex_pdb,
        })

    return {
        "run": data,
        "poses": poses,
    }

@app.get("/api/docking/runs/{run_id}/analyze/{pose_idx}")
def analyze_run_pose(run_id: str, pose_idx: int, cutoff: float = 4.0):
    run_dir = RUNS_DIR / run_id
    rec_file = run_dir / "receptor.pdbqt"
    poses_file = run_dir / "poses.pdbqt"

    if not rec_file.exists() or not poses_file.exists():
        return JSONResponse(status_code=404, content={"error": "Run not found"})

    receptor_text = rec_file.read_text(encoding="utf-8")
    models = _split_pdbqt_models(poses_file.read_text(encoding="utf-8"))

    if pose_idx < 0 or pose_idx >= len(models):
        return JSONResponse(status_code=404, content={"error": "Pose not found"})

    out = analyze_pose_contacts_from_files(
        receptor_pdbqt_text=receptor_text,
        ligand_pose_pdbqt_text=models[pose_idx],
        cutoff_res=float(cutoff),
    )
    return out



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

@app.post("/api/load_ligand")
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
