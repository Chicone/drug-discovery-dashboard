from fastapi import FastAPI, Query, UploadFile, File, Form
from fastapi.responses import Response, JSONResponse
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
from rdkit.Chem import QED
from typing import Optional
from pathlib import Path
import json
from typing import Any
from functools import lru_cache
import requests
from functools import lru_cache

import numpy as np
from fastapi.responses import JSONResponse

PDBe_UNIPROT_MAPPING_URL = (
    "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"
)

@lru_cache(maxsize=32)
def _fetch_pdbe_uniprot_mapping(pdb_id: str) -> dict:
    url = PDBe_UNIPROT_MAPPING_URL.format(
        pdb_id=pdb_id.lower(),
    )
    r = requests.get(url, timeout=20.0)
    r.raise_for_status()
    return r.json()




def build_pdb_author_to_uniprot(pdbe_json: dict, pdb_id: str, accession: str) -> dict[int, int]:
    """
    Returns {author_resi: uniprot_resi} for a given PDB and UniProt accession.
    This is the mapping you want if your viewer uses PDB residue numbers (atom.resi).
    """
    pdb_id = pdb_id.lower()
    root = pdbe_json.get(pdb_id, {})
    uniprot = root.get("UniProt", {})
    acc_block = uniprot.get(accession, {})
    mappings = acc_block.get("mappings", [])

    out: dict[int, int] = {}

    for m in mappings:
        unp_start = int(m["unp_start"])
        unp_end = int(m["unp_end"])

        start_author = m["start"].get("author_residue_number")
        end_author = m["end"].get("author_residue_number")

        start_internal = int(m["start"]["residue_number"])
        end_internal = int(m["end"]["residue_number"])

        # Prefer author numbering (what your PDB/3Dmol uses)
        if start_author is not None:
            start_author = int(start_author)
        else:
            # fallback: derive author numbering using offset from internal numbering
            start_author = start_internal

        if end_author is not None:
            end_author = int(end_author)
        else:
            # derive end_author using the same offset inferred at the start
            offset = start_author - start_internal
            end_author = end_internal + offset

        # sanity: lengths must match
        expected_len = unp_end - unp_start
        got_len = end_author - start_author
        if expected_len != got_len:
            # If this triggers, something unusual is happening (insertions etc.)
            # You can choose to raise, log, or still do best-effort mapping.
            pass

        for i in range(unp_end - unp_start + 1):
            out[start_author + i] = unp_start + i

    return out



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

    # Minimum distance per residue
    res_min_dist: dict[tuple[str, str, str], float] = {}

    for li in range(d2.shape[0]):
        for ri in range(d2.shape[1]):
            dist = float(np.sqrt(d2[li, ri]))
            if dist > cutoff_res:
                continue
            r = rec[ri]
            key = (r["chain_id"], r["res_seq"], r["res_name"])
            if key not in res_min_dist or dist < res_min_dist[key]:
                res_min_dist[key] = dist

    # Residues sorted by sequence (for contiguity inspection)
    def _res_seq_int(x: str) -> int:
        try:
            return int(x)
        except Exception:
            return 10 ** 9

    residues_seq = [
        {
            "chain": c,
            "res_seq": rs,
            "res_name": rn,
            "min_dist": round(res_min_dist[(c, rs, rn)], 2),
        }
        for (c, rs, rn) in sorted(
            res_min_dist.keys(),
            key=lambda k: (k[0], _res_seq_int(k[1]), k[1])
        )
    ]

    # (1) residues within cutoff_res
    cutoff2 = cutoff_res * cutoff_res
    close = d2 <= cutoff2
    rec_idx = np.where(close)[1]

    residues_set: set[tuple[str, str, str]] = set()
    for j in rec_idx.tolist():
        a = rec[j]
        residues_set.add((a["chain_id"], a["res_seq"], a["res_name"]))

    residues = [
        {
            "chain": c,
            "res_seq": rs,
            "res_name": rn,
            "min_dist": round(res_min_dist[(c, rs, rn)], 2),
        }
        for (c, rs, rn) in sorted(
            res_min_dist.keys(),
            key=lambda k: res_min_dist[k]
        )
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
        "residues_seq": residues_seq,  # sorted by sequence
        "counts": {
            "hbond_candidates": hbond_candidates,
            "hydrophobic": hydrophobic,
            "polar": polar,
        },
    }

app = FastAPI()



@app.get("/api/pdb/{pdb_id}/pdb_to_uniprot")
def pdb_to_uniprot_map(
    pdb_id: str,
    accession: str = "P29274",
):
    pdbe_json = _fetch_pdbe_uniprot_mapping(pdb_id)
    mapping = build_pdb_author_to_uniprot(pdbe_json, pdb_id, accession)

    return {
        "pdb_id": pdb_id.lower(),
        "accession": accession,
        "pdb_to_uniprot": {str(k): v for k, v in mapping.items()},
    }



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

from fastapi import Query  # add near your other fastapi imports
@app.get("/api/docking/runs")
def list_docking_runs(limit: int = Query(0, ge=0)):
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
            scores = data.get("results", {}).get("scores", [])
            best = None
            numeric_scores = [s for s in scores if isinstance(s, (int, float))]
            if numeric_scores:
                best = min(numeric_scores)

            runs.append(
                {
                    "run_id": data.get("run_id"),
                    "created_at": data.get("created_at"),
                    "ligand": data.get("ligand"),
                    "receptor": data.get("receptor"),
                    "best_score": best,
                }
            )
        except Exception:
            continue

    runs.sort(key=lambda r: r.get("created_at") or "", reverse=True)

    # limit=0 means "no limit" (backward compatible)
    if limit and limit > 0:
        return runs[:limit]
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


# ================================================================
# MOLECULAR DYNAMICS (MVP) ENDPOINTS
# ================================================================
import subprocess
import zipfile

MD_RUNS_DIR = Path(__file__).resolve().parent / "data" / "md_runs"


def _md_job_dir(job_id: str) -> Path:
    return MD_RUNS_DIR / job_id


def _safe_mkdir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _write_json(p: Path, obj: dict) -> None:
    p.write_text(json.dumps(obj, indent=2), encoding="utf-8")


def _read_json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def _md_status_from_dir(job_dir: Path) -> str:
    status_file = job_dir / "status.json"
    if status_file.exists():
        try:
            return _read_json(status_file).get("status", "unknown")
        except Exception:
            return "unknown"

    # Fallback inference
    if (job_dir / "out" / "done.marker").exists():
        return "done"
    if (job_dir / "out" / "error.marker").exists():
        return "error"
    return "running"


def _launch_md_docker(job_dir: Path) -> None:
    """
    Launch MD runner container asynchronously.
    Expects an image named: md-runner:latest

    Container contract:
      - Reads:  /job/input/protein.pdb
      - Reads:  /job/input/params.json
      - Writes: /job/log.txt
      - Writes: /job/status.json (optional but recommended)
      - Writes: /job/out/* (outputs)
    """
    # IMPORTANT: This runs Docker from the host. For Docker-in-Docker,
    # you would need the Docker socket mounted into the backend container.
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{str(job_dir)}:/job",
        "md-runner:latest",
    ]
    subprocess.Popen(
        cmd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


@app.post("/api/md/jobs")
async def create_md_job(
    protein_pdb: UploadFile = File(...),
    preset: str = Form("cg_popc_50ns"),

    # NEW:
    scenario: str = Form("protein_only"),
    orthosteric_ligand: Optional[UploadFile] = File(None),
    allosteric_pose: Optional[UploadFile] = File(None),
):
    try:
        allowed_scenarios = {
            "protein_only",
            "protein_plus_orthosteric",
            "protein_plus_orthosteric_plus_allosteric",
        }
        if scenario not in allowed_scenarios:
            return JSONResponse(
                status_code=400,
                content={"error": f"Invalid scenario: {scenario}"},
            )

        # Validate scenario requirements
        if scenario == "protein_plus_orthosteric" and orthosteric_ligand is None:
            return JSONResponse(
                status_code=400,
                content={"error": "Scenario requires orthosteric_ligand"},
            )

        if scenario == "protein_plus_orthosteric_plus_allosteric":
            if orthosteric_ligand is None or allosteric_pose is None:
                return JSONResponse(
                    status_code=400,
                    content={
                        "error": (
                            "Scenario requires orthosteric_ligand "
                            "and allosteric_pose"
                        )
                    },
                )

        MD_RUNS_DIR.mkdir(parents=True, exist_ok=True)

        job_id = str(uuid.uuid4())
        job_dir = _md_job_dir(job_id)
        job_dir.mkdir(exist_ok=False)

        input_dir = job_dir / "input"
        out_dir = job_dir / "out"
        _safe_mkdir(input_dir)
        _safe_mkdir(out_dir)

        # --- Save protein ---
        pdb_path = input_dir / "protein.pdb"
        protein_pdb.file.seek(0)
        with open(pdb_path, "wb") as f:
            shutil.copyfileobj(protein_pdb.file, f)

        # --- Save optional ligands ---
        orth_path = None
        if orthosteric_ligand is not None:
            orth_name = Path(orthosteric_ligand.filename).name
            orth_path = input_dir / orth_name
            orthosteric_ligand.file.seek(0)
            with open(orth_path, "wb") as f:
                shutil.copyfileobj(orthosteric_ligand.file, f)

        allo_path = None
        if allosteric_pose is not None:
            allo_name = Path(allosteric_pose.filename).name
            allo_path = input_dir / allo_name
            allosteric_pose.file.seek(0)
            with open(allo_path, "wb") as f:
                shutil.copyfileobj(allosteric_pose.file, f)

        created_at = datetime.now(
            ZoneInfo("Europe/Paris")
        ).isoformat(timespec="seconds")

        run_record = {
            "job_id": job_id,
            "created_at": created_at,
            "preset": preset,
            "scenario": scenario,
            "protein_filename": protein_pdb.filename,
            "orthosteric_filename": orthosteric_ligand.filename
            if orthosteric_ligand is not None else None,
            "allosteric_pose_filename": allosteric_pose.filename
            if allosteric_pose is not None else None,
        }
        _write_json(job_dir / "run.json", run_record)

        # Initialize status + log so frontend can read immediately
        _write_json(job_dir / "status.json", {"status": "queued"})
        (job_dir / "log.txt").write_text("", encoding="utf-8")

        # Save params for container (future-proof for your md-runner)
        params = {
            "preset": preset,
            "scenario": scenario,
            "files": {
                "protein_pdb": "protein.pdb",
                "orthosteric_ligand": orth_path.name if orth_path else None,
                "allosteric_pose": allo_path.name if allo_path else None,
            },
        }
        _write_json(input_dir / "params.json", params)

        # Launch container
        _launch_md_docker(job_dir)

        # Mark as running
        _write_json(job_dir / "status.json", {"status": "running"})

        return JSONResponse({"job_id": job_id})

    except Exception as e:
        return JSONResponse(status_code=500, content={"error": str(e)})



@app.get("/api/md/jobs")
def list_md_jobs():
    if not MD_RUNS_DIR.exists():
        return []

    jobs = []
    for job_dir in MD_RUNS_DIR.iterdir():
        if not job_dir.is_dir():
            continue
        run_json = job_dir / "run.json"
        if not run_json.exists():
            continue
        try:
            data = _read_json(run_json)
            data["status"] = _md_status_from_dir(job_dir)
            jobs.append(data)
        except Exception:
            continue

    jobs.sort(key=lambda r: r.get("created_at") or "", reverse=True)
    return jobs


@app.get("/api/md/jobs/{job_id}")
def get_md_job(job_id: str):
    job_dir = _md_job_dir(job_id)
    run_json = job_dir / "run.json"
    if not run_json.exists():
        return JSONResponse(status_code=404, content={"error": "Job not found"})

    data = _read_json(run_json)
    data["status"] = _md_status_from_dir(job_dir)
    return data


@app.get("/api/md/jobs/{job_id}/log")
def get_md_job_log(job_id: str):
    job_dir = _md_job_dir(job_id)
    log_path = job_dir / "log.txt"
    if not log_path.exists():
        return JSONResponse(status_code=404, content={"error": "Log not found"})
    return Response(log_path.read_text(encoding="utf-8"), media_type="text/plain")


@app.get("/api/md/jobs/{job_id}/files")
def list_md_job_files(job_id: str):
    job_dir = _md_job_dir(job_id)
    out_dir = job_dir / "out"
    if not out_dir.exists():
        return JSONResponse(status_code=404, content={"error": "Job not found"})

    files = []
    for p in out_dir.rglob("*"):
        if p.is_file():
            rel = p.relative_to(out_dir).as_posix()
            files.append({
                "name": rel,
                "size": p.stat().st_size,
            })

    files.sort(key=lambda x: x["name"])
    return {"files": files}


@app.get("/api/md/jobs/{job_id}/download")
def download_md_job(job_id: str):
    job_dir = _md_job_dir(job_id)
    out_dir = job_dir / "out"
    if not out_dir.exists():
        return JSONResponse(status_code=404, content={"error": "Job not found"})

    zip_path = job_dir / "out.zip"
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for p in out_dir.rglob("*"):
            if p.is_file():
                z.write(p, arcname=p.relative_to(out_dir).as_posix())

    return FileResponse(
        path=str(zip_path),
        filename=f"md_{job_id}.zip",
        media_type="application/zip",
    )

