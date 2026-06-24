# backend/pipelines/ligand_cg_builder.py
#
# Robust Martini-3 ligand CG builder using Auto-MartiniM3.
#
# Inputs:
#   - SMILES: chemistry template (bond orders, aromaticity)
#   - Docked complex PDB: ligand coordinates (HETATM pose)
#
# Outputs:
#   - orthosteric_cg.pdb   (CG beads at the docked pose)
#   - Orthosteric.itp      (Martini 3 topology for that ligand)
#
# Requirements (same conda env as your pipeline):
#   - rdkit
#   - auto_martiniM3  (python -m auto_martiniM3 must work)
#
# Notes:
#   - Auto-MartiniM3 restricts names to <= 5 characters (README).
#     We'll use "ORT" by default.
#
from __future__ import annotations

import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem



@dataclass
class AtomRec:
    pdb_serial: int
    name: str
    resname: str
    x: float
    y: float
    z: float
    element: str


def _run(cmd: List[str], cwd: Path) -> None:
    print("\n$ " + " ".join(cmd))
    subprocess.run(cmd, cwd=str(cwd), check=True)


def _guess_element(atom_name: str) -> str:
    # Very small heuristic; PDB element column is often missing in docked files.
    # "CL", "BR" etc must be handled.
    s = re.sub(r"[^A-Za-z]", "", atom_name).upper()
    if s.startswith(("CL", "BR")):
        return s[:2]
    return s[:1] if s else "C"


def load_ligand_atoms_from_complex(
    complex_pdb: Path,
    resname: Optional[str] = None,
) -> List[AtomRec]:
    atoms: List[AtomRec] = []

    for line in complex_pdb.read_text().splitlines():
        if not line.startswith("HETATM"):
            continue

        rn = line[17:20].strip()
        if resname is not None and rn != resname:
            continue

        serial = int(line[6:11])
        name = line[12:16].strip()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

        elem = line[76:78].strip().upper()
        if elem == "H":   # skip hydrogens
            continue

        elem = line[76:78].strip().upper()
        if not elem:
            elem = _guess_element(name)

        atoms.append(AtomRec(serial, name, rn, x, y, z, elem))

    if not atoms:
        raise RuntimeError(
            f"No ligand HETATM records found in {complex_pdb} "
            f"(resname filter={resname!r})."
        )
    return atoms




def rdkit_mol_from_smiles_and_pose_mapped(
    smiles: str,
    complex_pdb: Path,
    ligand_resname_in_pdb: str | None = None,
) -> Chem.Mol:
    """
    Map docked PDB coordinates onto the SMILES atom order without relying on
    RDKit bond inference or template bond-order transfer.

    Strategy:
      - parse heavy atoms + coords from PDB
      - infer a PDB bond graph using covalent radii + distance thresholds
      - use SMILES bond graph as ground truth
      - find an isomorphism mapping (SMILES atom -> PDB atom)
      - copy PDB coords onto a SMILES-derived RDKit mol (SMILES atom order)
    """

    # -----------------------------
    # A) Read PDB heavy atoms
    # -----------------------------
    lig_res = ligand_resname_in_pdb or "UNL"

    pdb_atoms = load_ligand_atoms_from_complex(complex_pdb, resname=lig_res)
    # load_ligand_atoms_from_complex already skips H
    n_pdb = len(pdb_atoms)
    if n_pdb == 0:
        raise RuntimeError(
            f"No heavy atoms found for ligand resname={lig_res!r} "
            f"in {complex_pdb}"
        )

    pdb_elems = [a.element.upper() for a in pdb_atoms]
    pdb_xyz = np.asarray([[a.x, a.y, a.z] for a in pdb_atoms], dtype=float)

    # -----------------------------
    # B) Build SMILES heavy mol
    # -----------------------------
    smi = Chem.MolFromSmiles(smiles)
    if smi is None:
        raise RuntimeError(f"Invalid SMILES: {smiles!r}")

    smi = Chem.RemoveHs(smi)
    n_smi = smi.GetNumAtoms()
    if n_smi != n_pdb:
        raise RuntimeError(
            "SMILES and docked pose have different heavy atom counts: "
            f"SMILES={n_smi} PDB={n_pdb}."
        )

    smi_elems = [a.GetSymbol().upper() for a in smi.GetAtoms()]

    # -----------------------------
    # C) Build adjacency graphs
    # -----------------------------
    smi_adj = _adjacency_from_rdkit_mol(smi)
    pdb_adj = _adjacency_from_geometry(pdb_elems, pdb_xyz)

    # Strong early sanity check: degree multisets per element should match
    if not _quick_degree_signature_ok(smi_elems, smi_adj, pdb_elems, pdb_adj):
        raise RuntimeError(
            "PDB-inferred bond graph signature does not match SMILES "
            "(element/degree mismatch). This usually means the PDB coords "
            "are missing atoms, have wrong elements, or the bond inference "
            "threshold is off for this ligand."
        )

    # -----------------------------
    # D) Find mapping: SMILES -> PDB
    # -----------------------------
    mapping = _find_isomorphism_mapping(
        smi_elems=smi_elems,
        smi_adj=smi_adj,
        pdb_elems=pdb_elems,
        pdb_adj=pdb_adj,
        pdb_xyz=pdb_xyz,
    )

    # mapping[i_smi] = j_pdb
    if mapping is None:
        raise RuntimeError(
            "Failed to find a consistent atom mapping from SMILES to PDB. "
            "This means either:\n"
            "  - the inferred PDB bond graph is wrong, or\n"
            "  - the SMILES is not the same connectivity as the docked pose.\n"
            "We can tighten/loosen the bond threshold or add extra pruning."
        )

    # -----------------------------
    # E) Create output mol (SMILES chemistry, docked coords)
    # -----------------------------
    out = Chem.Mol(smi)
    conf = Chem.Conformer(out.GetNumAtoms())

    for i_smi, j_pdb in enumerate(mapping):
        x, y, z = pdb_xyz[int(j_pdb)]
        conf.SetAtomPosition(i_smi, Chem.rdGeometry.Point3D(x, y, z))

    out.RemoveAllConformers()
    out.AddConformer(conf, assignId=True)
    return out


def _adjacency_from_rdkit_mol(m: Chem.Mol) -> list[set[int]]:
    n = m.GetNumAtoms()
    adj = [set() for _ in range(n)]
    for b in m.GetBonds():
        i = b.GetBeginAtomIdx()
        j = b.GetEndAtomIdx()
        adj[i].add(j)
        adj[j].add(i)
    return adj


def _adjacency_from_geometry(
    elems: list[str],
    xyz: np.ndarray,
    *,
    fudge: float = 0.45,
    max_bond: float = 1.85,
) -> list[set[int]]:
    """
    Infer covalent bonds from distances.

    - fudge: extra Angstrom added to sum of covalent radii
    - max_bond: hard cap to prevent absurd long "bonds"
    """
    cov = {
        "H": 0.31,
        "C": 0.76,
        "N": 0.71,
        "O": 0.66,
        "F": 0.57,
        "P": 1.07,
        "S": 1.05,
        "CL": 1.02,
        "BR": 1.20,
        "I": 1.39,
    }

    n = len(elems)
    adj = [set() for _ in range(n)]

    def rad(e: str) -> float:
        e = e.upper()
        if e in cov:
            return cov[e]
        # handle two-letter halogens if element column got truncated
        if e.startswith("CL"):
            return cov["CL"]
        if e.startswith("BR"):
            return cov["BR"]
        return 0.77  # fallback ~ carbon

    for i in range(n):
        ri = rad(elems[i])
        xi = xyz[i]
        for j in range(i + 1, n):
            rj = rad(elems[j])
            d = float(np.linalg.norm(xi - xyz[j]))
            thr = ri + rj + fudge
            if d <= thr and d <= max_bond:
                adj[i].add(j)
                adj[j].add(i)

    return adj


def _quick_degree_signature_ok(
    smi_elems: list[str],
    smi_adj: list[set[int]],
    pdb_elems: list[str],
    pdb_adj: list[set[int]],
) -> bool:
    """
    Compare multiset of (element, degree) between SMILES and inferred PDB graph.
    If this fails, mapping is very unlikely to work.
    """
    def sig(elems, adj):
        s = []
        for i, e in enumerate(elems):
            s.append((e, len(adj[i])))
        return sorted(s)

    return sig(smi_elems, smi_adj) == sig(pdb_elems, pdb_adj)


def _find_isomorphism_mapping(
    *,
    smi_elems: list[str],
    smi_adj: list[set[int]],
    pdb_elems: list[str],
    pdb_adj: list[set[int]],
    pdb_xyz: np.ndarray,
) -> list[int] | None:
    """
    Backtracking graph isomorphism with pruning:
      - element match
      - degree match
      - neighbor consistency for already assigned atoms
      - optional geometry tie-break (helps with symmetric cases)
    """
    n = len(smi_elems)

    # Candidate pdb atoms for each smi atom
    smi_deg = [len(smi_adj[i]) for i in range(n)]
    pdb_deg = [len(pdb_adj[j]) for j in range(n)]

    candidates: list[list[int]] = []
    for i in range(n):
        e = smi_elems[i]
        di = smi_deg[i]
        c = [j for j in range(n)
             if pdb_elems[j] == e and pdb_deg[j] == di]
        if not c:
            return None
        candidates.append(c)

    # Order SMILES atoms: hardest first (few candidates, high degree)
    order = list(range(n))
    order.sort(key=lambda i: (len(candidates[i]), -smi_deg[i]))

    mapping = [-1] * n
    used = [False] * n

    # Precompute PDB pairwise distances for symmetry-breaking
    pdb_dist = _pairwise_dist_matrix(pdb_xyz)

    # Simple SMILES topological distances (for tie-break only)
    smi_topo = _topo_dist_matrix(smi_adj)

    def ok_partial(i_smi: int, j_pdb: int) -> bool:
        # Check consistency with already-mapped neighbors
        for nb in smi_adj[i_smi]:
            jb = mapping[nb]
            if jb == -1:
                continue
            if jb not in pdb_adj[j_pdb]:
                return False

        # Optional: preserve some distance relations to break symmetries
        # For each already-mapped atom k, compare relative distances
        for k in range(n):
            jk = mapping[k]
            if jk == -1 or k == i_smi:
                continue

            # If SMILES graph says they are far in topology,
            # they should not be covalently close in PDB.
            if smi_topo[i_smi][k] >= 4:
                if pdb_dist[j_pdb][jk] < 1.25:
                    return False

        return True

    def dfs(pos: int) -> bool:
        if pos == n:
            return True

        i = order[pos]

        # Try candidates, but in a geometry-informed order to stabilize
        cand = candidates[i]
        cand_sorted = sorted(
            cand,
            key=lambda j: _geom_score_for_choice(i, j, mapping,
                                                 pdb_dist, smi_topo),
        )

        for j in cand_sorted:
            if used[j]:
                continue
            if not ok_partial(i, j):
                continue

            mapping[i] = j
            used[j] = True

            if dfs(pos + 1):
                return True

            mapping[i] = -1
            used[j] = False

        return False

    if dfs(0):
        # Return in SMILES atom order
        return mapping

    return None


def _pairwise_dist_matrix(xyz: np.ndarray) -> np.ndarray:
    n = xyz.shape[0]
    d = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            dij = float(np.linalg.norm(xyz[i] - xyz[j]))
            d[i, j] = dij
            d[j, i] = dij
    return d


def _topo_dist_matrix(adj: list[set[int]]) -> list[list[int]]:
    """
    All-pairs shortest path lengths in an unweighted graph.
    """
    n = len(adj)
    dist = [[10**9] * n for _ in range(n)]
    for i in range(n):
        dist[i][i] = 0
        q = [i]
        head = 0
        while head < len(q):
            u = q[head]
            head += 1
            for v in adj[u]:
                if dist[i][v] > dist[i][u] + 1:
                    dist[i][v] = dist[i][u] + 1
                    q.append(v)
    return dist


def _geom_score_for_choice(
    i_smi: int,
    j_pdb: int,
    mapping: list[int],
    pdb_dist: np.ndarray,
    smi_topo: list[list[int]],
) -> float:
    """
    Lower is better. Used only to order candidate tries.
    Penalize choices that make topologically close atoms far in PDB.
    """
    score = 0.0
    for k_smi, k_pdb in enumerate(mapping):
        if k_pdb == -1:
            continue

        td = smi_topo[i_smi][k_smi]
        gd = pdb_dist[j_pdb][k_pdb]

        # If topologically adjacent, prefer physically close
        if td == 1:
            score += abs(gd - 1.40)
        # If far in graph, prefer not super close
        elif td >= 4:
            score += max(0.0, 1.20 - gd) * 5.0

    return score





def write_sdf_with_pose(mol: Chem.Mol, out_sdf: Path) -> None:
    w = Chem.SDWriter(str(out_sdf))
    if w is None:
        raise RuntimeError(f"Failed to open SDF writer for {out_sdf}")
    w.write(mol)
    w.close()


def gro_to_pdb(gro_path: Path, pdb_path: Path,
               resname: str = "ORT", chain_id: str = "L", resid: int = 1) -> None:
    """
    Strict, PDB-spec-compliant .gro -> .pdb converter for Martini CG.
    Produces fully valid fixed-column ATOM records.
    """

    lines = gro_path.read_text().splitlines()
    if len(lines) < 3:
        raise RuntimeError(f"Invalid .gro: {gro_path}")

    # Number of atoms
    try:
        n_atoms = int(lines[1].strip())
    except ValueError:
        raise RuntimeError(f"Invalid atom count line in {gro_path}")

    atom_lines = lines[2:2 + n_atoms]

    out = []
    serial = 1

    for l in atom_lines:
        # .gro format: resnr(5) resname(5) atomname(5) atomnr(5) x y z
        atomname = l[10:15].strip()

        # Convert nm → Å
        x = float(l[20:28]) * 10.0
        y = float(l[28:36]) * 10.0
        z = float(l[36:44]) * 10.0

        # Strict PDB ATOM format
        line = (
            f"ATOM  "
            f"{serial:5d} "
            f"{atomname:<4s}"
            f"{resname:>4s}"
            f"{chain_id:>2s}"
            f"{resid:4d}    "
            f"{x:8.3f}"
            f"{y:8.3f}"
            f"{z:8.3f}"
            f"{1.00:6.2f}"
            f"{0.00:6.2f}"
            f"{'':10s}"
            f"{'C':>2s}"
        )
        out.append(line)
        serial += 1

    out.append("END")
    pdb_path.write_text("\n".join(out) + "\n")



def build_ligand(
    smiles: str,
    complex_pdb: Path,
    cg_pdb_out: Path,
    itp_out: Path,
    *,
    ligand_resname_in_pdb: Optional[str] = None,
    cg_name: str = "ORT",
    workdir: Optional[Path] = None,
) -> None:
    """
    Build Martini-3 CG ligand topology + coordinates at the docked pose.

    Steps:
      1) Build RDKit mol in SMILES atom order with DOCKED coordinates (Å)
      2) Write SDF (Å) for auto_martiniM3
      3) Run auto_martiniM3 -> produces CG .gro (nm) + .itp
      4) Compute rigid transform from RDKit heavy atoms (Å) -> docked heavy atoms (Å)
      5) Apply that transform to CG beads (nm -> Å -> transform -> nm)
      6) Export aligned CG PDB
    """
    if len(cg_name) > 5:
        raise ValueError("cg_name must be <= 5 characters (Auto-MartiniM3).")

    complex_pdb = Path(complex_pdb)
    cg_pdb_out = Path(cg_pdb_out)
    itp_out = Path(itp_out)

    if workdir is None:
        workdir = itp_out.parent
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    # 1) RDKit mol: SMILES chemistry + docked coords (Å), in SMILES atom order
    mol = rdkit_mol_from_smiles_and_pose_mapped(
        smiles=smiles,
        complex_pdb=complex_pdb,
        ligand_resname_in_pdb=ligand_resname_in_pdb,
    )

    # 2) Write SDF (auto_martiniM3 reads coords from here)
    sdf_path = workdir / f"{cg_name}_pose.sdf"
    write_sdf_with_pose(mol, sdf_path)

    # 3) Run Auto-MartiniM3 -> {cg_name}.itp and {cg_name}.gro (nm)
    _run(
        [
            sys.executable, "-m", "auto_martiniM3",
            "--sdf", str(sdf_path),
            "--mol", cg_name,
        ],
        cwd=workdir,
    )

    gen_itp = workdir / f"{cg_name}.itp"
    gen_gro = workdir / f"{cg_name}.gro"

    if not gen_itp.exists():
        raise RuntimeError(f"Auto-MartiniM3 did not produce {gen_itp}")
    if not gen_gro.exists():
        raise RuntimeError(f"Auto-MartiniM3 did not produce {gen_gro}")

    # 4) Save topology where your pipeline expects it
    itp_out.write_text(gen_itp.read_text(), encoding="utf-8")

    relax_martini_ligand_itp(itp_out)

    # ==============================================================
    # 5) ALIGN CG BEADS TO THE DOCKED POSE (robust, no atom ordering)
    # ==============================================================

    # A) Parse bead -> atom groups from the .itp produced by Auto-MartiniM3
    bead_groups = parse_bead_groups_from_itp(gen_itp)

    # B) Target bead positions (Å): centroids of the bead-mapped atoms
    #    from the RDKit mol (already at docked pose coordinates)
    target_beads_A = bead_centroids_from_rdkit(mol, bead_groups)

    # C) Current bead positions (Å): from generated GRO (nm -> Å)
    curr_beads_nm = gro_bead_coords_nm(gen_gro)
    curr_beads_A = curr_beads_nm * 10.0

    if target_beads_A.shape != curr_beads_A.shape:
        raise RuntimeError(
            "Bead alignment failed: target beads "
            f"{target_beads_A.shape} vs GRO beads {curr_beads_A.shape}. "
            "This usually means bead-group parsing failed or Auto-MartiniM3 "
            "changed its output format."
        )

    # D) Compute rigid transform mapping current beads -> target beads (Å -> Å)
    #    kabsch_rigid_transform(P, Q) returns R,t such that R@Q + t ~= P
    R, t = kabsch_rigid_transform(target_beads_A, curr_beads_A)

    # E) Apply transform to GRO coordinates (nm -> Å -> transform -> nm)
    aligned_gro = workdir / f"{cg_name}_aligned.gro"
    apply_rigid_to_gro(gen_gro, aligned_gro, R, t)

    # 6) Export aligned CG PDB
    gro_to_pdb(aligned_gro, cg_pdb_out, resname=cg_name, chain_id="L")

    print(f"[ligand] wrote: {itp_out}")
    print(f"[ligand] wrote: {cg_pdb_out}")
    print(f"[ligand] wrote: {aligned_gro}")

    # 6) Export aligned CG PDB
    gro_to_pdb(aligned_gro, cg_pdb_out, resname=cg_name, chain_id="L")

    print(f"[ligand] wrote: {itp_out}")
    print(f"[ligand] wrote: {cg_pdb_out}")



def parse_bead_groups_from_itp(itp_path: Path) -> List[List[int]]:
    """
    Parse Auto-MartiniM3 .itp and extract bead -> atom-index groups from
    the '; atoms:' comment in the [atoms] section.

    Returns: groups[bead_index_0based] = [atom_indices]
    """
    txt = itp_path.read_text(encoding="utf-8", errors="replace")
    lines = txt.splitlines()

    in_atoms = False
    groups: List[List[int]] = []

    for line in lines:
        s = line.strip()

        if s.startswith("["):
            in_atoms = (s.lower() == "[atoms]")
            continue

        if not in_atoms:
            continue

        if not s or s.startswith(";"):
            continue

        # Example line (trimmed):
        # 1  TN3  1 ORTH N01 1 0 36 ; C#N ; atoms: N0, C1, ; ...
        if "; atoms:" not in line:
            continue

        # Left part has bead id; we only need ordering, but keep it robust.
        left = line.split(";", 1)[0].strip()
        parts = left.split()
        if len(parts) < 1:
            continue

        try:
            bead_id = int(parts[0])  # 1-based in itp
        except ValueError:
            continue

        # Extract tokens after 'atoms:' until the next ';' (or end).
        rhs = line.split("; atoms:", 1)[1]
        rhs = rhs.split(";", 1)[0]

        # Tokens like "N0, C1, C14, N15"
        toks = [t.strip() for t in rhs.replace(",", " ").split()]
        idxs: List[int] = []
        for t in toks:
            # Take trailing integer after element letters: C17 -> 17
            num = ""
            for ch in t[::-1]:
                if ch.isdigit():
                    num = ch + num
                else:
                    break
            if num:
                idxs.append(int(num))

        # Ensure list size
        while len(groups) < bead_id:
            groups.append([])

        groups[bead_id - 1] = idxs

    if not groups or any(len(g) == 0 for g in groups):
        raise RuntimeError(
            f"Failed to parse bead groups from {itp_path}. "
            f"Expected '; atoms:' comments in [atoms]."
        )

    return groups


def _detect_0_or_1_based(atom_indices: List[int], n_atoms: int) -> int:
    """
    Decide whether indices are 0-based or 1-based.
    Returns offset to add to raw index to get 0-based (0 or -1).
    """
    mx = max(atom_indices)
    mn = min(atom_indices)

    # If any 0 appears, it's almost certainly 0-based.
    if mn == 0:
        return 0

    # If max equals n_atoms, likely 1-based (1..n_atoms)
    if mx == n_atoms:
        return -1

    # If max equals n_atoms-1, likely 0-based
    if mx == n_atoms - 1:
        return 0

    # Fallback: assume 0-based (safer than shifting into negatives)
    return 0


def bead_centroids_from_rdkit(
    mol: Chem.Mol,
    bead_groups: List[List[int]],
) -> np.ndarray:
    """
    Compute reference bead centroids (Angstrom) from the RDKit conformer,
    using bead->atom groups extracted from the .itp comments.
    """
    conf = mol.GetConformer()
    n_atoms = mol.GetNumAtoms()

    flat = [i for g in bead_groups for i in g]
    offset = _detect_0_or_1_based(flat, n_atoms)

    cents = []
    for g in bead_groups:
        pts = []
        for raw_idx in g:
            idx = raw_idx + offset
            if idx < 0 or idx >= n_atoms:
                raise RuntimeError(
                    f"Atom index {raw_idx} (offset {offset}) out of range "
                    f"for RDKit mol with {n_atoms} atoms."
                )
            p = conf.GetAtomPosition(idx)
            pts.append([p.x, p.y, p.z])
        pts = np.asarray(pts, dtype=float)
        cents.append(pts.mean(axis=0))

    return np.asarray(cents, dtype=float)


def gro_bead_coords_nm(gro_path: Path) -> np.ndarray:
    """
    Read GRO and return coordinates as (n_beads, 3) in nm.
    """
    lines = gro_path.read_text().splitlines()
    n = int(lines[1].strip())
    atom_lines = lines[2:2 + n]

    xyz = np.zeros((n, 3), dtype=float)
    for i, l in enumerate(atom_lines):
        xyz[i, 0] = float(l[20:28])
        xyz[i, 1] = float(l[28:36])
        xyz[i, 2] = float(l[36:44])

    return xyz


def kabsch_rigid_transform(P: np.ndarray, Q: np.ndarray):
    """
    Find R,t such that (R @ Q.T).T + t best matches P (least squares).
    P, Q are (N,3) arrays in same units.
    """
    if P.shape != Q.shape or P.shape[1] != 3:
        raise ValueError("P and Q must be (N,3) arrays with same shape.")

    Pc = P.mean(axis=0)
    Qc = Q.mean(axis=0)

    X = P - Pc
    Y = Q - Qc

    H = Y.T @ X
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T

    # Fix improper rotation (reflection)
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    t = Pc - (R @ Qc)
    return R, t


def apply_rigid_to_gro(
    gro_in: Path,
    gro_out: Path,
    R: np.ndarray,
    t: np.ndarray,
) -> None:
    """
    Apply rigid transform to GRO coordinates:
      x_A = R @ x_A + t_A   (in Angstrom)
    GRO is nm, so convert nm->A, transform, then A->nm.
    """
    lines = gro_in.read_text().splitlines()
    n = int(lines[1].strip())
    atom_lines = lines[2:2 + n]
    box_line = lines[2 + n] if len(lines) > 2 + n else None

    out_atom_lines = []
    for l in atom_lines:
        x_nm = float(l[20:28])
        y_nm = float(l[28:36])
        z_nm = float(l[36:44])

        xA = np.array([x_nm, y_nm, z_nm], dtype=float) * 10.0
        xA2 = (R @ xA) + t
        x_nm2 = xA2 / 10.0

        out_atom_lines.append(
            f"{l[:20]}"
            f"{x_nm2[0]:8.3f}{x_nm2[1]:8.3f}{x_nm2[2]:8.3f}"
            f"{l[44:]}"
        )

    out_lines = [lines[0], lines[1]] + out_atom_lines
    if box_line is not None:
        out_lines.append(box_line)

    gro_out.write_text("\n".join(out_lines) + "\n")


def debug_compare_smiles_vs_pdb(smiles_mol, pdb_mol):
    print("\n====== DEBUG: SMILES vs PDB comparison ======")

    print(f"SMILES heavy atoms: {smiles_mol.GetNumAtoms()}")
    print(f"PDB heavy atoms:    {pdb_mol.GetNumAtoms()}")

    print("\nSMILES atoms:")
    for i, a in enumerate(smiles_mol.GetAtoms()):
        print(f"  {i:2d}: {a.GetSymbol()}")

    print("\nPDB atoms:")
    for i, a in enumerate(pdb_mol.GetAtoms()):
        print(f"  {i:2d}: {a.GetSymbol()}")

    print("\nTrying substructure match (SMILES -> PDB):")
    m = pdb_mol.GetSubstructMatch(smiles_mol)
    print("Match:", m)

    print("\nTrying reverse match (PDB -> SMILES):")
    m2 = smiles_mol.GetSubstructMatch(pdb_mol)
    print("Reverse match:", m2)

    print("\nSMILES adjacency:")
    for b in smiles_mol.GetBonds():
        print(f" {b.GetBeginAtomIdx()} - {b.GetEndAtomIdx()}  type={b.GetBondType()}")

    print("\nPDB adjacency:")
    for b in pdb_mol.GetBonds():
        print(f" {b.GetBeginAtomIdx()} - {b.GetEndAtomIdx()}  type={b.GetBondType()}")

    print("====== END DEBUG ======\n")


def _force_pdb_element_column(line: str, elem: str) -> str:
    elem = elem.strip().upper()[:2]
    line = line.rstrip("\n")
    if len(line) < 78:
        line = line.ljust(78)
    # element column is 77-78 (0-based 76:78)
    return line[:76] + f"{elem:>2s}" + line[78:] + "\n"


def relax_martini_ligand_itp(itp_path, cap_fc=2500.0):
    lines = Path(itp_path).read_text().splitlines()
    out = []
    in_bonds = False

    for line in lines:
        stripped = line.strip()

        # Enter the bonds section IF the line contains "[bonds]"
        if stripped.lower().startswith("[bonds]"):
            in_bonds = True
            out.append(line)
            continue

        # If we hit a new section, exit bond block
        if in_bonds and stripped.startswith("[") and not stripped.lower().startswith("[bonds]"):
            in_bonds = False
            out.append(line)
            continue

        # If inside [bonds], attempt modification
        if in_bonds:
            # Skip comments
            if stripped.startswith(";") or stripped == "":
                out.append(line)
                continue

            # Tokenize
            toks = stripped.split()
            # Expected at least: i j funct length fc
            if len(toks) >= 5 and toks[0].isdigit() and toks[1].isdigit():
                try:
                    fc = float(toks[4])
                except ValueError:
                    out.append(line)
                    continue

                # Replace if above threshold
                if fc > cap_fc:
                    toks[4] = f"{cap_fc:.2f}"
                    # Write back as valid GROMACS line
                    new_line = "   " + "   ".join(toks)
                    out.append(new_line)
                    continue

        # Default: write unchanged
        out.append(line)

    # Save file
    Path(itp_path).write_text("\n".join(out) + "\n")