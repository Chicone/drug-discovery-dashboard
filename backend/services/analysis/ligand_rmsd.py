from pathlib import Path
import mdtraj as md
import numpy as np

def detect_backbone_atoms(traj):
    """Detect whether to use CA (atomistic) or BB (Martini)."""
    top = traj.topology

    ca = top.select("name CA")
    if len(ca) > 0:
        return ca

    bb = top.select("name BB")
    if len(bb) > 0:
        return bb

    raise ValueError("No CA or BB atoms found for protein alignment.")

def load_ligand_residue(system_gro: Path) -> str:
    """
    Auto-detect ligand residue name.
    Currently ORT for your Martini ligands.
    """
    with open(system_gro, "r") as f:
        for line in f:
            if "ORT" in line:
                return "ORT"
    raise ValueError("Ligand residue 'ORT' not found in system.gro")

def compute_ligand_rmsd_single(job_dir: Path, ligand_resname="ORT"):
    out_dir = job_dir / "out"

    for candidate in ["md.gro", "npt.gro", "em.gro"]:
        system = out_dir / candidate
        if system.exists():
            break

    traj = None

    preferred = [
        "md_fit.xtc",
        "md_centered.xtc",
        "md_nojump.xtc",
        "md_whole.xtc",
        "md.xtc",
    ]

    for name in preferred:
        p = out_dir / name
        if p.exists() and p.stat().st_size > 0:
            traj = p
            print(f"[RMSD] Using processed trajectory: {p.name}")
            break

    if traj is None:
        import re

        part_files = list(out_dir.glob("md.part*.xtc"))

        if not part_files:
            raise FileNotFoundError(
                f"No trajectory found in {out_dir} "
                f"({', '.join(preferred)} or md.partXXXX.xtc)"
            )

        def part_index(path):
            m = re.search(r"\.part(\d+)\.xtc", path.name)
            return int(m.group(1)) if m else -1

        traj = max(part_files, key=part_index)
        print(f"[RMSD] Using continuation trajectory: {traj.name}")

    print("SYSTEM:", system)
    print("SYSTEM EXISTS:", system.exists())
    print("TRAJ:", traj)
    print("TRAJ EXISTS:", traj.exists())

    if not system.exists():
        return {"error": f"Missing file: {system}"}
    if not traj.exists():
        return {"error": f"Missing file: {traj}"}


    # Load full trajectory with topology
    try:
        t = md.load(str(traj), top=str(system))
    except Exception as e:
        return {
            "error": f"Failed to load trajectory with topology.\n"
                     f"traj={traj}\n"
                     f"top={system}\n"
                     f"{type(e).__name__}: {e}"
        }
    print("DEBUG: md.txc not found", traj)

    # t = md.load(str(traj), top=str(system))
    # Reference = first frame of trajectory (Trajectory object)
    ref = t[0:1]
    ref.xyz = ref.xyz.copy()

    # Backbone atoms for alignment
    backbone = detect_backbone_atoms(t)

    # Align entire trajectory using the frozen protein reference
    t.superpose(t, 0, atom_indices=backbone)

    ligand_idx = t.topology.select(f"resname {ligand_resname}")
    if len(ligand_idx) == 0:
        return {"error": f"Ligand '{ligand_resname}' not found"}

    print("DEBUG:",
          "frames =", len(t),
          "ligand atoms =", len(ligand_idx),
          "times =", t.time[:5])

    # RMSD vs frame 0 of same trajectory
    import numpy as np

    coords0 = t.xyz[0, ligand_idx, :]  # (N, 3)
    diff = (t.xyz[:, ligand_idx, :] - coords0[None, :, :]) * 10.0
    rmsd = np.sqrt((diff * diff).sum(axis=2).mean(axis=1))

    # Extract true per-frame times from MDTraj (ps)
    times_ps = t.time.astype(float).tolist()

    # MD timestep metadata
    timestep_val = getattr(t, "timestep", None)
    timestep_val = float(timestep_val) if timestep_val is not None else None


    return {
        "job_id": job_dir.name,
        "frames": int(len(rmsd)),
        "timestep_ps": timestep_val,
        "times_ps": times_ps,  # <<< NEW
        "rmsd": rmsd.tolist(),
    }


def compute_ligand_rmsd_multi(job_dirs: list[Path]):
    """
    Compute ligand RMSD for several replicas and aggregate statistics.
    """

    results = []
    rmsd_arrays = []

    for jd in job_dirs:
        r = compute_ligand_rmsd_single(jd)
        results.append(r)

        if "rmsd" in r:
            rmsd_arrays.append(np.array(r["rmsd"]))

    if len(rmsd_arrays) == 0:
        return {"error": "No successful RMSD computations"}

    # Normalize lengths (pad shorter trajectories)
    max_len = max(arr.size for arr in rmsd_arrays)
    padded = np.array([
        np.pad(arr, (0, max_len - arr.size), "edge")
        for arr in rmsd_arrays
    ])

    return {
        "replicas": results,
        "aggregate": {
            "mean": padded.mean(axis=0).tolist(),
            "std": padded.std(axis=0).tolist(),
            "min": padded.min(axis=0).tolist(),
            "max": padded.max(axis=0).tolist(),
        },
    }
