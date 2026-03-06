from pathlib import Path

def generate_backward_mapping(itp_file: Path, output_map: Path):

    lines = itp_file.read_text().splitlines()

    beads = []
    atom_mapping = []
    residue = None
    atom_index = 1

    in_atoms = False

    for line in lines:

        stripped = line.strip()

        if stripped.lower().startswith("[atoms]"):
            in_atoms = True
            continue

        if stripped.startswith("[") and not stripped.lower().startswith("[atoms]"):
            in_atoms = False

        if not in_atoms:
            continue

        if not stripped or stripped.startswith(";"):
            continue

        parts = stripped.split()

        if not parts[0].isdigit():
            continue

        bead = parts[4]
        residue = parts[3]

        if bead not in beads:
            beads.append(bead)

        if "atoms:" in line:

            atoms_part = line.split("atoms:")[1]
            atoms = [a.strip() for a in atoms_part.split(",") if a.strip()]

            for atom in atoms:
                atom_mapping.append((atom_index, atom, bead))
                atom_index += 1

    with open(output_map, "w") as f:

        f.write("[molecule]\n")
        f.write(f"{residue}\n\n")

        f.write("[martini]\n")
        f.write(" ".join(beads) + "\n\n")

        f.write("[mapping]\n")
        f.write("charmm36\n\n")

        f.write("[atoms]\n")

        for idx, atom, bead in atom_mapping:
            f.write(f"{idx:3d} {atom:4s} {bead}\n")