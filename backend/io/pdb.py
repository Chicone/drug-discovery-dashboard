from pathlib import Path

def save_complex_pdb(
    protein_pdb: Path,
    ligand_pdb: Path,
    output_pdb: Path,
    ligand_chain: str = "Z",
):
    """
    Merge a protein PDB and a ligand PDB into a single complex PDB.

    - Protein is written unchanged
    - Ligand atoms are appended
    - Ligand is assigned to a separate chain
    """

    def _is_atom_line(line: str) -> bool:
        return line.startswith(("ATOM", "HETATM"))

    with open(output_pdb, "w") as out:

        # 1) write protein as-is
        with open(protein_pdb) as f:
            for line in f:
                if line.startswith("END"):
                    continue
                out.write(line)

        out.write("TER\n")

        # 2) write ligand
        with open(ligand_pdb) as f:
            for line in f:
                if not _is_atom_line(line):
                    continue

                # enforce chain ID
                line = (
                    line[:21]
                    + ligand_chain
                    + line[22:]
                )
                out.write(line)

        out.write("END\n")
