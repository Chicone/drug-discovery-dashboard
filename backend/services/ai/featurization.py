from rdkit import Chem
import torch
from torch_geometric.data import Data


def mol_to_graph(smiles: str) -> Data:
    # Convert SMILES → RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    # Build node features (per atom)
    atom_features = []
    for atom in mol.GetAtoms():
        atom_features.append([
            atom.GetAtomicNum(),      # element type
            atom.GetDegree(),         # number of bonds
            atom.GetFormalCharge(),   # charge
            int(atom.GetIsAromatic()) # aromatic flag
        ])

    x = torch.tensor(atom_features, dtype=torch.float)

    # Build edges (bond connections, bidirectional)
    edges = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        edges.append([i, j])
        edges.append([j, i])

    # Handle molecules with no bonds
    if edges:
        edge_index = torch.tensor(edges, dtype=torch.long).t().contiguous()
    else:
        edge_index = torch.empty((2, 0), dtype=torch.long)

    # PyG graph object
    return Data(x=x, edge_index=edge_index)