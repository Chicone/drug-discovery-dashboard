from torch_geometric.utils.smiles import from_smiles


def mol_to_graph(smiles: str):
    # Use PyG's SMILES featurizer so inference matches MoleculeNet training.
    data = from_smiles(smiles)

    if data.x is None or data.x.numel() == 0:
        raise ValueError(f"Invalid SMILES: {smiles}")

    return data