import torch

from backend.services.ai.featurization import mol_to_graph
from backend.services.ai.gnn_model import SimpleMolecularGNN


# Load model (untrained for now)
_model = SimpleMolecularGNN()
_model.eval()


def predict_property(smiles: str) -> float:
    # Convert SMILES → graph
    data = mol_to_graph(smiles)

    # Forward pass (no gradients needed)
    with torch.no_grad():
        pred = _model(data)

    # Return scalar
    return float(pred.item())