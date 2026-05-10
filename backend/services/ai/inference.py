from pathlib import Path

import torch

from backend.services.ai.featurization import mol_to_graph
from backend.services.ai.train_esol import ESOLGNN


MODEL_PATH = Path("backend/services/ai/checkpoints/esol_gnn.pt")

_model = None
_metadata = None


def load_model():
    global _model, _metadata

    if _model is not None:
        return _model

    if not MODEL_PATH.exists():
        raise FileNotFoundError(
            f"Missing model checkpoint: {MODEL_PATH}. "
            "Run: python -m backend.services.ai.train_esol"
        )

    checkpoint = torch.load(MODEL_PATH, map_location="cpu")

    _model = ESOLGNN(in_channels=checkpoint["in_channels"])
    _model.load_state_dict(checkpoint["model_state_dict"])
    _model.eval()

    _metadata = checkpoint
    return _model


def predict_property(smiles: str) -> float:
    model = load_model()
    data = mol_to_graph(smiles)

    with torch.no_grad():
        pred = model(data)

    return float(pred.item())