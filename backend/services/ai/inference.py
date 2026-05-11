from pathlib import Path

import torch

from backend.services.ai.featurization import mol_to_graph
from backend.services.ai.train_esol import ESOLGNN
from backend.services.ai.train_tox21 import Tox21GNN
from backend.services.ai.train_ld50 import LD50GNN
from backend.services.ai.train_bbbp import BBBPGNN
from backend.services.ai.train_clintox import ClinToxGNN


MODEL_PATH = Path("backend/services/ai/checkpoints/esol_gnn.pt")
TOX21_MODEL_PATH = Path("backend/services/ai/checkpoints/tox21_gnn.pt")
LD50_MODEL_PATH = Path("backend/services/ai/checkpoints/ld50_gnn.pt")
BBBP_MODEL_PATH = Path("backend/services/ai/checkpoints/bbbp_gnn.pt")
CLINTOX_MODEL_PATH = Path("backend/services/ai/checkpoints/clintox_gnn.pt")

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


_tox21_model = None
_tox21_metadata = None


def load_tox21_model():
    global _tox21_model, _tox21_metadata

    if _tox21_model is not None:
        return _tox21_model, _tox21_metadata

    if not TOX21_MODEL_PATH.exists():
        raise FileNotFoundError(
            f"Missing model checkpoint: {TOX21_MODEL_PATH}. "
            "Run: python -m backend.services.ai.train_tox21"
        )

    checkpoint = torch.load(
        TOX21_MODEL_PATH,
        map_location="cpu",
        weights_only=True,
    )

    model = Tox21GNN(
        in_channels=checkpoint["in_channels"],
        out_channels=checkpoint["out_channels"],
    )
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()

    _tox21_model = model
    _tox21_metadata = checkpoint

    return _tox21_model, _tox21_metadata


def predict_tox21(smiles: str) -> dict:
    model, metadata = load_tox21_model()
    data = mol_to_graph(smiles)

    with torch.no_grad():
        logits = model(data)
        probabilities = torch.sigmoid(logits).squeeze(0)

    return {
        task: float(prob)
        for task, prob in zip(metadata["tasks"], probabilities)
    }


_ld50_model = None


def load_ld50_model():
    global _ld50_model

    if _ld50_model is not None:
        return _ld50_model

    if not LD50_MODEL_PATH.exists():
        raise FileNotFoundError(
            f"Missing model checkpoint: {LD50_MODEL_PATH}. "
            "Run: python -m backend.services.ai.train_ld50"
        )

    checkpoint = torch.load(
        LD50_MODEL_PATH,
        map_location="cpu",
        weights_only=True,
    )

    model = LD50GNN(in_channels=checkpoint["in_channels"])
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()

    _ld50_model = model
    return _ld50_model

def predict_ld50(smiles: str) -> float:
    model = load_ld50_model()
    data = mol_to_graph(smiles)

    with torch.no_grad():
        pred = model(data)

    return float(pred.item())

_bbbp_model = None
_bbbp_metadata = None


def load_bbbp_model():
    global _bbbp_model, _bbbp_metadata

    if _bbbp_model is not None:
        return _bbbp_model, _bbbp_metadata

    checkpoint = torch.load(
        BBBP_MODEL_PATH,
        map_location="cpu",
        weights_only=True,
    )

    model = BBBPGNN(
        in_channels=checkpoint["in_channels"],
        out_channels=checkpoint["out_channels"],
    )
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()

    _bbbp_model = model
    _bbbp_metadata = checkpoint

    return _bbbp_model, _bbbp_metadata


def predict_bbbp(smiles: str) -> dict:
    model, metadata = load_bbbp_model()
    data = mol_to_graph(smiles)

    with torch.no_grad():
        logits = model(data)
        prob = torch.sigmoid(logits).item()

    return {
        "BBBP": float(prob)
    }

_clintox_model = None
_clintox_metadata = None


def load_clintox_model():
    global _clintox_model, _clintox_metadata

    if _clintox_model is not None:
        return _clintox_model, _clintox_metadata

    checkpoint = torch.load(
        CLINTOX_MODEL_PATH,
        map_location="cpu",
        weights_only=True,
    )

    model = ClinToxGNN(
        in_channels=checkpoint["in_channels"],
        out_channels=checkpoint["out_channels"],
    )
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()

    _clintox_model = model
    _clintox_metadata = checkpoint

    return _clintox_model, _clintox_metadata


def predict_clintox(smiles: str) -> dict:
    model, metadata = load_clintox_model()
    data = mol_to_graph(smiles)

    with torch.no_grad():
        logits = model(data)
        probs = torch.sigmoid(logits).squeeze(0)

    return {
        task: float(prob)
        for task, prob in zip(metadata["tasks"], probs)
    }

