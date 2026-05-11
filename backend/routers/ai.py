from fastapi import APIRouter
from pydantic import BaseModel

from backend.services.ai.inference import (
    predict_property,
    predict_tox21,
    predict_ld50,
    predict_bbbp,
    predict_clintox,
)

router = APIRouter()

class GNNPredictRequest(BaseModel):
    smiles: str

@router.post("/gnn/predict")
def predict_gnn(req: GNNPredictRequest):
    prediction = predict_property(req.smiles)

    return {
        "smiles": req.smiles,
        "model": "esol_gnn",
        "prediction": prediction,
    }

@router.post("/tox21/predict")
def predict_tox21_endpoint(req: GNNPredictRequest):
    predictions = predict_tox21(req.smiles)

    max_prob = max(predictions.values())

    if max_prob >= 0.7:
        risk = "high"
    elif max_prob >= 0.4:
        risk = "medium"
    else:
        risk = "low"

    return {
        "smiles": req.smiles,
        "model": "tox21_gnn",
        "prediction_type": "multi-label toxicity classification",
        "overall_risk": risk,
        "max_probability": max_prob,
        "toxicity_probabilities": predictions,
    }

@router.post("/ld50/predict")
def predict_ld50_endpoint(req: GNNPredictRequest):
    prediction = predict_ld50(req.smiles)

    return {
        "smiles": req.smiles,
        "model": "ld50_gnn",
        "prediction_type": "acute toxicity regression",
        "prediction": prediction,
        "unit": "log(1/(mol/kg))",
    }

@router.post("/bbbp/predict")
def predict_bbbp_endpoint(req: GNNPredictRequest):
    predictions = predict_bbbp(req.smiles)
    probability = predictions["BBBP"]

    label = "likely permeable" if probability >= 0.5 else "likely non-permeable"

    return {
        "smiles": req.smiles,
        "model": "bbbp_gnn",
        "prediction_type": "blood-brain barrier permeability classification",
        "probability": probability,
        "label": label,
    }


@router.post("/clintox/predict")
def predict_clintox_endpoint(req: GNNPredictRequest):
    predictions = predict_clintox(req.smiles)

    ct_tox = predictions.get("CT_TOX", 0.0)

    if ct_tox >= 0.7:
        risk = "high"
    elif ct_tox >= 0.4:
        risk = "medium"
    else:
        risk = "low"

    return {
        "smiles": req.smiles,
        "model": "clintox_gnn",
        "prediction_type": "clinical toxicity classification",
        "clinical_toxicity_signal": risk,
        "ct_tox_probability": ct_tox,
        "probabilities": predictions,
    }