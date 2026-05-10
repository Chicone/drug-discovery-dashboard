from fastapi import APIRouter
from pydantic import BaseModel

from backend.services.ai.inference import predict_property


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