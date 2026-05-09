def predict_property(smiles: str) -> float:
    # Placeholder. Next step: RDKit -> graph -> PyTorch Geometric model.
    if not smiles.strip():
        raise ValueError("Empty SMILES")

    return 0.0