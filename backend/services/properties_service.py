# backend/services/properties_service.py

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
from rdkit.Chem import QED


def compute_properties(smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return {"error": "Invalid SMILES"}

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hba = Lipinski.NumHAcceptors(mol)
    hbd = Lipinski.NumHDonors(mol)
    rb = Lipinski.NumRotatableBonds(mol)
    qed_score = round(QED.qed(mol), 3)

    props = {
        "MolecularWeight": mw,
        "LogP": logp,
        "TPSA": tpsa,
        "NumHAcceptors": hba,
        "NumHDonors": hbd,
        "RotatableBonds": rb,
        "QED": qed_score,
    }

    insights = []
    heavy_atoms = mol.GetNumHeavyAtoms()

    # -------------------------
    # SMALL MOLECULE CHECK
    # -------------------------
    if mw < 150 or heavy_atoms < 10:
        insights.append(
            "Molecule is very small (<150 Da or <10 heavy atoms). "
            "Often too small for selective drug-like behavior."
        )

        if logp < 0:
            insights.append("LogP is low/negative (polar compound).")
        elif logp > 5:
            insights.append("LogP unusually high for such a small molecule.")

        props["insights"] = insights
        return props

    # -------------------------
    # HIGH MW
    # -------------------------
    if mw > 700:
        insights.append(
            "High molecular weight (>700 Da); beyond classical small-molecule space."
        )

    # -------------------------
    # LIPINSKI
    # -------------------------
    violations = sum(
        [
            mw > 500,
            logp > 5,
            hba > 10,
            hbd > 5,
        ]
    )

    if 150 <= mw <= 700:
        if violations == 0:
            insights.append("Lipinski OK — no violations.")
        elif violations == 1:
            insights.append("1 Lipinski violation (still may be OK).")
        else:
            insights.append(f"{violations} Lipinski violations — may reduce oral absorption.")
    else:
        insights.append("Lipinski rules not informative for this MW range.")

    # -------------------------
    # TPSA (polarity)
    # -------------------------
    if tpsa > 140:
        insights.append("High polarity (>140 Å²) — poor passive permeability.")
    elif 60 <= tpsa <= 140:
        insights.append("Moderate polarity (60–140 Å²).")
    else:
        insights.append("Low polarity (<60 Å²) — good permeability.")

    # -------------------------
    # LOGP
    # -------------------------
    if logp < 0:
        insights.append("Very polar (logP < 0).")
    elif logp <= 5:
        insights.append("LogP in acceptable range (0–5).")
    else:
        insights.append("High lipophilicity (logP > 5).")

    # -------------------------
    # H-Donors/Acceptors
    # -------------------------
    if hba > 10 or hbd > 5:
        insights.append("High number of HBA/HBD — may reduce permeability.")
    else:
        insights.append("HBA/HBD within typical drug-like range.")

    # -------------------------
    # FLEXIBILITY
    # -------------------------
    if rb > 10:
        insights.append("High flexibility (>10 rotatable bonds) can reduce bioavailability.")
    else:
        insights.append("Acceptable molecular flexibility.")

    props["insights"] = insights
    return props
