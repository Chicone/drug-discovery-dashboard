from fastapi import FastAPI, Query
from fastapi.responses import Response
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors, AllChem
from rdkit.Chem import QED

app = FastAPI()


@app.get("/api/properties")
def get_properties(smiles: str = Query(..., description="SMILES string")):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return {"error": "Invalid SMILES"}

    # --- base property computation ---
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

    # ================================================================
    # 1️⃣  SIZE AND COMPLEXITY CHECKS
    # ================================================================
    heavy_atoms = mol.GetNumHeavyAtoms()

    if mw < 150 or heavy_atoms < 10:
        insights.append(
            "Molecule is very small (<150 Da or <10 heavy atoms). "
            "Such compounds (e.g., ethanol, acetone) can be biologically active "
            "but usually lack the size and complexity required for selective, drug-like behaviour."
        )

        if logp < 0:
            insights.append("LogP is low/negative, indicating a polar compound.")
        elif logp > 5:
            insights.append("LogP is unusually high for such a small molecule.")
        props["insights"] = insights
        return props

    # ================================================================
    # 2️⃣  HIGH-MW REGION
    # ================================================================
    if mw > 700:
        insights.append(
            "Molecular weight is high (>700 Da), moving it beyond classical small-molecule space; "
            "oral absorption may require special design strategies."
        )

    # ================================================================
    # 3️⃣  LIPINSKI SUMMARY (apply only for 150–700 Da)
    # ================================================================
    violations = 0
    if mw > 500:
        violations += 1
    if logp > 5:
        violations += 1
    if hba > 10:
        violations += 1
    if hbd > 5:
        violations += 1

    if 150 <= mw <= 700:
        if violations == 0:
            insights.append(
                "Physicochemical properties are consistent with Lipinski's oral drug-like region."
            )
        elif violations == 1:
            insights.append(
                "One Lipinski parameter is outside the usual range; may still be acceptable depending on the compound series."
            )
        else:
            insights.append(
                f"{violations} Lipinski parameters are outside the usual range; "
                "consider optimization for oral bioavailability."
            )
    else:
        insights.append(
            "Compound lies outside the size range where Lipinski's rules are typically informative."
        )

    # ================================================================
    # 4️⃣  POLARITY / PERMEABILITY
    # ================================================================
    if tpsa > 140:
        insights.append("TPSA is high (>140 Å²), associated with poor passive permeability.")
    elif 60 <= tpsa <= 140:
        insights.append("TPSA is moderate (60–140 Å²), often acceptable for oral drugs.")
    else:
        insights.append("TPSA is low (<60 Å²), generally favourable for permeability.")

    # ================================================================
    # 5️⃣  LIPOPHILICITY (LogP)
    # ================================================================
    if logp < 0:
        insights.append("LogP is low/negative; compound is polar and may need transporters.")
    elif 0 <= logp <= 5:
        insights.append("LogP is within the broadly acceptable 0–5 range.")
    else:
        insights.append("LogP is high (>5), which may impair solubility and increase off-target risk.")

    # ================================================================
    # 6️⃣  HYDROGEN BONDING
    # ================================================================
    if hba > 10 or hbd > 5:
        insights.append(
            "Hydrogen bond donors/acceptors exceed typical Lipinski limits; "
            "this may reduce membrane permeability."
        )
    else:
        insights.append("Hydrogen bonding parameters are within typical drug-like limits.")

    # ================================================================
    # 7️⃣  FLEXIBILITY
    # ================================================================
    if rb > 10:
        insights.append(
            "Rotatable bond count is high (>10), which can reduce bioavailability and binding selectivity."
        )
    else:
        insights.append("Molecular flexibility is within a commonly acceptable range.")

    props["insights"] = insights
    return props


# ================================================================
# 3D COORDINATE ENDPOINT
# ================================================================
@app.get("/api/mol3d")
def mol3d(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return {"error": "Invalid SMILES"}
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
    mol_block = Chem.MolToMolBlock(mol)
    return Response(mol_block, media_type="text/plain")
