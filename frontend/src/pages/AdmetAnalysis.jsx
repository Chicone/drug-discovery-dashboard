import { useState } from "react";
import {
  Paper,
  Typography,
  TextField,
  Button,
  Stack,
  Table,
  TableBody,
  TableCell,
  TableRow,
  Alert,
  MenuItem,
  Chip,
  Divider,
} from "@mui/material";

const API_BASE = "http://localhost:8001";

const exampleMolecules = [
  { name: "Ethanol", smiles: "CCO" },
  { name: "Caffeine", smiles: "Cn1cnc2c1c(=O)n(C)c(=O)n2C" },
  { name: "Aspirin", smiles: "CC(=O)Oc1ccccc1C(=O)O" },
  { name: "Ibuprofen", smiles: "CC(C)Cc1ccc(cc1)C(C)C(=O)O" },
  { name: "Cyanide ion", smiles: "[C-]#N" },
];

function interpretSolubility(logS) {
  if (logS > -2) return { label: "High solubility", color: "success" };
  if (logS > -4) return { label: "Moderate solubility", color: "warning" };
  return { label: "Low solubility", color: "error" };
}

function riskColor(risk) {
  if (risk === "high") return "error";
  if (risk === "medium") return "warning";
  return "success";
}

function interpretLd50(score) {
  if (score >= 3) {
    return { label: "High acute-toxicity signal", color: "error" };
  }
  if (score >= 2) {
    return { label: "Moderate acute-toxicity signal", color: "warning" };
  }
  return { label: "Low predicted LD50-model signal", color: "success" };
}

function AdmetAnalysis() {
  const [smiles, setSmiles] = useState("");
  const [model, setModel] = useState("esol");
  const [results, setResults] = useState(null);
  const [error, setError] = useState(null);

  const predictAdmet = async () => {
    setError(null);
    setResults(null);

    const endpoint =
      model === "tox21"
        ? `${API_BASE}/api/ai/tox21/predict`
        : model === "ld50"
        ? `${API_BASE}/api/ai/ld50/predict`
        : model === "bbbp"
        ? `${API_BASE}/api/ai/bbbp/predict`
        : model === "clintox"
        ? `${API_BASE}/api/ai/clintox/predict`
        : `${API_BASE}/api/ai/gnn/predict`;

    const response = await fetch(endpoint, {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({ smiles }),
    });

    if (!response.ok) {
      setError("Prediction failed. Check the SMILES or backend.");
      return;
    }

    const data = await response.json();

    if (model === "esol") {
      const prediction = Number(data.prediction);
      const interpretation = interpretSolubility(prediction);

      setResults({
        type: "esol",
        model: data.model,
        smiles: data.smiles,
        property: "ESOL aqueous solubility",
        prediction: prediction.toFixed(3),
        units: "logS mol/L",
        interpretation,
      });
    } else if (model === "tox21") {
      setResults({
        type: "tox21",
        model: data.model,
        smiles: data.smiles,
        property: "Tox21 multi-label toxicity",
        overallRisk: data.overall_risk,
        maxProbability: Number(data.max_probability).toFixed(3),
        toxicityProbabilities: data.toxicity_probabilities,
      });
    } else if (model === "ld50") {
      const prediction = Number(data.prediction);
      const interpretation = interpretLd50(prediction);

      setResults({
        type: "ld50",
        model: data.model,
        smiles: data.smiles,
        property: "LD50 acute toxicity",
        prediction: prediction.toFixed(3),
        units: data.unit || "log(1/(mol/kg))",
        interpretation,
      });
    } else if (model === "bbbp") {
      setResults({
        type: "bbbp",
        model: data.model,
        smiles: data.smiles,
        property: "Blood-brain barrier permeability",
        probability: Number(data.probability).toFixed(3),
        label: data.label,
      });
    } else if (model === "clintox") {
      setResults({
        type: "clintox",
        model: data.model,
        smiles: data.smiles,
        property: "ClinTox clinical toxicity",
        clinicalToxicitySignal: data.clinical_toxicity_signal,
        ctToxProbability: Number(data.ct_tox_probability).toFixed(3),
        probabilities: data.probabilities,
      });
    }
  };

  return (
    <Paper sx={{ p: 3, background: "#1e1e1e" }}>
      <Typography variant="h5" gutterBottom>
        ⚗️ AI ADMET Prediction
      </Typography>

      <Typography variant="body1" sx={{ mb: 2 }}>
        Predict ADMET-related molecular properties from SMILES using AI models.
      </Typography>

      <Stack spacing={2} sx={{ mb: 3 }}>
        <TextField
          select
          label="Prediction model"
          value={model}
          onChange={(e) => setModel(e.target.value)}
          sx={{ width: "360px" }}
        >
          <MenuItem value="esol">ESOL solubility GNN</MenuItem>
          <MenuItem value="tox21">Tox21 toxicity GNN</MenuItem>
          <MenuItem value="ld50">LD50 acute toxicity GNN</MenuItem>
          <MenuItem value="bbbp">BBBP permeability GNN</MenuItem>
          <MenuItem value="clintox">ClinTox clinical toxicity GNN</MenuItem>
        </TextField>

        <Stack direction="row" spacing={2}>
          <TextField
            label="SMILES"
            variant="outlined"
            value={smiles}
            onChange={(e) => setSmiles(e.target.value)}
            sx={{ input: { color: "white" }, width: "420px" }}
          />

          <Button variant="contained" onClick={predictAdmet}>
            Predict
          </Button>
        </Stack>

        <Stack direction="row" spacing={1}>
          {exampleMolecules.map((mol) => (
            <Button
              key={mol.name}
              size="small"
              variant="outlined"
              onClick={() => setSmiles(mol.smiles)}
            >
              {mol.name}
            </Button>
          ))}
        </Stack>
      </Stack>

      <Divider sx={{ mb: 2 }} />

      <Typography variant="subtitle1" gutterBottom>
        Model status
      </Typography>

      <Stack direction="row" spacing={1} sx={{ mb: 2 }}>
        <Chip label="ESOL GNN" color={model === "esol" ? "success" : "default"} />
        <Chip label="Tox21 GNN" color={model === "tox21" ? "success" : "default"} />
        <Chip label="LD50 GNN" color={model === "ld50" ? "success" : "default"} />
        <Chip label="BBBP GNN" color={model === "bbbp" ? "success" : "default"} />
        <Chip
          label="ClinTox GNN"
          color={model === "clintox" ? "success" : "default"}
        />
      </Stack>

      {error && (
        <Alert severity="error" sx={{ mb: 2 }}>
          {error}
        </Alert>
      )}

      {results && results.type === "esol" && (
        <>
          <Typography variant="subtitle1" gutterBottom>
            Prediction result
          </Typography>

          <Chip
            label={results.interpretation.label}
            color={results.interpretation.color}
            sx={{ mb: 2 }}
          />

          <Table>
            <TableBody>
              <TableRow>
                <TableCell sx={{ color: "white", fontWeight: "bold" }}>
                  Model
                </TableCell>
                <TableCell sx={{ color: "white" }}>{results.model}</TableCell>
              </TableRow>

              <TableRow>
                <TableCell sx={{ color: "white", fontWeight: "bold" }}>
                  Input SMILES
                </TableCell>
                <TableCell sx={{ color: "white" }}>{results.smiles}</TableCell>
              </TableRow>

              <TableRow>
                <TableCell sx={{ color: "white", fontWeight: "bold" }}>
                  Property
                </TableCell>
                <TableCell sx={{ color: "white" }}>
                  {results.property}
                </TableCell>
              </TableRow>

              <TableRow>
                <TableCell sx={{ color: "white", fontWeight: "bold" }}>
                  Prediction
                </TableCell>
                <TableCell sx={{ color: "white" }}>
                  {results.prediction} {results.units}
                </TableCell>
              </TableRow>
            </TableBody>
          </Table>
        </>
      )}

      {results && results.type === "tox21" && (
        <>
          <Typography variant="subtitle1" gutterBottom>
            Toxicity prediction result
          </Typography>

          <Alert severity="warning" sx={{ mb: 2 }}>
            This model predicts only specific toxicity pathways (Tox21 assays).
            It does NOT represent overall toxicity or acute poisoning risk.
          </Alert>

          <Stack direction="row" spacing={1} sx={{ mb: 2 }}>
            <Chip
              label={`Tox21 signal: ${results.overallRisk}`}
              color={riskColor(results.overallRisk)}
            />
            <Chip label={`Max probability: ${results.maxProbability}`} />
          </Stack>

          <Table>
            <TableBody>
              {Object.entries(results.toxicityProbabilities).map(
                ([target, probability]) => (
                  <TableRow key={target}>
                    <TableCell sx={{ color: "white", fontWeight: "bold" }}>
                      {target}
                    </TableCell>
                    <TableCell sx={{ color: "white" }}>
                      {Number(probability).toFixed(3)}
                    </TableCell>
                  </TableRow>
                )
              )}
            </TableBody>
          </Table>
        </>
      )}

      {results && results.type === "ld50" && (
        <>
          <Typography variant="subtitle1" gutterBottom>
            Acute toxicity prediction result
          </Typography>

          <Alert severity="warning" sx={{ mb: 2 }}>
            This LD50 model is experimental and should not be used for safety,
            handling, medical, or regulatory decisions.
          </Alert>

          <Chip
            label={results.interpretation.label}
            color={results.interpretation.color}
            sx={{ mb: 2 }}
          />

          <Table>
            <TableBody>
              <TableRow>
                <TableCell sx={{ color: "white", fontWeight: "bold" }}>
                  Model
                </TableCell>
                <TableCell sx={{ color: "white" }}>{results.model}</TableCell>
              </TableRow>

              <TableRow>
                <TableCell sx={{ color: "white", fontWeight: "bold" }}>
                  Input SMILES
                </TableCell>
                <TableCell sx={{ color: "white" }}>{results.smiles}</TableCell>
              </TableRow>

              <TableRow>
                <TableCell sx={{ color: "white", fontWeight: "bold" }}>
                  Property
                </TableCell>
                <TableCell sx={{ color: "white" }}>
                  {results.property}
                </TableCell>
              </TableRow>

              <TableRow>
                <TableCell sx={{ color: "white", fontWeight: "bold" }}>
                  Prediction
                </TableCell>
                <TableCell sx={{ color: "white" }}>
                  {results.prediction} {results.units}
                </TableCell>
              </TableRow>
            </TableBody>
          </Table>
        </>
      )}

      {results && results.type === "bbbp" && (
        <>
          <Typography variant="subtitle1" gutterBottom>
            BBBP permeability prediction result
          </Typography>

          <Alert severity="warning" sx={{ mb: 2 }}>
            This model predicts a BBBP dataset signal, not guaranteed real
            CNS exposure or pharmacokinetic behaviour.
          </Alert>

          <Chip
            label={results.label}
            color={results.probability >= 0.5 ? "success" : "warning"}
            sx={{ mb: 2 }}
          />

          <Table>
            <TableBody>
              <TableRow>
                <TableCell sx={{ color: "white", fontWeight: "bold" }}>
                  Model
                </TableCell>
                <TableCell sx={{ color: "white" }}>{results.model}</TableCell>
              </TableRow>

              <TableRow>
                <TableCell sx={{ color: "white", fontWeight: "bold" }}>
                  Input SMILES
                </TableCell>
                <TableCell sx={{ color: "white" }}>{results.smiles}</TableCell>
              </TableRow>

              <TableRow>
                <TableCell sx={{ color: "white", fontWeight: "bold" }}>
                  Property
                </TableCell>
                <TableCell sx={{ color: "white" }}>
                  {results.property}
                </TableCell>
              </TableRow>

              <TableRow>
                <TableCell sx={{ color: "white", fontWeight: "bold" }}>
                  Probability
                </TableCell>
                <TableCell sx={{ color: "white" }}>
                  {results.probability}
                </TableCell>
              </TableRow>
            </TableBody>
          </Table>
        </>
      )}

      {results && results.type === "clintox" && (
        <>
          <Typography variant="subtitle1" gutterBottom>
            ClinTox prediction result
          </Typography>

          <Alert severity="warning" sx={{ mb: 2 }}>
            This model predicts ClinTox dataset signals. It is not clinical,
            regulatory, or medical evidence.
          </Alert>

          <Stack direction="row" spacing={1} sx={{ mb: 2 }}>
            <Chip
              label={`Clinical toxicity signal: ${results.clinicalToxicitySignal}`}
              color={riskColor(results.clinicalToxicitySignal)}
            />
            <Chip label={`CT_TOX probability: ${results.ctToxProbability}`} />
          </Stack>

          <Table>
            <TableBody>
              <TableRow>
                <TableCell sx={{ color: "white", fontWeight: "bold" }}>
                  Model
                </TableCell>
                <TableCell sx={{ color: "white" }}>{results.model}</TableCell>
              </TableRow>

              <TableRow>
                <TableCell sx={{ color: "white", fontWeight: "bold" }}>
                  Input SMILES
                </TableCell>
                <TableCell sx={{ color: "white" }}>{results.smiles}</TableCell>
              </TableRow>

              <TableRow>
                <TableCell sx={{ color: "white", fontWeight: "bold" }}>
                  Property
                </TableCell>
                <TableCell sx={{ color: "white" }}>
                  {results.property}
                </TableCell>
              </TableRow>

              {Object.entries(results.probabilities).map(
                ([target, probability]) => (
                  <TableRow key={target}>
                    <TableCell sx={{ color: "white", fontWeight: "bold" }}>
                      {target}
                    </TableCell>
                    <TableCell sx={{ color: "white" }}>
                      {Number(probability).toFixed(3)}
                    </TableCell>
                  </TableRow>
                )
              )}
            </TableBody>
          </Table>
        </>
      )}
    </Paper>
  );
}

export default AdmetAnalysis;