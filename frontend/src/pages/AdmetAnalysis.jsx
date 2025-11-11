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
} from "@mui/material";

function AdmetAnalysis() {
  const [smiles, setSmiles] = useState("");
  const [results, setResults] = useState(null);

  const predictAdmet = async () => {
    const mock = {
      "Drug-likeness (Lipinski)": "Yes",
      "Predicted logP": 2.8,
      "Predicted solubility": "Moderate",
      "Toxicity risk": "Low",
    };
    setResults(mock);
  };

  return (
    <Paper sx={{ p: 3, background: "#1e1e1e" }}>
      <Typography variant="h5" gutterBottom>
        ⚗️ ADMET Prediction
      </Typography>
      <Typography variant="body1" sx={{ mb: 2 }}>
        Evaluate absorption, distribution, metabolism, excretion, and toxicity properties.
      </Typography>

      <Stack direction="row" spacing={2} sx={{ mb: 2 }}>
        <TextField
          label="SMILES"
          variant="outlined"
          value={smiles}
          onChange={(e) => setSmiles(e.target.value)}
          sx={{ input: { color: "white" }, width: "300px" }}
        />
        <Button variant="contained" onClick={predictAdmet}>
          Predict
        </Button>
      </Stack>

      {results && (
        <Table>
          <TableBody>
            {Object.entries(results).map(([key, value]) => (
              <TableRow key={key}>
                <TableCell sx={{ color: "white" }}>{key}</TableCell>
                <TableCell sx={{ color: "white" }}>{value}</TableCell>
              </TableRow>
            ))}
          </TableBody>
        </Table>
      )}
    </Paper>
  );
}

export default AdmetAnalysis;
