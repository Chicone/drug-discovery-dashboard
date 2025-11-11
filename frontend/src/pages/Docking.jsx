import { Paper, Typography, Button, Stack, TextField } from "@mui/material";

function Docking() {
  return (
    <Paper sx={{ p: 3, background: "#1e1e1e" }}>
      <Typography variant="h5" gutterBottom>
        🧠 Molecular Docking
      </Typography>
      <Typography variant="body1" sx={{ mb: 2 }}>
        Predict how a small molecule binds to a target protein. Upload structures and
        analyze binding poses.
      </Typography>

      <Stack direction="row" spacing={2} sx={{ mb: 2 }}>
        <TextField
          label="Ligand SMILES or file"
          variant="outlined"
          sx={{ input: { color: "white" }, width: "300px" }}
        />
        <Button variant="contained" disabled>
          Upload Protein (coming soon)
        </Button>
        <Button variant="contained" color="secondary" disabled>
          Run Docking
        </Button>
      </Stack>

      <Typography variant="body2" sx={{ color: "#aaa" }}>
        Planned features:
        <ul>
          <li>Integrate AutoDock Vina or RDKit docking backend</li>
          <li>Visualize binding poses and scores</li>
          <li>Estimate binding energy and interaction hotspots</li>
        </ul>
      </Typography>
    </Paper>
  );
}

export default Docking;
