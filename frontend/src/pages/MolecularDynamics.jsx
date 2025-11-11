import { Paper, Typography, Button, Stack } from "@mui/material";

function MolecularDynamics() {
  return (
    <Paper sx={{ p: 3, background: "#1e1e1e" }}>
      <Typography variant="h5" gutterBottom>
        💫 Molecular Dynamics
      </Typography>
      <Typography variant="body1" sx={{ mb: 2 }}>
        This module will allow you to visualize and analyze molecular dynamics simulations.
        You’ll be able to explore atomistic motions, ligand stability, and conformational
        changes over time.
      </Typography>

      <Stack direction="row" spacing={2}>
        <Button variant="contained" color="secondary" disabled>
          Upload trajectory (coming soon)
        </Button>
        <Button variant="outlined" disabled>
          Run short simulation
        </Button>
      </Stack>

      <Typography variant="body2" sx={{ mt: 3, color: "#aaa" }}>
        🔬 Feature roadmap:
        <ul>
          <li>Integrate 3Dmol.js or NGL.js for trajectory visualization</li>
          <li>Connect to backend (OpenMM / GROMACS) for simulations</li>
          <li>Display RMSD, energy, and hydrogen bond plots</li>
        </ul>
      </Typography>
    </Paper>
  );
}

export default MolecularDynamics;
