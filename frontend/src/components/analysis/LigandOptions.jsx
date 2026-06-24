import { Box, Typography, Checkbox, FormControlLabel } from "@mui/material";

export default function LigandOptions() {
  return (
    <Box sx={{ mb: 3 }}>
      <Typography sx={{ fontWeight: 600, mb: 1 }}>Ligand Stability</Typography>

      <FormControlLabel
        control={<Checkbox />}
        label="COM Drift"
        sx={{ color: "white" }}
      />

      <FormControlLabel
        control={<Checkbox />}
        label="RMSD of ligand only"
        sx={{ color: "white" }}
      />

      <FormControlLabel
        control={<Checkbox />}
        label="Orientation metric"
        sx={{ color: "white" }}
      />
    </Box>
  );
}
