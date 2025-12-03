import { Paper, Typography, Button, Stack } from "@mui/material";

function Reports() {
  return (
    <Paper sx={{ p: 3, background: "#1e1e1e" }}>
      <Typography variant="h5" gutterBottom>
        📊 Reports
      </Typography>
      <Typography variant="body1" sx={{ mb: 2 }}>
        Generate, export, and review analysis summaries for your virtual screening and
        computational experiments.
      </Typography>

      <Stack direction="row" spacing={2}>
        <Button variant="contained" color="secondary" disabled>
          Export to PDF (coming soon)
        </Button>
        <Button variant="outlined" disabled>
          Export to CSV
        </Button>
      </Stack>

      <Typography variant="body2" sx={{ mt: 3, color: "#aaa" }}>
        🧠 Planned features:
        <ul>
          <li>Integrate property and docking results into formatted reports</li>
          <li>Support PDF and Excel exports</li>
          <li>Generate visual plots (RMSD, activity, etc.)</li>
        </ul>
      </Typography>
    </Paper>
  );
}

export default Reports;
