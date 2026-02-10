// frontend/src/pages/Analysis.jsx
import { useState, useEffect } from "react";
import {
  Box,
  Paper,
  Typography,
  Button,
  Stack,
  MenuItem,
  Select,
  FormControl,
  InputLabel,
} from "@mui/material";
import AnalysisPlot from "../components/analysis/AnalysisPlot";

export default function Analysis() {
  const [jobs, setJobs] = useState([]);
  const [selectedJobs, setSelectedJobs] = useState([]);
  const [plotData, setPlotData] = useState(null);
  const [loading, setLoading] = useState(false);

console.log("plotData NOW:", plotData);


  // Fetch available MD jobs
  useEffect(() => {
    fetch("/api/md/jobs")
      .then((r) => r.json())
      .then((jobs) => setJobs(jobs))
      .catch((e) => console.error("Failed to load MD jobs", e));
  }, []);

const runRMSD = async () => {
  if (selectedJobs.length === 0) return;

  setLoading(true);
  setPlotData(null);

  const params = selectedJobs.map((j) => `job_ids=${j}`).join("&");

  try {
    const r = await fetch(`/api/analysis/ligand_rmsd?${params}`);
    const json = await r.json();

    console.log("RMSD RESULT:", json);

    setPlotData(json);  // 👈 IMPORTANT
  } catch (err) {
    console.error(err);
  }

  setLoading(false);
};

const [plotMode, setPlotMode] = useState("aggregate");
// "aggregate" | "individual"



  return (
    <Paper sx={{ p: 3, background: "#1e1e1e", color: "white" }}>
      <Typography variant="h5" gutterBottom>
        📈 Ligand RMSD Analysis
      </Typography>

      {/* Job selection */}
      <Stack direction="row" spacing={2} sx={{ mt: 2, mb: 3 }}>
        <FormControl sx={{ minWidth: 240 }}>
          <InputLabel sx={{ color: "#ccc" }}>Select MD Jobs</InputLabel>
          <Select
            multiple
            value={selectedJobs}
            label="Select MD Jobs"
            onChange={(e) => setSelectedJobs(e.target.value)}
            sx={{ color: "white" }}
          >
            {jobs.map((j) => (
              <MenuItem key={j.job_id} value={j.job_id}>
                {j.job_id}
              </MenuItem>
            ))}
          </Select>
        </FormControl>

        <Button
          variant="contained"
          color="secondary"
          onClick={runRMSD}
          disabled={selectedJobs.length === 0}
        >
          Compute RMSD
        </Button>
      </Stack>

      <FormControl sx={{ minWidth: 200 }}>
        <InputLabel sx={{ color: "#ccc" }}>Plot Mode</InputLabel>
        <Select
          value={plotMode}
          label="Plot Mode"
          onChange={(e) => setPlotMode(e.target.value)}
          sx={{ color: "white" }}
        >
          <MenuItem value="aggregate">Aggregate (mean/min/max)</MenuItem>
          <MenuItem value="individual">Individual traces</MenuItem>
        </Select>
      </FormControl>

      <AnalysisPlot
        title="Ligand RMSD"
        plotData={plotData}
        plotMode={plotMode}
      />

    </Paper>
  );
}
