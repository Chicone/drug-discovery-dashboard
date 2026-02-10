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
import Plot from "react-plotly.js";

export default function Analysis() {
  const [jobs, setJobs] = useState([]);
  const [selectedJobs, setSelectedJobs] = useState([]);
  const [analysisData, setAnalysisData] = useState(null);
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

const dt_ps = Number(
  plotData?.replicas?.find((r) => r?.timestep_ps != null)?.timestep_ps
);

const x_ns =
  plotData?.aggregate?.mean?.map((_, i) =>
    Number.isFinite(dt_ps) ? (i * dt_ps) / 1000.0 : i
  ) ?? [];


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

      {/* Plot */}
     {plotData && plotData.aggregate && (
          <Box sx={{ mt: 4 }}>
            <Plot
              data={[
                { x: x_ns, y: plotData.aggregate.mean, type: "scatter", mode: "lines", name: "Mean RMSD" },
                { x: x_ns, y: plotData.aggregate.min,  type: "scatter", mode: "lines", name: "Min" },
                { x: x_ns, y: plotData.aggregate.max,  type: "scatter", mode: "lines", name: "Max" },
              ]}
              layout={{
                title: "Ligand RMSD over Time",
                paper_bgcolor: "#121212",
                plot_bgcolor: "#121212",
                font: { color: "white" },
                margin: { t: 60, r: 40, b: 80, l: 60 },
                xaxis: { title: { text: "Time (ns)" }, showline: true, linecolor: "white" },
                yaxis: { title: { text: "RMSD (Å)" }, showline: true, linecolor: "white" },
              }}
              style={{ width: "100%", height: "500px" }}
              config={{ displayModeBar: false }}
            />
          </Box>
        )}

    </Paper>
  );
}
