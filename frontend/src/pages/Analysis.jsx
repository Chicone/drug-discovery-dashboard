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
  ToggleButton,
  ToggleButtonGroup
} from "@mui/material";
import AnalysisPlot from "../components/analysis/AnalysisPlot";

export default function Analysis() {
  const [jobs, setJobs] = useState([]);
  const [selectedJobs, setSelectedJobs] = useState([]);
  const [plotData, setPlotData] = useState(null);
  const [loading, setLoading] = useState(false);
  const [whiteBG, setWhiteBG] = useState(false);

  const [metric, setMetric] = useState("rmsd");


console.log("plotData NOW:", plotData);


  // Fetch available MD jobs
  useEffect(() => {
    fetch("/api/md/jobs")
      .then((r) => r.json())
      .then((jobs) => setJobs(jobs))
      .catch((e) => console.error("Failed to load MD jobs", e));
  }, []);

const runAnalysis = async () => {
  if (selectedJobs.length === 0) return;

  setLoading(true);
  setPlotData(null);

  const params = selectedJobs.map((j) => `job_ids=${j}`).join("&");

  const endpoint =
    metric === "rmsd"
      ? "/api/analysis/ligand_rmsd"
      : "/api/analysis/ligand_com_distance";

  try {
    const r = await fetch(`${endpoint}?${params}`);
    const json = await r.json();

    console.log("ANALYSIS RESULT:", json);
    setPlotData(json);
  } catch (err) {
    console.error(err);
  }

  setLoading(false);
};


  const [plotMode, setPlotMode] = useState("individual");
  // "aggregate" | "individual"


  return (
    <Paper sx={{ p: 3, background: "#1e1e1e", color: "white" }}>
        <Typography variant="h5" gutterBottom>
          📈 {metric === "rmsd" ? "Ligand RMSD Analysis" : "Ligand–Pocket COM Analysis"}
        </Typography>
           <Box
          sx={{
            mt: 3,
            mb: 4,
            p: 3,
            borderRadius: 2,
            background: "#252525",
          }}
        >

          {/* Row 1: Job select + button */}
          <Stack direction="row" spacing={2} alignItems="center" sx={{ mb: 3 }}>
            <FormControl sx={{ minWidth: 280 }}>
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
              onClick={runAnalysis}
              disabled={selectedJobs.length === 0}
              sx={{ height: 56 }}
            >
              {metric === "rmsd" ? "Compute RMSD" : "Compute COM"}
            </Button>
          </Stack>

          {/* Row 2: Metric + Plot mode */}
          <Stack direction="row" spacing={3} alignItems="center" sx={{ mb: 2 }}>

            <ToggleButtonGroup
              value={metric}
              exclusive
              onChange={(e, newValue) => {
                if (newValue !== null) {
                  setMetric(newValue);
                  setPlotData(null);   // 🔥 CRITICAL
                }
              }}
              size="small"
            >
              <ToggleButton value="rmsd">RMSD</ToggleButton>
              <ToggleButton value="com">COM Distance</ToggleButton>
            </ToggleButtonGroup>

            <FormControl sx={{ minWidth: 200 }}>
              <InputLabel sx={{ color: "#ccc" }}>Plot Mode</InputLabel>
              <Select
                value={plotMode}
                label="Plot Mode"
                onChange={(e) => setPlotMode(e.target.value)}
                sx={{ color: "white" }}
                size="small"
              >
                <MenuItem value="individual">Individual</MenuItem>
                <MenuItem value="aggregate">Aggregate</MenuItem>
              </Select>
            </FormControl>

            <Box sx={{ display: "flex", alignItems: "center" }}>
              <label style={{ color: "#ccc", cursor: "pointer" }}>
                <input
                  type="checkbox"
                  checked={whiteBG}
                  onChange={() => setWhiteBG(!whiteBG)}
                  style={{ marginRight: 8 }}
                />
                White export
              </label>
            </Box>

          </Stack>

        </Box>

    <AnalysisPlot
      title={
        metric === "rmsd"
      ? "Ligand RMSD"
      : "Ligand–Pocket COM Distance"
    }

        plotData={plotData}
        plotMode={plotMode}
        whiteBG={whiteBG}

      />

    </Paper>
  );
}
