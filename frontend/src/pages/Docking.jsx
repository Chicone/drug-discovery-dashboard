import { useState } from "react";
import {
  Box,
  Paper,
  Typography,
  Button,
  Stack,
  LinearProgress,
  Divider,
} from "@mui/material";

import MoleculeViewer from "../MoleculeViewer";
import Grid from "@mui/material/Grid";
import TextField from "@mui/material/TextField";

function Docking() {
  const [receptorFile, setReceptorFile] = useState(null);
  const [ligandFile, setLigandFile] = useState(null);
  const [pdbText, setPdbText] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState("");
  const [poses, setPoses] = useState([]);
  const [selectedPoseIdx, setSelectedPoseIdx] = useState(0);
  const [dockParams, setDockParams] = useState({
    num_modes: 10,
    exhaustiveness: 8,
    cpu: 0,

    center_x: 0,
    center_y: 0,
    center_z: 0,

    size_x: 20,
    size_y: 20,
    size_z: 20,
  });

  const updateDockParam = (key, value) => {
    setDockParams((p) => ({ ...p, [key]: value }));
  };


  const handleDocking = async () => {
    if (!receptorFile || !ligandFile) {
      setError("Please upload both receptor and ligand files.");
      return;
    }

    setLoading(true);
    setError("");
    setPdbText("");

    const formData = new FormData();
    formData.append("receptor", receptorFile);
    formData.append("ligand", ligandFile);
    formData.append("num_modes", dockParams.num_modes);
    formData.append("exhaustiveness", dockParams.exhaustiveness);
    formData.append("cpu", dockParams.cpu);
    formData.append("center_x", dockParams.center_x);
    formData.append("center_y", dockParams.center_y);
    formData.append("center_z", dockParams.center_z);
    formData.append("size_x", dockParams.size_x);
    formData.append("size_y", dockParams.size_y);
    formData.append("size_z", dockParams.size_z);

    try {
    const res = await fetch("/api/dock_vina", {
      method: "POST",
      body: formData,
    });

    if (!res.ok) {
      const errText = await res.text();
      throw new Error(errText || "Docking failed");
    }

    const data = await res.json();

    // data.poses is an array of { mode, score, pdb }
    if (!data.poses || data.poses.length === 0) {
      setError("Docking completed but returned no poses.");
      return;
    }

    setPoses(data.poses);
    setSelectedPoseIdx(0);
    setPdbText(data.poses[0].pdb);

    } catch (err) {
      console.error("Docking error:", err);
      setError(err.message);
    } finally {
      setLoading(false);
    }
  };

  return (
    <Paper
      sx={{
        p: 3,
        background: "#1e1e1e",
        color: "white",
        maxWidth: "100%",
        overflowX: "hidden",
      }}
    >
      <Typography variant="h5" gutterBottom>
        🧠 Molecular Docking
      </Typography>

      <Typography variant="body1" sx={{ mb: 2, color: "#aaa" }}>
        Predict how a small molecule binds to a target protein.
        Upload structures and analyze binding poses.
      </Typography>

<Grid container spacing={2} sx={{ mb: 2 }} alignItems="flex-start">
  {/* Left control panel */}
    <Grid item xs={12} md={4} lg={3}>
    <Paper
      elevation={3}
      sx={{
        p: 2,
        borderRadius: 3,
        background: "#2a2a2a",
        width: "100%",
        maxWidth: 300,
        mx: "auto",
      }}
    >
      <Typography variant="h6" sx={{ mb: 2 }}>
        Docking setup
      </Typography>

      {/* Uploads */}
      <Stack spacing={1} sx={{ mb: 1.5 }}>
        <Button
          variant="contained"
          component="label"
          color="secondary"
          fullWidth
          size="small"
          sx={{ py: 0.8 }}
        >
          Upload Receptor (.pdb)
          <input
            type="file"
            accept=".pdb"
            hidden
            onChange={(e) => setReceptorFile(e.target.files[0])}
          />
        </Button>

        <Button
          variant="contained"
          component="label"
          color="secondary"
          fullWidth
          size="small"
          sx={{ py: 0.8 }}
        >
          Upload Ligand (.smi or .pdb)
          <input
            type="file"
            accept=".smi,.pdb"
            hidden
            onChange={(e) => setLigandFile(e.target.files[0])}
          />
        </Button>
      </Stack>

      <Divider sx={{ my: 2 }} />

      {/* Parameters */}
<Typography variant="subtitle1" sx={{ mb: 0.5 }}>
  Parameters
</Typography>

<Box sx={{ maxWidth: 150, mx: "left" }}>
  <Stack spacing={1}>
    <TextField
      label="Number of modes"
      type="number"
      fullWidth
      size="small"
      margin="dense"
      value={dockParams.num_modes}
      onChange={(e) =>
        setDockParams((p) => ({
          ...p,
          num_modes: Number(e.target.value),
        }))
      }
    />

    <TextField
      label="Exhaustiveness"
      type="number"
      fullWidth
      size="small"
      margin="dense"
      value={dockParams.exhaustiveness}
      onChange={(e) =>
        setDockParams((p) => ({
          ...p,
          exhaustiveness: Number(e.target.value),
        }))
      }
    />

    <TextField
      label="CPU cores (0 = auto)"
      type="number"
      fullWidth
      size="small"
      margin="dense"
      value={dockParams.cpu}
      onChange={(e) =>
        setDockParams((p) => ({
          ...p,
          cpu: Number(e.target.value),
        }))
      }
    />
  </Stack>
</Box>

      <Divider sx={{ my: 1.5 }} />

<Typography variant="subtitle2"  sx={{ mb: 0.5 }}>
  Docking box (Å)
</Typography>

<Box
  sx={{
    width: "100%",
    maxWidth: "100%",
    display: "grid",
    gridTemplateColumns: "repeat(3, minmax(0, 1fr))",
    gap: 1,
  }}
>
  <TextField
    label="Center X"
    type="number"
    size="small"
    margin="dense"
    sx={{ minWidth: 0 }}
    value={dockParams.center_x}
    onChange={(e) =>
      setDockParams((p) => ({ ...p, center_x: Number(e.target.value) }))
    }
  />

  <TextField
    label="Center Y"
    type="number"
    size="small"
    margin="dense"
    sx={{ minWidth: 0 }}
    value={dockParams.center_y}
    onChange={(e) =>
      setDockParams((p) => ({ ...p, center_y: Number(e.target.value) }))
    }
  />

  <TextField
    label="Center Z"
    type="number"
    size="small"
    margin="dense"
    sx={{ minWidth: 0 }}
    value={dockParams.center_z}
    onChange={(e) =>
      setDockParams((p) => ({ ...p, center_z: Number(e.target.value) }))
    }
  />

  <TextField
    label="Size X"
    type="number"
    size="small"
    margin="dense"
    sx={{ minWidth: 0 }}
    inputProps={{ min: 1 }}
    value={dockParams.size_x}
    onChange={(e) =>
      setDockParams((p) => ({ ...p, size_x: Number(e.target.value) }))
    }
  />

  <TextField
    label="Size Y"
    type="number"
    size="small"
    margin="dense"
    sx={{ minWidth: 0 }}
    inputProps={{ min: 1 }}
    value={dockParams.size_y}
    onChange={(e) =>
      setDockParams((p) => ({ ...p, size_y: Number(e.target.value) }))
    }
  />

  <TextField
    label="Size Z"
    type="number"
    size="small"
    margin="dense"
    sx={{ minWidth: 0 }}
    inputProps={{ min: 1 }}
    value={dockParams.size_z}
    onChange={(e) =>
      setDockParams((p) => ({ ...p, size_z: Number(e.target.value) }))
    }
  />
</Box>

      <Divider sx={{ my: 1.5 }} />


      {/* Run */}
      <Button
        variant="contained"
        color="primary"
        fullWidth
        size="large"
        onClick={handleDocking}
        disabled={loading || !receptorFile || !ligandFile}
      >
        Run Docking
      </Button>
    </Paper>
  </Grid>

  {/* Right side (optional area for future viewer/logs) */}
  <Grid item xs={12} md={7} lg={8}>
    {/* Keep empty for now, or add logs/preview later */}
    <Grid item xs={12} md={8} lg={9}>
  {pdbText && (
    <Box>
      <Typography variant="h6" sx={{ mb: 1 }}>
        Docked Complex Visualization
      </Typography>
      <MoleculeViewer pdbText={pdbText} />
    </Box>
  )}

  {poses.length > 0 && (
    <Box sx={{ mt: 3 }}>
      <Typography variant="subtitle1" sx={{ mb: 1 }}>
        Docking poses
      </Typography>

      <Box sx={{ display: "flex", gap: 1, flexWrap: "wrap" }}>
        {poses.map((pose, i) => (
          <Button
            key={pose.mode ?? i}
            size="small"
            variant={i === selectedPoseIdx ? "contained" : "outlined"}
            onClick={() => {
              setSelectedPoseIdx(i);
              setPdbText(pose.pdb);
            }}
          >
            Pose {pose.mode}
            {typeof pose.score === "number" ? ` (${pose.score.toFixed(2)})` : ""}
          </Button>
        ))}
      </Box>
    </Box>
  )}
</Grid>

  </Grid>
</Grid>

      {loading && (
        <Box sx={{ mt: 3 }}>
          <Typography variant="body2" sx={{ mb: 1 }}>
            Running AutoDock Vina...
          </Typography>
          <LinearProgress />
        </Box>
      )}

      {error && (
        <Typography color="error" sx={{ mt: 2 }}>
          {error}
        </Typography>
      )}





      {!pdbText && !loading && (
        <Typography variant="body2" sx={{ color: "#aaa", mt: 3 }}>
          Upload a receptor and ligand to begin docking.
        </Typography>
      )}
    </Paper>
  );
}

export default Docking;
