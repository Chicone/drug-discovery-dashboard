import { useState } from "react";
import {
  Box,
  Paper,
  Typography,
  Button,
  Stack,
  LinearProgress,
} from "@mui/material";
import MoleculeViewer from "../MoleculeViewer";

function Docking() {
  const [receptorFile, setReceptorFile] = useState(null);
  const [ligandFile, setLigandFile] = useState(null);
  const [pdbText, setPdbText] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState("");

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

    try {
    const res = await fetch("/api/dock_vina", {
      method: "POST",
      body: formData,
    });

    if (!res.ok) {
      const errText = await res.text();
      throw new Error(errText || "Docking failed");
    }

    // Decode FileResponse properly as text
    const arrayBuffer = await res.arrayBuffer();
    const decoder = new TextDecoder("utf-8");
    const pdb = decoder.decode(arrayBuffer);

    // Confirm the PDB contains ATOM lines
    if (!pdb.includes("ATOM") && !pdb.includes("HETATM")) {
      console.warn("⚠️ Docking returned empty or invalid PDB:", pdb.slice(0, 200));
      setError("Docking completed but no valid PDB content received.");
    } else {
      setPdbText(pdb);
    }

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

      <Stack
        direction={{ xs: "column", sm: "row" }}
        spacing={2}
        alignItems="center"
        sx={{ mb: 2 }}
      >
        <Button variant="contained" component="label" color="secondary">
          Upload Receptor (.pdb)
          <input
            type="file"
            accept=".pdb"
            hidden
            onChange={(e) => setReceptorFile(e.target.files[0])}
          />
        </Button>

        <Button variant="contained" component="label" color="secondary">
          Upload Ligand (.smi or .pdb)
          <input
            type="file"
            accept=".smi,.pdb"
            hidden
            onChange={(e) => setLigandFile(e.target.files[0])}
          />
        </Button>

        <Button
          variant="contained"
          color="primary"
          onClick={handleDocking}
          disabled={loading}
        >
          Run Docking
        </Button>
      </Stack>

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

      {pdbText && (
        <Box sx={{ mt: 4 }}>
          <Typography variant="h6" sx={{ mb: 1 }}>
            Docked Complex Visualization
          </Typography>
          <MoleculeViewer pdbText={pdbText} />
        </Box>
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
