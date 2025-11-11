import { useState, useEffect } from "react";
import {
  Box,
  Paper,
  Typography,
  TextField,
  Button,
  Stack,
  Divider,
  LinearProgress,
} from "@mui/material";
import MoleculeViewer from "../MoleculeViewer";

function MolecularDesign() {
  const [smiles, setSmiles] = useState("");
  const [data, setData] = useState(null);
  const [error, setError] = useState(null);
  const [RDKit, setRDKit] = useState(null);
  const [molName, setMolName] = useState("");

  useEffect(() => {
    const loadRDKit = async () => {
      const RDKitModule = await window.initRDKitModule();
      setRDKit(RDKitModule);
    };
    if (!window.initRDKitModule) {
      const script = document.createElement("script");
      script.src = "https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js";
      script.onload = loadRDKit;
      document.body.appendChild(script);
    } else {
      loadRDKit();
    }
  }, []);

  const fetchProperties = async () => {
    setError(null);
    setData(null);
    try {
      const res = await fetch(`/api/properties?smiles=${encodeURIComponent(smiles)}`);
      const json = await res.json();
      if (json.error) setError(json.error);
      else setData(json);
    } catch (e) {
      setError("Could not connect to backend");
    }
  };

  const renderMolecule = () => {
    if (!RDKit || !smiles) return null;
    try {
      const mol = RDKit.get_mol(smiles);
      const svg = mol.get_svg();
      mol.delete();
      return (
        <Box
          dangerouslySetInnerHTML={{ __html: svg }}
          sx={{
            border: "1px solid rgba(255,255,255,0.1)",
            borderRadius: 1,
            mt: 2,
            p: 1,
            backgroundColor: "#1a1a1a",
          }}
        />
      );
    } catch {
      return (
        <Typography variant="body2" color="error" sx={{ mt: 1 }}>
          Invalid SMILES for visualization
        </Typography>
      );
    }
  };

  return (
    <Paper
      sx={{
        p: 3,
        backgroundColor: "#1e1e1e",
        color: "white",
        maxWidth: "100%",
        overflowX: "hidden",
      }}
    >
      <Typography variant="h5" gutterBottom>
        🧪 Molecular Design
      </Typography>

      <Typography variant="body2" sx={{ mb: 2, color: "#aaa" }}>
        Enter a SMILES string or upload a .smi / .txt file to analyze molecular
        properties and visualize the structure.
      </Typography>

      <Stack direction={{ xs: "column", sm: "row" }} spacing={2} alignItems="center">
        <TextField
          variant="outlined"
          label="SMILES"
          placeholder="e.g. CCO"
          size="small"
          value={smiles}
          onChange={(e) => setSmiles(e.target.value)}
          sx={{
            width: { xs: "100%", sm: 300 },
            input: { color: "white" },
            label: { color: "#aaa" },
            "& .MuiOutlinedInput-root": {
              "& fieldset": { borderColor: "#555" },
              "&:hover fieldset": { borderColor: "#888" },
              "&.Mui-focused fieldset": { borderColor: "white" },
            },
          }}
        />

        <Button
          variant="contained"
          component="label"
          color="secondary"
          sx={{ whiteSpace: "nowrap" }}
        >
          Upload File
          <input
            type="file"
            accept=".smi,.txt"
            hidden
            onChange={(e) => {
              const file = e.target.files[0];
if (!file) return;
const reader = new FileReader();
reader.onload = (ev) => {
  const text = ev.target.result;
  const firstLine = text.split("\n").find((l) => l.trim().length > 0);
  if (firstLine) {
    const smilesFromFile = firstLine.split(/\s+/)[0];
    setSmiles(smilesFromFile);
    setMolName(file.name.replace(/\.[^/.]+$/, "")); // 👈 use filename as molecule name
  }
};
reader.readAsText(file);


            }}
          />
        </Button>

        <Button
          variant="contained"
          color="primary"
          onClick={fetchProperties}
          sx={{ whiteSpace: "nowrap" }}
        >
          Analyze
        </Button>
      </Stack>

      {error && (
        <Typography color="error" sx={{ mt: 2 }}>
          {error}
        </Typography>
      )}

      {data && (
        <Box sx={{ mt: 4 }}>
          <Divider sx={{ mb: 2, borderColor: "rgba(255,255,255,0.1)" }} />
          <Typography variant="h6" gutterBottom>
            Molecular Properties
          </Typography>

          {/* Unified rectangular layout */}
          <Box
            sx={{
              mt: 3,
              display: "flex",
              flexDirection: { xs: "column", md: "row" },
              alignItems: "stretch",
              justifyContent: "center",
              border: "1px solid rgba(255,255,255,0.1)",
              borderRadius: 2,
              overflow: "hidden",
            }}
          >
            {/* LEFT: 2D + properties + QED bar */}
            <Box
              sx={{
                flex: 1,
                display: "flex",
                flexDirection: "column",
                justifyContent: "flex-start",
                alignItems: "center",
                p: 2,
                borderRight: {
                  xs: "none",
                  md: "1px solid rgba(255,255,255,0.08)",
                },
                boxSizing: "border-box",
                gap: 2,
              }}
            >
            {molName && (
          <Typography
            variant="h6" // slightly larger than subtitle1
            sx={{
              fontWeight: 700,
              color: "primary.main",
              mb: 1,
              letterSpacing: 0.8,
              textTransform: "uppercase", // 👈 converts to all caps
            }}
          >
            {molName}
          </Typography>
        )}

              {renderMolecule()}

              {/* Properties table (without QED) */}
              <Box
                component="table"
                sx={{
                  fontSize: "0.9rem",
                  borderCollapse: "collapse",
                  width: "100%",
                  maxWidth: 320,
                  "& td, & th": {
                    border: "1px solid rgba(255,255,255,0.1)",
                    padding: "4px 6px",
                  },
                  "& th": { color: "#ccc" },
                }}
              >
                <tbody>
                  {Object.entries(data)
                    .filter(([key]) => !["insights", "QED"].includes(key))
                    .map(([key, value]) => (
                      <tr key={key}>
                        <td>
                          <strong>{key}</strong>
                        </td>
                        <td>{value}</td>
                      </tr>
                    ))}
                </tbody>
              </Box>

              {/* QED score bar */}
              {data?.QED !== undefined && (
                <Box sx={{ width: "100%", maxWidth: 320, mt: 1 }}>
                  <Typography
                    variant="subtitle2"
                    sx={{ color: "primary.main", fontWeight: 600, mb: 0.5 }}
                  >
                    QED Score
                  </Typography>
                  <LinearProgress
                    variant="determinate"
                    value={data.QED * 100}
                    sx={{
                      height: 8,
                      borderRadius: 2,
                      "& .MuiLinearProgress-bar": {
                        backgroundColor:
                          data.QED > 0.7
                            ? "#4caf50"
                            : data.QED > 0.4
                            ? "#ffb300"
                            : "#f44336",
                      },
                    }}
                  />
                  <Typography
                    variant="body2"
                    sx={{ color: "text.secondary", mt: 0.5, fontSize: "0.9rem" }}
                  >
                    {data.QED.toFixed(2)}{" "}
                    {data.QED > 0.7
                      ? "– Highly drug-like"
                      : data.QED > 0.4
                      ? "– Moderately drug-like"
                      : "– Low drug-likeness"}
                  </Typography>
                </Box>
              )}
            </Box>

            {/* RIGHT: 3D Viewer */}
            <Box
              sx={{
                flex: 1.2,
                display: "flex",
                justifyContent: "center",
                alignItems: "center",
                p: 2,
                boxSizing: "border-box",
              }}
            >
              <MoleculeViewer smiles={smiles} />
            </Box>
          </Box>

          {/* Design Insights BELOW the rectangle */}
          {data?.insights && (
            <Box
              sx={{
                mt: 3,
                background: "rgba(255,255,255,0.05)",
                borderRadius: 1,
                border: "1px solid rgba(255,255,255,0.1)",
                p: 2,
              }}
            >
              <Typography
                variant="h6"
                sx={{ fontWeight: 600, mb: 1, color: "primary.main" }}
              >
                Design Insights
              </Typography>
              {data.insights.map((line, i) => (
                <Typography
                  key={i}
                  variant="body2"
                  sx={{ color: "text.secondary", mb: 0.5, fontSize: "0.9rem" }}
                >
                  • {line}
                </Typography>
              ))}
            </Box>
          )}
        </Box>
      )}
    </Paper>
  );
}

export default MolecularDesign;
