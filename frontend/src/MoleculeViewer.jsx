import { useEffect, useRef, useState } from "react";
import {
  Box,
  FormControl,
  InputLabel,
  MenuItem,
  Select,
} from "@mui/material";

function MoleculeViewer({ smiles, pdbFile }) {
  const viewerRef = useRef(null);
  const [style, setStyle] = useState("stick"); // for SMILES

  useEffect(() => {
    if (!window.$3Dmol) return;

    const element = viewerRef.current;
    element.innerHTML = "";
    const viewer = window.$3Dmol.createViewer(element, {
      backgroundColor: "#111111",
    });

    // ---- Case 1: PDB file uploaded ----
    if (pdbFile) {
      const reader = new FileReader();
      reader.onload = (e) => {
        const pdbData = e.target.result;
        viewer.addModel(pdbData, "pdb");
        viewer.setStyle({}, { cartoon: { color: "spectrum" } });
        viewer.zoomTo();
        viewer.render();
      };
      reader.readAsText(pdbFile);
      return;
    }

    // ---- Case 2: SMILES (default) ----
    if (!smiles) return;

    fetch(`/api/mol3d?smiles=${encodeURIComponent(smiles)}`)
      .then((r) => r.text())
      .then((molData) => {
        if (!molData) return;
        viewer.addModel(molData, "mol");

        if (style === "ballstick") {
          viewer.setStyle({}, {
            stick: { radius: 0.15 },
            sphere: { scale: 0.25 },
          });
        } else if (style === "sphere") {
          viewer.setStyle({}, { sphere: { scale: 0.3 } });
        } else {
          viewer.setStyle({}, { [style]: {} });
        }

        viewer.zoomTo();
        viewer.render();
      });
  }, [smiles, pdbFile, style]);

  // hide style dropdown when showing a PDB file
  const showStyleSelector = !pdbFile;

  return (
    <Box
      sx={{
        display: "flex",
        flexDirection: "column",
        alignItems: "flex-start",
        width: "100%",
      }}
    >
      {showStyleSelector && (
        <FormControl
          variant="outlined"
          size="small"
          sx={{
            mb: 1,
            minWidth: 140,
            "& .MuiInputLabel-root": { color: "#aaa" },
            "& .MuiOutlinedInput-root": {
              "& fieldset": { borderColor: "#555" },
              "&:hover fieldset": { borderColor: "#888" },
              "&.Mui-focused fieldset": { borderColor: "white" },
            },
          }}
        >
          <InputLabel>Style</InputLabel>
          <Select
            value={style}
            onChange={(e) => setStyle(e.target.value)}
            label="Style"
            sx={{ color: "white" }}
          >
            <MenuItem value="stick">Sticks</MenuItem>
            <MenuItem value="sphere">Balls (Space-fill)</MenuItem>
            <MenuItem value="ballstick">Ball-and-Stick</MenuItem>
            <MenuItem value="line">Wireframe</MenuItem>
            <MenuItem value="cross">Cross</MenuItem>
            <MenuItem value="cartoon">Cartoon</MenuItem>
          </Select>
        </FormControl>
      )}

      <div
        ref={viewerRef}
        style={{
          width: "600px",
          height: "600px",
          border: "1px solid #555",
          borderRadius: "6px",
          background: "#111",
          marginTop: "0.3rem",
          position: "relative",
          overflow: "hidden",
        }}
      />
    </Box>
  );
}

export default MoleculeViewer;
