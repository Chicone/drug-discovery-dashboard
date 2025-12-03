import { useEffect, useRef, useState } from "react";
import { Box, FormControl, InputLabel, MenuItem, Select, Typography } from "@mui/material";

/**
 * Props:
 *  - smiles?: string
 *  - pdbText?: string         // raw PDB text
 *  - pdbFile?: File           // File object from <input type=file>
 *  - pdbUrl?: string          // URL to fetch a PDB from
 */
export default function MoleculeViewer({ smiles, pdbText, pdbFile, pdbUrl }) {
  const viewerRef = useRef(null);
  const [style, setStyle] = useState("stick");
  const [err, setErr] = useState("");

  // helper to render PDB text (protein + ligand)
  const renderPdbString = (pdbStr) => {
    if (!window.$3Dmol || !viewerRef.current) return;
    setErr("");

    const aminoAcids = new Set([
      "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
      "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    ]);

    const el = viewerRef.current;
    el.innerHTML = "";
    const viewer = window.$3Dmol.createViewer(el, { backgroundColor: "#111" });

    let model;
    try {
      model = viewer.addModel(pdbStr, "pdb");
    } catch (e) {
      console.error("3Dmol addModel failed:", e);
      setErr("Could not parse PDB file.");
      return;
    }

    if (!model) {
      setErr("No model returned from PDB.");
      return;
    }

    const atoms = model.selectedAtoms();
    if (!atoms || atoms.length === 0) {
      setErr("No atoms found in the PDB.");
      return;
    }

    const hasProtein = atoms.some((a) => aminoAcids.has(a.resn));
    const hasLigand = atoms.some((a) => !aminoAcids.has(a.resn));

//     if (hasProtein && hasLigand) {
//       // Protein → cartoon
//       viewer.setStyle(
//         { resn: Array.from(aminoAcids) },
//         { cartoon: { color: "spectrum", opacity: 0.85 } }
//       );
//
//       // Ligand (e.g. UNL) → atoms only, no bonds
//       viewer.setStyle(
//         { not: { resn: Array.from(aminoAcids) } },
//         { sphere: { radius: 0.6 } }   // you can adjust radius (0.5–0.8)
//       );


    if (hasProtein && hasLigand) {
      viewer.setStyle(
        { resn: Array.from(aminoAcids) },
        { cartoon: { color: "spectrum" } }
      );
      viewer.setStyle(
        { not: { resn: Array.from(aminoAcids) } },
        { stick: { colorscheme: "greenCarbon", radius: 0.2 } }
      );
    } else if (hasProtein) {
      viewer.setStyle({}, { cartoon: { color: "spectrum" } });
    } else {
      viewer.setStyle({}, { stick: { colorscheme: "greenCarbon", radius: 0.2 } });
    }

    viewer.zoomTo();
    viewer.render();
    setTimeout(() => viewer.resize(), 0);
  };

  // main rendering effect
  useEffect(() => {
    if (!window.$3Dmol || !viewerRef.current) return;

    // Case 1: explicit PDB text (e.g., from docking)
    if (pdbText && pdbText.trim().length > 0) {
      renderPdbString(pdbText);
      return;
    }

    // Case 2: uploaded PDB file (legacy compatibility)
    if (pdbFile instanceof File) {
      const reader = new FileReader();
      reader.onload = (e) => renderPdbString(e.target.result);
      reader.onerror = () => setErr("Failed to read uploaded PDB file.");
      reader.readAsText(pdbFile);
      return;
    }

    // Case 3: SMILES mode
    if (smiles) {
      const el = viewerRef.current;
      el.innerHTML = "";
      const viewer = window.$3Dmol.createViewer(el, { backgroundColor: "#111" });

      fetch(`/api/mol3d?smiles=${encodeURIComponent(smiles)}`)
        .then((r) => r.text())
        .then((molData) => {
          if (!molData) return;
          viewer.addModel(molData, "mol");
          viewer.setStyle({}, { stick: { colorscheme: "cyanCarbon" } });
          viewer.zoomTo();
          viewer.render();
          setTimeout(() => viewer.resize(), 0);
        })
        .catch(() => setErr("Failed to load 3D structure from SMILES."));
    }
  }, [pdbText, pdbFile, smiles]);

  const showStyleSelector = !(pdbText || pdbFile || pdbUrl);

  return (
    <Box sx={{ display: "flex", flexDirection: "column", alignItems: "flex-start", width: "100%" }}>
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
          <Select value={style} onChange={(e) => setStyle(e.target.value)} label="Style" sx={{ color: "white" }}>
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

      {err && (
        <Typography sx={{ mt: 1 }} color="error">
          {err}
        </Typography>
      )}
    </Box>
  );
}
