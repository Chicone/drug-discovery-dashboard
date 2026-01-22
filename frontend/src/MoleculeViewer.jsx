import { useEffect, useRef, useState } from "react";
import { Box, FormControl, InputLabel, MenuItem, Select, Typography } from "@mui/material";

const AMINO_ACIDS_ARR = [
  "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
  "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
];

const AMINO_ACIDS_SET = new Set(AMINO_ACIDS_ARR);


/**
 * Props:
 *  - smiles?: string
 *  - pdbText?: string         // raw PDB text
 *  - pdbFile?: File           // File object from <input type=file>
 *  - pdbUrl?: string          // URL to fetch a PDB from
 */
export default function MoleculeViewer({ smiles, pdbText, pdbFile, pdbUrl }) {
  const [manualHelixMap, setManualHelixMap] = useState(null); // object or null
  const [structureKey, setStructureKey] = useState("");       // e.g. "2YDV"
  const [useManualMap, setUseManualMap] = useState(false);
  const viewerRef = useRef(null);
  const viewerObjRef = useRef(null);
  const hasRenderedOnceRef = useRef(false);

  const [style, setStyle] = useState("stick");
  const [err, setErr] = useState("");

  // --- Helix selection state ---
  const [helices, setHelices] = useState([]);
  const [selectedHelixKeys, setSelectedHelixKeys] = useState(new Set());

  const helicesComputedRef = useRef(false);

  const onManualMapFile = async (file) => {
    if (!file) return;
    try {
      const text = await file.text();
      const obj = JSON.parse(text);
      setManualHelixMap(obj);
      setErr("");
    } catch (e) {
      console.error("Manual helix map parse error:", e);
      setErr("Failed to parse manual helix map JSON.");
    }
  };

  const buildHelicesFromManualMap = (key) => {
    if (!manualHelixMap || !key) return null;

    const entry = manualHelixMap[key.toUpperCase()];
    if (!entry) return null;

    const chain = entry.chain || "A";
    const helicesOut = [];

    // Expect keys like "TM1", "TM2", ... "TM7", optionally "H8"
    const ordered = Object.keys(entry)
      .filter((k) => k !== "chain")
      .sort((a, b) => {
        // TM1..TM7 then H8
        const ra = parseInt(a.replace(/\D/g, ""), 10);
        const rb = parseInt(b.replace(/\D/g, ""), 10);
        return (ra || 99) - (rb || 99) || a.localeCompare(b);
      });

    for (const name of ordered) {
      const range = entry[name];
      if (!Array.isArray(range) || range.length !== 2) continue;

      const start = Number(range[0]);
      const end = Number(range[1]);
      if (!Number.isFinite(start) || !Number.isFinite(end)) continue;

      const lo = Math.min(start, end);
      const hi = Math.max(start, end);

      helicesOut.push({
        key: name, // use "TM1" etc as stable key
        label: `${name} — ${chain}:${lo}-${hi} (${hi - lo + 1} res)`,
        chain,
        start: lo,
        end: hi,
      });
    }

    return helicesOut.length > 0 ? helicesOut : null;
  };


  const toggleHelix = (key) => {
    setSelectedHelixKeys((prev) => {
      const next = new Set(prev);
      next.has(key) ? next.delete(key) : next.add(key);
      return next;
    });
  };

  const clearHelices = () => setSelectedHelixKeys(new Set());


  useEffect(() => {
    if (!window.$3Dmol || !viewerRef.current) return;

    // Create once
    viewerObjRef.current = window.$3Dmol.createViewer(viewerRef.current, {
      backgroundColor: "#111",
    });

    // optional: resize on next tick
    setTimeout(() => viewerObjRef.current?.resize(), 0);

    return () => {
      viewerObjRef.current = null;
      hasRenderedOnceRef.current = false;
    };
  }, []);


  // helper to render PDB text (protein + ligand)
  const renderPdbString = (pdbStr) => {
    if (!window.$3Dmol || !viewerRef.current || !viewerObjRef.current) return;
    setErr("");

    const viewer = viewerObjRef.current;

    // Save current view (orientation/zoom/translation)
    const prevView = hasRenderedOnceRef.current ? viewer.getView() : null;

    // Update model without recreating viewer
    viewer.removeAllModels();

    let model;
    try {
      model = viewer.addModel(pdbStr, "pdb", { computeStruct: true });
    } catch (e) {
      console.error("3Dmol addModel failed:", e);
      setErr("Could not parse PDB file.");
      return;
    }

    if (!model) {
      setErr("No model returned from PDB.");
      return;
    }

    // Compute helices ONCE per PDB (manual-only)
    if (!helicesComputedRef.current) {
      const manual = useManualMap ? buildHelicesFromManualMap(structureKey) : null;

      setHelices(Array.isArray(manual) ? manual : []);
      helicesComputedRef.current = true;
    }



    const atoms = model.selectedAtoms();
    if (!atoms || atoms.length === 0) {
      setErr("No atoms found in the PDB.");
      return;
    }

    const hasProtein = atoms.some((a) => AMINO_ACIDS_SET.has(a.resn));
    const hasLigand = atoms.some((a) => !AMINO_ACIDS_SET.has(a.resn));


    // --- SAFE styling block (never blanks the viewer) ---

    const safeSetStyle = (sel, sty) => {
      try {
        viewer.setStyle(sel, sty);
      } catch (e) {
        console.error("viewer.setStyle failed", sel, sty, e);
        setErr(`Style error: ${e?.message || e}`);
      }
    };

    if (hasProtein && hasLigand) {
      // Protein as cartoon
      safeSetStyle(
        { resn: AMINO_ACIDS_ARR },
        { cartoon: { color: "spectrum" } }
      );

      // Ligand / non-protein as sticks
      safeSetStyle(
        { not: { resn: AMINO_ACIDS_ARR } },
        { stick: { colorscheme: "greenCarbon", radius: 0.2 } }
      );

    } else if (hasProtein) {
      // Protein only
      safeSetStyle(
        {},
        { cartoon: { color: "spectrum" } }
      );

    } else {
      // Ligand only
      safeSetStyle(
        {},
        { stick: { colorscheme: "greenCarbon", radius: 0.2 } }
      );
    }


    // --- Override selected helices to sticks ---
    if (selectedHelixKeys.size > 0 && helices.length > 0) {
      for (const h of helices) {
        if (!selectedHelixKeys.has(h.key)) continue;

        const start = parseInt(h.start, 10);
        const end = parseInt(h.end, 10);
        if (!Number.isFinite(start) || !Number.isFinite(end)) continue;

        const resis = [];
        for (let i = Math.min(start, end); i <= Math.max(start, end); i++) {
          resis.push(i);
        }

        safeSetStyle(
          { chain: h.chain, resi: resis },
          { stick: { radius: 0.25 } }
        );
      }
    }


    // Restore previous camera view if we have one, otherwise zoomTo once
    if (prevView) {
      viewer.setView(prevView);
    } else {
      viewer.zoomTo();
    }

    viewer.render();
    hasRenderedOnceRef.current = true;
    setTimeout(() => viewer.resize(), 0);
  };


  // main rendering effect
  useEffect(() => {
    helicesComputedRef.current = false;
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
      const viewer = viewerObjRef.current;
      if (!viewer) return;

      const prevView = hasRenderedOnceRef.current ? viewer.getView() : null;
      viewer.removeAllModels();


      fetch(`/api/mol3d?smiles=${encodeURIComponent(smiles)}`)
        .then((r) => r.text())
        .then((molData) => {
        if (!molData) return;

        // Preserve camera if we already rendered something before
        const prevView2 = hasRenderedOnceRef.current ? viewer.getView() : null;

        viewer.addModel(molData, "mol");
        viewer.setStyle({}, { stick: { colorscheme: "cyanCarbon" } });

        if (prevView2) {
          viewer.setView(prevView2);
        } else {
          viewer.zoomTo();
        }

        viewer.render();
        hasRenderedOnceRef.current = true;
        setTimeout(() => viewer.resize(), 0);
      })

        .catch(() => setErr("Failed to load 3D structure from SMILES."));
    }
  }, [pdbText, pdbFile, smiles, selectedHelixKeys, useManualMap, structureKey, manualHelixMap]);

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



      <div style={{ width: "600px", marginBottom: "8px", color: "white" }}>
        <label style={{ display: "block", marginBottom: "6px" }}>
          <input
            type="checkbox"
            checked={useManualMap}
            onChange={(e) => setUseManualMap(e.target.checked)}
            style={{ marginRight: "8px" }}
          />
          Use manual TM helix map (optional)
        </label>

        <div style={{ display: "flex", gap: "10px", alignItems: "center" }}>
          <input
            type="file"
            accept=".json,application/json"
            onChange={(e) => onManualMapFile(e.target.files?.[0] || null)}
          />

          <input
            type="text"
            placeholder='Structure key, e.g. "2YDV"'
            value={structureKey}
            onChange={(e) => setStructureKey(e.target.value)}
            style={{ width: "220px" }}
          />
        </div>

        <div style={{ color: "#aaa", fontSize: "0.85em", marginTop: "6px" }}>
          JSON should define ranges per structure key (example below).
        </div>
      </div>


      <div style={{ color: "white", background: "#333", padding: "6px", marginBottom: "6px", fontSize: "0.9em" }}>
        Helices loaded (manual): {helices.length}
      </div>



      {helices.length > 0 && (
      <div style={{ width: "600px", marginBottom: "8px" }}>
        <details>
          <summary style={{ cursor: "pointer", color: "#ccc" }}>
            Helices ({helices.length}) — show selected as sticks
          </summary>

          <div style={{ paddingLeft: "10px", marginTop: "6px" }}>
            <button onClick={clearHelices}>Clear</button>

            {helices.map((h) => (
              <div key={h.key}>
                <label>
                  <input
                    type="checkbox"
                    checked={selectedHelixKeys.has(h.key)}
                    onChange={() => toggleHelix(h.key)}
                  />
                  {h.label}
                </label>
              </div>
            ))}
          </div>
        </details>
      </div>
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
