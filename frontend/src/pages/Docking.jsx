import { useState, useEffect } from "react";
import {
  Box,
  Paper,
  Typography,
  Button,
  Stack,
  LinearProgress,
  Divider,
    Tooltip,
} from "@mui/material";

import MoleculeViewer from "../MoleculeViewer";
import Grid from "@mui/material/Grid";
import TextField from "@mui/material/TextField";

function Docking() {
  const [residueView, setResidueView] = useState("distance"); // "distance" | "sequence"
  const [selectedRunId, setSelectedRunId] = useState(null);
  const [poseAnalysis, setPoseAnalysis] = useState(null);
  const [runHistory, setRunHistory] = useState([]);
  const [receptorFile, setReceptorFile] = useState(null);
  const [ligandFile, setLigandFile] = useState(null);
  const [pdbText, setPdbText] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState("");
  const [poses, setPoses] = useState([]);
  const [selectedPoseIdx, setSelectedPoseIdx] = useState(0);
  const [boxFile, setBoxFile] = useState(null);
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

  const [historyLimit, setHistoryLimit] = useState(1);


  const updateDockParam = (key, value) => {
    setDockParams((p) => ({ ...p, [key]: value }));
  };

    const applyBoxFromJsonText = (jsonText) => {
      const data = JSON.parse(jsonText);

      const c = data.box_center;
      const s = data.box_size;

      if (!Array.isArray(c) || c.length !== 3 || !Array.isArray(s) || s.length !== 3) {
        throw new Error("Invalid box JSON: expected box_center and box_size arrays of length 3.");
      }

      setDockParams((p) => ({
        ...p,
        center_x: Number(c[0]),
        center_y: Number(c[1]),
        center_z: Number(c[2]),
        size_x: Number(s[0]),
        size_y: Number(s[1]),
        size_z: Number(s[2]),
      }));
    };

    const loadRunHistory = async (limit = 1) => {
      try {
        const q = Number(limit) > 0 ? `?limit=${Number(limit)}` : "";
        const res = await fetch(`/api/docking/runs${q}`);
        if (!res.ok) {
          console.error("Failed to fetch run history:", await res.text());
          return;
        }
        const data = await res.json();
        setRunHistory(data);
      } catch (err) {
        console.error("Failed to load run history:", err);
      }
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
    loadRunHistory(historyLimit);


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

    useEffect(() => {
      loadRunHistory(historyLimit);
    }, [historyLimit]);



    const loadRun = async (runId) => {
      try {
        const res = await fetch(`/api/docking/runs/${runId}`);
        if (!res.ok) return;

        setSelectedRunId(runId);

        const data = await res.json();

        // 1) poses (this part already works for you)
        if (data.poses?.length) {
          setPoses(data.poses);
          setSelectedPoseIdx(0);
          setPdbText(data.poses[0].pdb);
        }

        setSelectedPoseIdx(0);
        setPdbText(data.poses[0].pdb);
        loadPoseAnalysis(runId, 0);


        // 2) run metadata (box + params)
        const run = data.run || {};

        // Box
        const box = run.box;
        if (box?.center?.length === 3 && box?.size?.length === 3) {
          setDockParams((p) => ({
            ...p,
            center_x: Number(box.center[0]),
            center_y: Number(box.center[1]),
            center_z: Number(box.center[2]),
            size_x: Number(box.size[0]),
            size_y: Number(box.size[1]),
            size_z: Number(box.size[2]),
          }));
        }

        // Vina params (support BOTH possible names)
        const vina = run.vina || run.vina_params || {};
        if (vina.exhaustiveness != null || vina.num_modes != null || vina.cpu != null) {
          setDockParams((p) => ({
            ...p,
            exhaustiveness: vina.exhaustiveness != null
              ? Number(vina.exhaustiveness)
              : p.exhaustiveness,
            num_modes: vina.num_modes != null
              ? Number(vina.num_modes)
              : p.num_modes,
            cpu: vina.cpu != null ? Number(vina.cpu) : p.cpu,
          }));
        }
      } catch (err) {
        console.error("Failed to load run:", err);
      }
    };

  const loadPoseAnalysis = async (runId, poseIdx) => {
    if (!runId && runId !== 0) return;

    try {
      const res = await fetch(
        `/api/docking/runs/${runId}/analyze/${poseIdx}?cutoff=4.0`
      );
      if (!res.ok) return;
      const data = await res.json();
      setPoseAnalysis(data);
    } catch (e) {
      console.error("Failed to load pose analysis:", e);
    }
  };

    const groupResiduesIntoRanges = (residuesSeq) => {
    // residuesSeq items: { chain, res_seq, res_name, min_dist }
    // We group contiguous residue numbers per chain into ranges.

    const toInt = (x) => {
      const n = Number.parseInt(String(x), 10);
      return Number.isFinite(n) ? n : null;
    };

    // Group by chain
    const byChain = new Map();
    (residuesSeq ?? []).forEach((r) => {
      const chain = r.chain ?? "?";
      if (!byChain.has(chain)) byChain.set(chain, []);
      byChain.get(chain).push(r);
    });

    const ranges = [];

    for (const [chain, list] of byChain.entries()) {
      // Sort by residue index numeric
      const sorted = [...list].sort((a, b) => {
        const ai = toInt(a.res_seq);
        const bi = toInt(b.res_seq);
        if (ai === null && bi === null) return String(a.res_seq).localeCompare(String(b.res_seq));
        if (ai === null) return 1;
        if (bi === null) return -1;
        return ai - bi;
      });

      let start = null;
      let prev = null;

      const flush = () => {
        if (!start) return;
        const sName = start.res_name;
        const sSeq = start.res_seq;
        const eName = prev.res_name;
        const eSeq = prev.res_seq;

        const label =
          String(sSeq) === String(eSeq)
            ? `${sName}${sSeq} (${chain})`
            : `${sName}${sSeq}–${eName}${eSeq} (${chain})`;

        ranges.push({
          key: `${chain}:${sSeq}-${eSeq}`,
          label,
          chain,
          start_seq: sSeq,
          end_seq: eSeq,
        });

        start = null;
        prev = null;
      };

      for (const r of sorted) {
        const i = toInt(r.res_seq);
        if (i === null) {
          // Non-integer residue id: treat as its own "range"
          if (start) flush();
          ranges.push({
            key: `${chain}:${r.res_seq}-${r.res_seq}`,
            label: `${r.res_name}${r.res_seq} (${chain})`,
            chain,
            start_seq: r.res_seq,
            end_seq: r.res_seq,
          });
          continue;
        }

        if (!start) {
          start = r;
          prev = r;
          continue;
        }

        const prevI = toInt(prev.res_seq);
        if (prevI !== null && i === prevI + 1) {
          // contiguous
          prev = r;
        } else {
          flush();
          start = r;
          prev = r;
        }
      }

      flush();
    }

    // Keep chains grouped but stable
    return ranges;
  };


  const residuesToShow =
  residueView === "sequence"
    ? (poseAnalysis?.residues_seq ?? [])
    : (poseAnalysis?.residues ?? []);

  const residueRanges =
    residueView === "sequence"
      ? groupResiduesIntoRanges(poseAnalysis?.residues_seq ?? [])
      : [];


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
  accept=".smi,.pdb,.sdf,.mol2"
  hidden
  onChange={async (e) => {
    const file = e.target.files?.[0];
    if (!file) return;

    setLigandFile(file);

    // 👇 THIS is where /api/load_ligand is called
    const fd = new FormData();
    fd.append("ligand", file);

    try {
      const res = await fetch("/api/load_ligand", {
        method: "POST",
        body: fd,
      });

      if (!res.ok) {
        console.error("load_ligand failed");
        return;
      }

      const data = await res.json();
      console.log("load_ligand response:", data);

      // 👇 THIS is where the form is auto-filled
      if (data.box?.center && data.box?.size) {
        setDockParams((p) => ({
          ...p,
          center_x: data.box.center[0],
          center_y: data.box.center[1],
          center_z: data.box.center[2],
          size_x: data.box.size[0],
          size_y: data.box.size[1],
          size_z: data.box.size[2],
        }));
      }
    } catch (err) {
      console.error("Failed to load ligand:", err);
    }
  }}
/>

        </Button>
      </Stack>

      <Divider sx={{ my: 2 }} />

      {/* Parameters */}
<Typography variant="subtitle1" sx={{ mb: 0.5 }}>
  Parameters
</Typography>

<Box sx={{ maxWidth: 120, mx: "left" }}>
  <Stack spacing={1}>
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
      label="Modes"
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
      label="CPUs (0 = auto)"
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

<Box sx={{ mb: 1.5 }}>
 <Tooltip title="Load, if available, a binding-site box from a box.json file to auto-fill the docking box fields" arrow>
  <Button
    variant="contained"
    component="label"
    color="secondary"
    fullWidth
    size="small"
    sx={{ py: 0.8 }}
  >
    Load docking box (JSON)
    <input
      type="file"
      accept=".json"
      hidden
      onChange={(e) => {
        const f = e.target.files?.[0];
        if (!f) return;

        const reader = new FileReader();
        reader.onload = () => {
          try {
            applyBoxFromJsonText(String(reader.result));
          } catch (err) {
            console.error(err);
            setError(err.message);
          }
        };
        reader.readAsText(f);
      }}
    />
  </Button>
 </Tooltip>
</Box>


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

  {poseAnalysis && (
    <Box sx={{ mt: 2 }}>
      <Typography variant="subtitle1" sx={{ mb: 1 }}>
        Pose interactions (within 4 Å)
      </Typography>

      <Typography variant="body2" sx={{ color: "#aaa" }}>
        H-bond candidates: {poseAnalysis.counts?.hbond_candidates ?? 0}
        {" | "}
        Hydrophobic: {poseAnalysis.counts?.hydrophobic ?? 0}
        {" | "}
        Polar: {poseAnalysis.counts?.polar ?? 0}
      </Typography>

      <Box sx={{ display: "flex", gap: 1, alignItems: "center", mt: 1 }}>
        <Typography variant="body2" sx={{ color: "#aaa" }}>
          Residues view:
        </Typography>

        <Button
          size="small"
          variant={residueView === "distance" ? "contained" : "outlined"}
          onClick={() => setResidueView("distance")}
        >
          Distance
        </Button>

        <Button
          size="small"
          variant={residueView === "sequence" ? "contained" : "outlined"}
          onClick={() => setResidueView("sequence")}
        >
          Sequence
        </Button>
      </Box>

      <Box sx={{ mt: 1, display: "flex", flexWrap: "wrap", gap: 1 }}>
        {residueView === "sequence"
          ? residueRanges.map((g) => (
              <Box
                key={g.key}
                sx={{
                  px: 1,
                  py: 0.5,
                  borderRadius: 1,
                  background: "#2a2a2a",
                  fontSize: "0.85rem",
                }}
              >
                {g.label}
              </Box>
            ))
          : residuesToShow.map((r) => (
              <Box
                key={`${r.chain}-${r.res_seq}-${r.res_name}`}
                sx={{
                  px: 1,
                  py: 0.5,
                  borderRadius: 1,
                  background: "#2a2a2a",
                  fontSize: "0.85rem",
                }}
              >
                {r.res_name}{r.res_seq} ({r.chain})
                {typeof r.min_dist === "number" ? ` – ${r.min_dist} Å` : ""}
              </Box>
            ))}
      </Box>

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
              if (selectedRunId) {
                loadPoseAnalysis(selectedRunId, i);
              }
            }}
          >
            Pose {pose.mode}
            {typeof pose.score === "number" ? ` (${pose.score.toFixed(2)})` : ""}
          </Button>
        ))}
      </Box>
    </Box>
  )}

  {runHistory.length > 0 && (
      <Box sx={{ mt: 3 }}>
        <Box sx={{ display: "flex", alignItems: "center", gap: 1, mb: 1 }}>
          <Typography variant="subtitle1" sx={{ flex: 1 }}>
            Docking run history
          </Typography>

          <TextField
            label="History entries"
            type="number"
            size="small"
            value={historyLimit}
            inputProps={{ min: 1 }}
            sx={{ width: 140 }}
            onChange={(e) => {
              const n = Number(e.target.value);
              if (!Number.isFinite(n)) return;
              setHistoryLimit(Math.max(1, Math.floor(n)));
            }}
          />
        </Box>


        <Stack spacing={0.5}>
          {runHistory.map((run) => (
            <Box
              key={run.run_id}
              sx={{
                p: 1,
                borderRadius: 1,
                background: "#2a2a2a",
                fontSize: "0.85rem",
                cursor: "pointer",
                "&:hover": { background: "#333" },
              }}
              onClick={() => loadRun(run.run_id)}
            >

              <strong>{run.ligand}</strong>
              {typeof run.best_score === "number" && (
                <> — best {run.best_score.toFixed(2)}</>
              )}
              <br />
              <span style={{ color: "#aaa" }}>{run.created_at}</span>
            </Box>
          ))}
        </Stack>
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
