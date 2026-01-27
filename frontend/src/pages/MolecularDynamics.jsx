import { useMemo, useRef, useState } from "react";
import {
  Paper,
  Typography,
  Button,
  Stack,
  Box,
  TextField,
  MenuItem,
  Divider,
} from "@mui/material";

function MolecularDynamics() {
  const proteinInputRef = useRef(null);
  const orthostericInputRef = useRef(null);
  const allostericInputRef = useRef(null);

  const [proteinFile, setProteinFile] = useState(null);

  // New: optional ligand inputs
  const [orthostericFile, setOrthostericFile] = useState(null);
  const [allostericPoseFile, setAllostericPoseFile] = useState(null);

  // New: scenario selector
  const [scenario, setScenario] = useState("protein_only");

  const [preset, setPreset] = useState("cg_popc_50ns");
  const [jobId, setJobId] = useState(null);
  const [status, setStatus] = useState(null);
  const [logText, setLogText] = useState("");
  const [files, setFiles] = useState([]);
  const [error, setError] = useState(null);
  const [isSubmitting, setIsSubmitting] = useState(false);

  const presets = useMemo(
    () => [
      { id: "cg_popc_50ns", label: "CG Martini: POPC, 50 ns" },
      { id: "cg_popc_200ns", label: "CG Martini: POPC, 200 ns" },
      { id: "cg_popc_chol_50ns", label: "CG Martini: POPC+CHOL, 50 ns" },
    ],
    []
  );

  const scenarios = useMemo(
    () => [
      { id: "protein_only", label: "Protein only" },
      { id: "protein_plus_orthosteric", label: "Protein + orthosteric ligand" },
      {
        id: "protein_plus_orthosteric_plus_allosteric",
        label: "Protein + orthosteric + allosteric (pose)",
      },
    ],
    []
  );

  async function createJob() {
    setError(null);

    if (!proteinFile) {
      setError("Please select a PDB file first.");
      return;
    }

    // Basic validation for scenarios
    if (
      scenario === "protein_plus_orthosteric" &&
      !orthostericFile
    ) {
      setError("Scenario requires an orthosteric ligand file.");
      return;
    }

    if (
      scenario === "protein_plus_orthosteric_plus_allosteric" &&
      (!orthostericFile || !allostericPoseFile)
    ) {
      setError(
        "Scenario requires both an orthosteric ligand file and an allosteric pose file."
      );
      return;
    }

    setIsSubmitting(true);
    setLogText("");
    setFiles([]);
    setJobId(null);
    setStatus("queued");

    try {
      const form = new FormData();
      form.append("protein_pdb", proteinFile);
      form.append("preset", preset);

      // New fields
      form.append("scenario", scenario);

      if (orthostericFile) {
        form.append("orthosteric_ligand", orthostericFile);
      }
      if (allostericPoseFile) {
        form.append("allosteric_pose", allostericPoseFile);
      }

      const res = await fetch("/api/md/jobs", {
        method: "POST",
        body: form,
      });

      if (!res.ok) {
        const txt = await res.text();
        throw new Error(txt || `Failed to create job (${res.status})`);
      }

      const data = await res.json();
      setJobId(data.job_id);
      setStatus("running");

      // Start polling status and logs
      pollJob(data.job_id);
    } catch (e) {
      setError(e.message || String(e));
      setStatus("error");
    } finally {
      setIsSubmitting(false);
    }
  }

  // FIX: stop polling based on the latest fetched status, not stale React state
  async function pollJob(id) {
    const intervalMs = 2000;
    let timer = null;

    async function tick() {
      let latestStatus = null;

      try {
        const sRes = await fetch(`/api/md/jobs/${id}`);
        if (sRes.ok) {
          const s = await sRes.json();
          latestStatus = s.status;
          setStatus(latestStatus);
        }

        const lRes = await fetch(`/api/md/jobs/${id}/log`);
        if (lRes.ok) {
          const txt = await lRes.text();
          setLogText(txt);
        }

        const fRes = await fetch(`/api/md/jobs/${id}/files`);
        if (fRes.ok) {
          const f = await fRes.json();
          setFiles(f.files || []);
        }
      } catch (e) {
        setError((prev) => prev || (e.message || String(e)));
      }

      if (latestStatus === "done" || latestStatus === "error") {
        if (timer) clearInterval(timer);
      }
    }

    await tick();
    timer = setInterval(tick, intervalMs);
  }

  function onPickProtein(e) {
    const f = e.target.files?.[0] || null;
    setProteinFile(f);
  }

  function onPickOrthosteric(e) {
    const f = e.target.files?.[0] || null;
    setOrthostericFile(f);
  }

  function onPickAllostericPose(e) {
    const f = e.target.files?.[0] || null;
    setAllostericPoseFile(f);
  }

  const canRun = !!proteinFile && !isSubmitting;

  return (
    <Paper
      sx={{
        p: 3,
        background: "#1e1e1e",
        maxWidth: "100%",
        overflowX: "hidden",
      }}
    >
      <Typography variant="h5" gutterBottom>
        💫 Molecular Dynamics
      </Typography>

      <Typography variant="body1" sx={{ mb: 2 }}>
        Upload a protein PDB, choose a preset and scenario, and launch a CG MD job
        (GROMACS + Martini in Docker).
      </Typography>

      <Stack spacing={2}>
        {/* Protein */}
        <Stack direction="row" spacing={2} alignItems="center">
          <Button
            variant="contained"
            color="secondary"
            onClick={() => proteinInputRef.current?.click()}
          >
            Select PDB
          </Button>
          <input
            ref={proteinInputRef}
            type="file"
            accept=".pdb,.ent"
            style={{ display: "none" }}
            onChange={onPickProtein}
          />
          <Typography variant="body2" sx={{ color: "#aaa" }}>
            {proteinFile ? proteinFile.name : "No file selected"}
          </Typography>
        </Stack>

        {/* Scenario */}
        <TextField
          select
          label="Scenario"
          value={scenario}
          onChange={(e) => setScenario(e.target.value)}
          size="small"
          sx={{ maxWidth: 520 }}
        >
          {scenarios.map((s) => (
            <MenuItem key={s.id} value={s.id}>
              {s.label}
            </MenuItem>
          ))}
        </TextField>

        {/* Orthosteric ligand */}
        <Stack direction="row" spacing={2} alignItems="center">
          <Button
            variant="outlined"
            color="secondary"
            onClick={() => orthostericInputRef.current?.click()}
            disabled={scenario === "protein_only"}
          >
            Select orthosteric ligand
          </Button>
          <input
            ref={orthostericInputRef}
            type="file"
            accept=".pdb,.sdf,.mol2,.pdbqt"
            style={{ display: "none" }}
            onChange={onPickOrthosteric}
          />
          <Typography variant="body2" sx={{ color: "#aaa" }}>
            {orthostericFile ? orthostericFile.name : "Optional (depends on scenario)"}
          </Typography>
        </Stack>

        {/* Allosteric pose */}
        <Stack direction="row" spacing={2} alignItems="center">
          <Button
            variant="outlined"
            color="secondary"
            onClick={() => allostericInputRef.current?.click()}
            disabled={scenario !== "protein_plus_orthosteric_plus_allosteric"}
          >
            Select allosteric pose
          </Button>
          <input
            ref={allostericInputRef}
            type="file"
            accept=".pdb,.pdbqt,.sdf,.mol2"
            style={{ display: "none" }}
            onChange={onPickAllostericPose}
          />
          <Typography variant="body2" sx={{ color: "#aaa" }}>
            {allostericPoseFile ? allostericPoseFile.name : "Optional (depends on scenario)"}
          </Typography>
        </Stack>

        {/* Preset */}
        <TextField
          select
          label="Preset"
          value={preset}
          onChange={(e) => setPreset(e.target.value)}
          size="small"
          sx={{ maxWidth: 360 }}
        >
          {presets.map((p) => (
            <MenuItem key={p.id} value={p.id}>
              {p.label}
            </MenuItem>
          ))}
        </TextField>

        <Stack direction="row" spacing={2}>
          <Button variant="outlined" onClick={createJob} disabled={!canRun}>
            Run simulation
          </Button>

          {jobId && (
            <Button
              variant="contained"
              color="primary"
              href={`/api/md/jobs/${jobId}/download`}
            >
              Download results
            </Button>
          )}
        </Stack>

        {error && (
          <Typography variant="body2" sx={{ color: "#ff8080" }}>
            {error}
          </Typography>
        )}

        <Divider sx={{ borderColor: "#333" }} />

        <Box>
          <Typography variant="body2" sx={{ color: "#aaa", mb: 1 }}>
            Status: {status || "idle"} {jobId ? `(job: ${jobId})` : ""}
          </Typography>

          <Typography variant="subtitle2" sx={{ mb: 1 }}>
            Logs
          </Typography>
          <Box
            sx={{
              background: "#111",
              border: "1px solid #333",
              borderRadius: 1,
              p: 2,
              fontFamily: "monospace",
              fontSize: 12,
              whiteSpace: "pre-wrap",
              maxHeight: 260,
              overflowY: "auto",
            }}
          >
            {logText || "No logs yet."}
          </Box>

          <Typography variant="subtitle2" sx={{ mt: 2, mb: 1 }}>
            Output files
          </Typography>
          <Box sx={{ color: "#aaa", fontSize: 13 }}>
            {files.length === 0 ? (
              <div>No output files yet.</div>
            ) : (
              <ul>
                {files.map((f) => (
                  <li key={f.name}>
                    {f.name} {f.size ? `(${f.size})` : ""}
                  </li>
                ))}
              </ul>
            )}
          </Box>
        </Box>
      </Stack>
    </Paper>
  );
}

export default MolecularDynamics;
