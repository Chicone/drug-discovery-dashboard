import { useEffect, useMemo, useRef, useState } from "react";
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

  // Optional ligand inputs
  const [orthostericFile, setOrthostericFile] = useState(null);
  const [allostericPoseFile, setAllostericPoseFile] = useState(null);

  // Scenario selector
  const [scenario, setScenario] = useState("protein_only");

  // Split presets: build vs run
  const [presetBuild, setPresetBuild] = useState("m3_popc_build");
  const [presetRun, setPresetRun] = useState("m3_popc_eq");

  const [jobId, setJobId] = useState(null);
  const [status, setStatus] = useState(null);
  const logRef = useRef(null);
  const logOffsetRef = useRef(0);
  //   const [logText, setLogText] = useState("");
  const [files, setFiles] = useState([]);
  const [error, setError] = useState(null);
  const [isSubmitting, setIsSubmitting] = useState(false);

  // Recent jobs
  const [recentJobs, setRecentJobs] = useState([]);
  const [historyLimit, setHistoryLimit] = useState(20);
  const [loadingJobs, setLoadingJobs] = useState(false);
  const [selectedJobId, setSelectedJobId] = useState(null);
  const [selectedParentJobId, setSelectedParentJobId] = useState(null);

  const [environment, setEnvironment] = useState("membrane");
  const [lastBuildJobId, setLastBuildJobId] = useState(null);
  const [mdDurationNs, setMdDurationNs] = useState(50);
  const logEndRef = useRef(null);

  const [showFiles, setShowFiles] = useState(false);

  // Persistent job locking
  const [lockedJobs, setLockedJobs] = useState(() => {
    try {
      const raw = localStorage.getItem("lockedJobs");
      if (raw) return new Set(JSON.parse(raw));
    } catch (e) {
      console.error("Failed to load locked jobs:", e);
    }
    return new Set();
  });

  const presetsAll = useMemo(
    () => [
      {
        id: "m3_popc_build",
        label: "Martini 3: POPC build (martinize + insane)",
      },
      { id: "m3_popc_em", label: "Martini 3: POPC + EM" },
      { id: "m3_popc_eq", label: "Martini 3: POPC + EM + NVT + NPT" },
      { id: "m3_popc_prod", label: "Martini 3: POPC production" },
    ],
    []
  );

  const buildPresets = useMemo(
    () =>
      presetsAll.filter(
        (p) => p.id.endsWith("_build") || p.id.endsWith("_em") || p.id.endsWith("_eq")
      ),
    [presetsAll]
  );

  const runPresets = useMemo(
    () => presetsAll.filter((p) => p.id.includes("m3_popc_prod")),
    [presetsAll]
  );

  const scenarios = useMemo(
    () => [
      { id: "protein_only", label: "Protein only" },
      { id: "protein_plus_orthosteric", label: "Protein + orthosteric ligand" },
      {
        id: "protein_plus_orthosteric_plus_allosteric",
        label: "Protein + orthosteric + allosteric ligand",
      },
      { id: "dimer_two_ligands", label: "Dimer + 2 orthosteric ligands" },
    ],
    []
  );

  const stopRunJob = async () => {
    if (!jobId) return;

    try {
      await fetch(`/api/md/jobs/${jobId}/stop`, {
        method: "POST",
      });
      setStatus("stopped");
    } catch (err) {
      setError("Failed to stop MD run");
    }
  };

  async function loadRecentJobs(limit = historyLimit) {
    setLoadingJobs(true);
    try {
      const n = Number(limit);
      const q = Number.isFinite(n) && n > 0 ? `?limit=${n}` : "";
      const res = await fetch(`/api/md/jobs${q}`);
      if (!res.ok) return;

      const data = await res.json();
      setRecentJobs(Array.isArray(data) ? data : []);
    } catch (e) {
      // Non fatal
    } finally {
      setLoadingJobs(false);
    }
  }

  useEffect(() => {
    loadRecentJobs(historyLimit);
  }, [historyLimit]);

  useEffect(() => {
    const t = setInterval(() => loadRecentJobs(historyLimit), 1000);
    return () => clearInterval(t);
  }, [historyLimit]);

  useEffect(() => {
    if (scenario === "protein_only") {
      setOrthostericFile(null);
      setAllostericPoseFile(null);
    } else if (scenario === "protein_plus_orthosteric") {
      setAllostericPoseFile(null);
    }
  }, [scenario]);

  useEffect(() => {
    const saved = localStorage.getItem("lastBuildJobId");
    if (saved) {
      setLastBuildJobId(saved);
    }
  }, []);



  useEffect(() => {
    try {
      localStorage.setItem("lockedJobs", JSON.stringify([...lockedJobs]));
    } catch (e) {
      console.error("Failed to save locked jobs:", e);
    }
  }, [lockedJobs]);


  function openJob(id) {
    if (!id) return;
    setError(null);
    setSelectedParentJobId(id);
    setJobId(id);
    setStatus("running");
    logOffsetRef.current = 0;
    setFiles([]);
    pollJob(id);
  }

  function workflowFromRunPreset(preset) {
    if (!preset.includes("_prod_")) {
      throw new Error(`Invalid run preset: ${preset}`);
    }
    return "run_md";
  }

  function appendLog(chunk) {
    const el = logRef.current;
    if (!el || !chunk) return;

    el.textContent += chunk;

    // Optional: keep only last X lines
    const MAX_LINES = 2000;
    let lines = el.textContent.split("\n");
    if (lines.length > MAX_LINES) {
      lines = lines.slice(lines.length - MAX_LINES);
      el.textContent = lines.join("\n");
    }

    el.scrollTop = el.scrollHeight; // auto scroll
  }


async function submitJob({ preset, workflow, parentJobId = null }) {
  setError(null);

  if (!proteinFile) {
    setError("Please select a PDB file first.");
    return;
  }

  if (
    scenario === "protein_plus_orthosteric_plus_allosteric" &&
    !allostericPoseFile
  ) {
    setError("Scenario requires an allosteric pose file.");
    return;
  }

  if (
    scenario !== "protein_only" &&
    !orthostericFile
  ) {
    setError("Scenario requires an orthosteric ligand from a docked complex.");
    return;
  }

  setIsSubmitting(true);
  logOffsetRef.current = 0;
  setFiles([]);
  setJobId(null);
  setStatus("queued");

  try {
    const form = new FormData();

    // Always required by backend
    form.append("protein_pdb", proteinFile);

    form.append("preset", preset);
    form.append("scenario", scenario);
    form.append("environment", environment);
    form.append("workflow", workflow);
    form.append("md_ns", mdDurationNs);

    if (parentJobId) {
      form.append("parent_job_id", parentJobId);
    }

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

    // 🔥 Store build job ID
    if (workflow === "build_and_equilibrate" || workflow === "build_only") {
      setLastBuildJobId(data.job_id);
      localStorage.setItem("lastBuildJobId", data.job_id);
    }


    pollJob(data.job_id);
    loadRecentJobs(historyLimit);
  } catch (e) {
    setError(e.message || String(e));
    setStatus("error");
  } finally {
    setIsSubmitting(false);
  }
}


async function createBuildJob() {
  const isBuildOnly = presetBuild.endsWith("_build");

  return submitJob({
    preset: presetBuild,
    workflow: isBuildOnly
      ? "build_only"
      : "build_and_equilibrate",
  });
}

async function createRunJob() {
  // Prevent double runs
  if (status === "running") {
    setError("A simulation is already running.");
    return;
  }

  const parentId = selectedParentJobId || lastBuildJobId;

  if (!parentId) {
      setError("Please select a job to start MD from.");
      return;
  }

  return submitJob({
    preset: presetRun,
    workflow: "run_md",
    parentJobId: parentId,
  });
}


  async function pollJob(id) {
    const intervalMs = 1000;
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

        const lRes = await fetch(
          `/api/md/jobs/${id}/log?offset=${logOffsetRef.current}`
        );

        if (lRes.ok) {
          const chunk = await lRes.text();
          const newOffset = parseInt(
            lRes.headers.get("X-Log-Offset") ||
              `${logOffsetRef.current + chunk.length}`,
            10
          );

          if (chunk.length > 0) {
            appendLog(chunk);
          }

          logOffsetRef.current = newOffset;
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

  const canSubmit = !!proteinFile && !isSubmitting;

  return (
    <Paper
      sx={{
        p: 3,
        background: "#1e1e1e",
        maxWidth: "100%",
        overflowX: "hidden",
        color: "white",
      }}
    >
      <Typography variant="h5" gutterBottom>
        💫 Molecular Dynamics
      </Typography>

      <Typography variant="body1" sx={{ mb: 2, color: "#aaa" }}>
        Build a Martini system, then run equilibration or production presets.
      </Typography>

      <Divider sx={{ borderColor: "#333", mb: 2 }} />

      <Stack direction={{ xs: "column", md: "row" }} spacing={2}>
        {/* Left: Recent jobs */}
        <Box sx={{ minWidth: { xs: "100%", md: 360 }, maxWidth: 440 }}>
          <Typography variant="subtitle1" sx={{ mb: 1 }}>
            Recent MD jobs
          </Typography>

          <Stack direction="row" spacing={1} sx={{ mb: 1, alignItems: "center" }}>
            <TextField
              label="Limit (0 = all)"
              type="number"
              size="small"
              value={historyLimit}
              onChange={(e) => {
                const n = Number(e.target.value);
                if (!Number.isFinite(n)) return;
                setHistoryLimit(Math.max(0, Math.floor(n)));
              }}
              sx={{ width: 160 }}
            />

            <Button
              variant="outlined"
              size="small"
              onClick={() => loadRecentJobs(historyLimit)}
              disabled={loadingJobs}
            >
              Refresh
            </Button>
          </Stack>

          <Box
            sx={{
              border: "1px solid #333",
              borderRadius: 1,
              p: 1,
              background: "#111",
              maxHeight: 520,
              overflowY: "auto",
            }}
          >
            {recentJobs.length === 0 ? (
              <Typography variant="body2" sx={{ color: "#aaa" }}>
                No jobs yet.
              </Typography>
            ) : (
              <Stack spacing={1}>
                {recentJobs.map((j) => {
                      const isActive = j.job_id === selectedParentJobId;

                      return (
                        <Box
                          key={j.job_id}
                          sx={{
                            p: 1,
                            borderRadius: 1,
                            background: isActive ? "#2a2a2a" : "#1a1a1a",
                            display: "flex",
                            justifyContent: "space-between",
                            alignItems: "center",
                            cursor: "pointer",
                          }}
                        >
                          {/* Left: job info */}
                          <Box onClick={() => openJob(j.job_id)}>
                            <Typography variant="body2" sx={{ color: "#ddd" }}>
                              {j.status || "?"} | {j.scenario || "md"} | {j.preset || "preset"}
                            </Typography>

                            <Typography variant="caption" sx={{ color: "#777", display: "block" }}>
                              {j.job_id}
                            </Typography>
                          </Box>

{/* Right section: lock + delete */}
<Stack direction="row" spacing={1} alignItems="center">
  {/* Lock toggle button */}
  <Button
    size="small"
    variant="outlined"
    sx={{
      minWidth: 32,
      padding: "2px 6px",
      color: lockedJobs.has(j.job_id) ? "#66bb6a" : "#aaa",
      borderColor: lockedJobs.has(j.job_id) ? "#66bb6a" : "#555",
      "&:hover": {
        borderColor: lockedJobs.has(j.job_id) ? "#81c784" : "#888",
      },
    }}
    onClick={(e) => {
      e.stopPropagation();
      setLockedJobs(prev => {
        const next = new Set(prev);
        if (next.has(j.job_id)) next.delete(j.job_id);
        else next.add(j.job_id);
        return next;
      });
    }}

  >
    {lockedJobs.has(j.job_id) ? "🔒" : "🔓"}
  </Button>

  {/* Delete button (disabled if locked) */}
  <Button
    size="small"
    color="error"
    disabled={lockedJobs.has(j.job_id) === true}
    onClick={async (e) => {
      e.stopPropagation();
      if (lockedJobs.has(j.job_id)) return; // extra safety
      await fetch(`/api/md/jobs/${j.job_id}`, { method: "DELETE" });
      loadRecentJobs(historyLimit);
    }}
  >
    ✕
  </Button>
</Stack>

                        </Box>
                      );
                    })}
              </Stack>
            )}
          </Box>
        </Box>

        {/* Right: Build and Run */}
        <Box sx={{ flex: 1, minWidth: 0 }}>
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
                accept=".pdb,.sdf,.mol2,.pdbqt,.smi,.smiles,.txt"
                style={{ display: "none" }}
                onChange={onPickOrthosteric}
              />
              <Typography variant="body2" sx={{ color: "#aaa" }}>
                {orthostericFile
                  ? orthostericFile.name
                  : "From a docked complex to extract orthosteric pose"}
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
                Select allosteric ligand
              </Button>
              <input
                ref={allostericInputRef}
                type="file"
                accept=".pdb,.pdbqt,.sdf,.mol2,.smi,.smiles,.txt"
                style={{ display: "none" }}
                onChange={onPickAllostericPose}
              />
              <Typography variant="body2" sx={{ color: "#aaa" }}>
                {allostericPoseFile
                  ? allostericPoseFile.name
                  : "From a docked complex to extract allosteric pose"}
              </Typography>
            </Stack>

            {/* Environment */}
            <TextField
              select
              label="Environment"
              value={environment}
              onChange={(e) => setEnvironment(e.target.value)}
              size="small"
              sx={{ maxWidth: 360 }}
            >
              <MenuItem value="membrane">
                Membrane + solvent (recommended for GPCR)
              </MenuItem>
              <MenuItem value="solvated_box">Solvated box (no membrane)</MenuItem>
            </TextField>

            <Divider sx={{ borderColor: "#333" }} />

            {/* Build system */}
            <Typography variant="subtitle1">Build system</Typography>
            <TextField
              select
              label="Build preset"
              value={presetBuild}
              onChange={(e) => setPresetBuild(e.target.value)}
              size="small"
              sx={{ maxWidth: 520 }}
            >
              {buildPresets.map((p) => (
                <MenuItem key={p.id} value={p.id}>
                  {p.label}
                </MenuItem>
              ))}
            </TextField>

            <Button
              variant="outlined"
              onClick={createBuildJob}
              disabled={!canSubmit}
              sx={{ maxWidth: 220 }}
            >
              Build system
            </Button>

            <Divider sx={{ borderColor: "#333" }} />

            {/* Run MD */}
            <Typography variant="subtitle1">Run MD</Typography>
            <TextField
              label="MD duration (ns)"
              type="number"
              value={mdDurationNs}
              onChange={(e) => setMdDurationNs(Number(e.target.value))}
              size="small"
              sx={{ maxWidth: 220 }}
              inputProps={{ min: 1, step: 10 }}
            />
            <TextField
              select
              label="Run preset"
              value={presetRun}
              onChange={(e) => setPresetRun(e.target.value)}
              size="small"
              sx={{ maxWidth: 520 }}
            >
              {runPresets.map((p) => (
                <MenuItem key={p.id} value={p.id}>
                  {p.label}
                </MenuItem>
              ))}
            </TextField>

<Stack direction="row" spacing={2}>
  <Button
    variant="outlined"
    onClick={createRunJob}
    disabled={!canSubmit || status === "running"}
  >
    Run MD
  </Button>

  <Button
    variant="outlined"
    color="error"
    onClick={stopRunJob}
    disabled={!jobId || status !== "running"}
  >
    Stop MD
  </Button>

  {jobId && status === "done" && (
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
               <pre
  ref={logRef}
  style={{
    background: "#111",
    border: "1px solid #333",
    borderRadius: "6px",
    padding: "12px",
    fontFamily: "monospace",
    fontSize: "12px",
    whiteSpace: "pre-wrap",
    maxHeight: "260px",
    overflowY: "auto",
    color: "#0f0",
  }}
>
</pre>



{/* Collapsible Output Files */}
<Box sx={{ mt: 2 }}>
  <Box
    onClick={() => setShowFiles((v) => !v)}
    sx={{
      display: "flex",
      justifyContent: "space-between",
      alignItems: "center",
      cursor: "pointer",
      userSelect: "none",
      background: "#111",
      border: "1px solid #333",
      borderRadius: 1,
      px: 2,
      py: 1,
    }}
  >
    <Typography variant="subtitle2">
      Output files
    </Typography>
    <Typography sx={{ fontSize: 18 }}>
      {showFiles ? "▾" : "▸"}
    </Typography>
  </Box>

  {showFiles && (
    <Box
      sx={{
        mt: 1,
        color: "#aaa",
        fontSize: 13,
        background: "#111",
        border: "1px solid #333",
        borderRadius: 1,
        p: 2,
        maxHeight: 250,
        overflowY: "auto",
      }}
    >
      {files.length === 0 ? (
        <div>No output files yet.</div>
      ) : (
        <ul style={{ margin: 0, paddingLeft: "20px" }}>
          {files.map((f) => (
            <li key={f.name}>
              {f.name} {f.size ? `(${f.size})` : ""}
            </li>
          ))}
        </ul>
      )}
    </Box>
  )}
</Box>

            </Box>
          </Stack>
        </Box>
      </Stack>
    </Paper>
  );
}

export default MolecularDynamics;
