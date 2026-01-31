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
  const [logText, setLogText] = useState("");
  const logOffsetRef = useRef(0);
  const [files, setFiles] = useState([]);
  const [error, setError] = useState(null);
  const [isSubmitting, setIsSubmitting] = useState(false);

  // Recent jobs
  const [recentJobs, setRecentJobs] = useState([]);
  const [historyLimit, setHistoryLimit] = useState(20);
  const [loadingJobs, setLoadingJobs] = useState(false);

  const [environment, setEnvironment] = useState("membrane");
  const [lastBuildJobId, setLastBuildJobId] = useState(null);

  const presetsAll = useMemo(
    () => [
      {
        id: "m3_popc_build",
        label: "Martini 3: POPC build (martinize + insane)",
      },
      { id: "m3_popc_em", label: "Martini 3: POPC + EM" },
      { id: "m3_popc_eq", label: "Martini 3: POPC + EM + NVT + NPT" },
      { id: "m3_popc_prod_50ns", label: "Martini 3: POPC production 50 ns" },
      { id: "m3_popc_prod_200ns", label: "Martini 3: POPC production 200 ns" },
      {
        id: "m3_popc_chol_prod_200ns",
        label: "Martini 3: POPC+CHOL production 200 ns",
      },
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
    () => presetsAll.filter((p) => p.id.includes("_prod_")),
    [presetsAll]
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
    const t = setInterval(() => loadRecentJobs(historyLimit), 5000);
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

  function openJob(id) {
    if (!id) return;
    setError(null);
    setJobId(id);
    setStatus("running");
    setLogText("");
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

  setIsSubmitting(true);
  setLogText("");
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
    if (workflow === "build_and_equilibrate") {
      setLastBuildJobId(data.job_id);
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
  if (!lastBuildJobId) {
    setError("Please build and equilibrate the system first.");
    return;
  }

  return submitJob({
    preset: presetRun,
    workflow: "run_md",
    parentJobId: lastBuildJobId,
  });
}


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
            setLogText((prev) => prev + chunk);
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
                      const isActive = j.job_id === jobId;

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

                          {/* Right: delete button */}
                          <Button
                            size="small"
                            color="error"
                            onClick={async (e) => {
                              e.stopPropagation(); // VERY IMPORTANT
                              await fetch(`/api/md/jobs/${j.job_id}`, { method: "DELETE" });
                              loadRecentJobs(historyLimit);
                            }}
                          >
                            ✕
                          </Button>
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
                  : "Optional. Backend may extract it from receptor PDB."}
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
                accept=".pdb,.pdbqt,.sdf,.mol2,.smi,.smiles,.txt"
                style={{ display: "none" }}
                onChange={onPickAllostericPose}
              />
              <Typography variant="body2" sx={{ color: "#aaa" }}>
                {allostericPoseFile
                  ? allostericPoseFile.name
                  : "Optional (depends on scenario)"}
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
                disabled={!canSubmit}
              >
                Run MD
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
        </Box>
      </Stack>
    </Paper>
  );
}

export default MolecularDynamics;
