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
  Tooltip
} from "@mui/material";


function MolecularDynamics() {
  const proteinInputRef = useRef(null);
  const orthostericInputRef = useRef(null);
  const allostericInputRef = useRef(null);
  const [proteinFile, setProteinFile] = useState(null);
  const [runComment, setRunComment] = useState("");
  const [jobComment, setJobComment] = useState("");
  const lastLogUpdateRef = useRef(Date.now());
  const finishingRef = useRef(false);
  const finishEmptyTicksRef = useRef(0);



  // Optional ligand inputs
  const [orthostericFile, setOrthostericFile] = useState(null);
  const [allostericPoseFile, setAllostericPoseFile] = useState(null);

  // cpd5 default SMILES
  const DEFAULT_SMILES = "n1c(N2CCCCC2)c(C#N)c(c3ccccc3)c(C#N)c(N)1";
  const [orthostericSmiles, setOrthostericSmiles] = useState(DEFAULT_SMILES);
  const [ligandCase, setLigandCase] = useState("default");

  const smilesInputRef = useRef(null);

  // Scenario selector
  const [scenario, setScenario] = useState("protein_plus_orthosteric_membrane");

  // Split presets: build vs run
  const [presetBuild, setPresetBuild] = useState("m3_popc_build");
  const [presetRun, setPresetRun] = useState("m3_popc_full");

  const [jobId, setJobId] = useState(null);
  const [status, setStatus] = useState(null);
  const logRef = useRef(null);
  const logOffsetRef = useRef(0);
  const logBufferRef = useRef([]);
  const pollTimerRef = useRef(null);

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

  const [lastBuildJobId, setLastBuildJobId] = useState(null);
  const [mdDurationNs, setMdDurationNs] = useState(5);
  const [numThreads, setNumThreads] = useState(1);

  const logEndRef = useRef(null);

  const [showFiles, setShowFiles] = useState(false);

  const [lastLogUpdate, setLastLogUpdate] = useState(Date.now());


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
      {
        id: "m3_popc_full",
        label: "Martini 3: POPC — Build + EQ + Production",
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
    () => presetsAll.filter((p) => p.id.includes("prod") || p.id.includes("full")),
    [presetsAll]
  );

const scenarios = useMemo(
  () => [
    { id: "protein_only_membrane", label: "Protein only + membrane" },
    { id: "protein_only_water", label: "Protein only in water" },

    { id: "protein_plus_orthosteric_membrane", label: "Protein + orthosteric + membrane" },
    { id: "protein_plus_orthosteric_water", label: "Protein + orthosteric in water" },

    { id: "ligand_water", label: "Ligand only in water" },
  ],
  []
);

const needsProtein = scenario !== "ligand_water";

const needsOrth = useMemo(
  () =>
    scenario === "protein_plus_orthosteric_membrane" ||
    scenario === "protein_plus_orthosteric_water" ||
    scenario === "ligand_water",
  [scenario]
);

const needsAllo = useMemo(
  () => false, // not implemented yet in your scenarios list
  [scenario]
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


  useEffect(() => {
    return () => {
      if (pollTimerRef.current) {
        clearInterval(pollTimerRef.current);
      }
    };
  }, []);

useEffect(() => {
  if (!needsOrth) {
    setOrthostericFile(null);
    if (orthostericInputRef.current) {
      orthostericInputRef.current.value = "";
    }
  }

  if (!needsAllo) {
    setAllostericPoseFile(null);
    if (allostericInputRef.current) {
      allostericInputRef.current.value = "";
    }
  }
}, [needsOrth, needsAllo]);



async function openJob(id) {
  if (!id) return;

  if (pollTimerRef.current) {
    clearInterval(pollTimerRef.current);
    pollTimerRef.current = null;
  }

  setError(null);
  setSelectedParentJobId(id);
  setJobId(id);

  // Clear previous log UI
  if (logRef.current) {
    logRef.current.textContent = "… waiting for new log lines (tail mode) …\n";
  }

  logBufferRef.current = [
    "… waiting for new log lines (tail mode) …"
  ];

  finishingRef.current = false;
  finishEmptyTicksRef.current = 0;
  lastLogUpdateRef.current = Date.now();

  // ✅ Proper tail initialization (NO 10**15)
  try {
    const res = await fetch(`/api/md/jobs/${id}/log?offset=0`);
    await res.text(); // consume body (important)

    const size = parseInt(
      res.headers.get("X-Log-Size") || "0",
      10
    );

    logOffsetRef.current = size;
  } catch (e) {
    console.error("Failed to initialize tail:", e);
    logOffsetRef.current = 0;
  }

  pollJob(id);
}

async function openChimera(jobId) {
  try {
    const res = await fetch(`/api/md/${jobId}/open-chimerax`, {
      method: "POST",
    });

    const data = await res.json();
    console.log("Opened:", data.trajectory);
  } catch (e) {
    setError("Failed to launch ChimeraX");
  }
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

  const cleaned = chunk.replace(/\r/g, "");
  const newLines = cleaned.split("\n");

  // Remove placeholder if present
  if (
    logBufferRef.current.length === 1 &&
    logBufferRef.current[0].includes("waiting for new log lines")
  ) {
    logBufferRef.current = [];
  }

  logBufferRef.current.push(...newLines);

  const MAX_LINES = 2000;
  if (logBufferRef.current.length > MAX_LINES) {
    logBufferRef.current = logBufferRef.current.slice(-MAX_LINES);
  }

  el.textContent = logBufferRef.current.join("\n");
  el.scrollTop = el.scrollHeight;
}


async function submitJob({ preset, workflow, parentJobId = null }) {
  setError(null);

  const form = new FormData();

  // Protein required unless ligand-only
  if (needsProtein && !proteinFile) {
    setError("Protein PDB required for this scenario.");
    return;
  }

  // Ligand scenarios require BOTH: pose PDB + SMILES
  if (needsOrth) {
    if (!orthostericFile) {
      setError(
        scenario === "ligand_water"
          ? "Ligand pose PDB is required."
          : "Orthosteric pose PDB is required."
      );
      return;
    }

    if (!orthostericSmiles?.trim()) {
      setError(
        scenario === "ligand_water"
          ? "Ligand SMILES is required."
          : "Orthosteric SMILES is required."
      );
      return;
    }

    if (ligandCase) {
      form.append("ligand_case", ligandCase);
    }

  }


  setIsSubmitting(true);
  logOffsetRef.current = 0;
  setFiles([]);
  setJobId(null);
  setStatus("queued");

  try {

    if (needsProtein) {
      form.append("protein_pdb", proteinFile);
    }

    form.append("preset", preset);
    form.append("scenario", scenario);
    form.append("workflow", workflow);
    form.append("md_ns", mdDurationNs);
    form.append("nt", numThreads);
    if (runComment?.trim()) {
      form.append("comment", runComment.trim());
    }

    if (parentJobId) {
      form.append("parent_job_id", parentJobId);
    }

    if (needsOrth) {
      if (orthostericFile) {
        form.append("orthosteric_ligand", orthostericFile);
      }

      if (orthostericSmiles?.trim()) {
        form.append("orthosteric_smiles", orthostericSmiles.trim());
      }
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
    openJob(data.job_id);


    // 🔥 Store build job ID
    if (workflow === "build_and_equilibrate" || workflow === "build_only") {
      setLastBuildJobId(data.job_id);
      localStorage.setItem("lastBuildJobId", data.job_id);
    }

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
//   // Prevent double runs
//   if (status === "running") {
//     setError("A simulation is already running.");
//     return;
//   }


  const parentId = selectedParentJobId || lastBuildJobId;

  if (!parentId) {
      setError("Please select a job to start MD from.");
      return;
  }

  // If full preset → start fresh (no parent), full workflow
  if (presetRun === "m3_popc_full") {
    return submitJob({
      preset: presetRun,
      workflow: "build_and_run_full",   // ← NEW workflow name you handle in backend
      parentJobId: null,
    });
  }

  // Normal production-only run
  return submitJob({
    preset: presetRun,
    workflow: "run_md",
    parentJobId: parentId,
  });
}


  async function pollJob(id) {

    const intervalMs = 1000;

    // 🔥 Clear any existing polling loop
    if (pollTimerRef.current) {
      clearInterval(pollTimerRef.current);
      pollTimerRef.current = null;
    }

    async function tick() {
      let latestStatus = null;

      try {
        const sRes = await fetch(`/api/md/jobs/${id}`);
        if (sRes.ok) {
          const s = await sRes.json();
          latestStatus = s.status;
          setStatus(latestStatus);
          if (s.comment !== undefined) {
              setJobComment(s.comment);
          }
        }

      const lRes = await fetch(
        `/api/md/jobs/${id}/log?offset=${logOffsetRef.current}`
      );

      if (lRes.ok) {
        const chunk = await lRes.text();

        const newOffset = parseInt(
          lRes.headers.get("X-Log-Offset") ??
            logOffsetRef.current + chunk.length,
          10
        );

        if (chunk.length > 0) {
          appendLog(chunk);
          lastLogUpdateRef.current = Date.now();
          setLastLogUpdate(lastLogUpdateRef.current);
          finishEmptyTicksRef.current = 0;
        } else if (finishingRef.current) {
          finishEmptyTicksRef.current += 1;
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



    const isFinished =
      latestStatus === "done" || latestStatus === "error";

    if (isFinished) {
      finishingRef.current = true;
    }

    if (finishingRef.current && finishEmptyTicksRef.current >= 2) {
      if (pollTimerRef.current) {
        clearInterval(pollTimerRef.current);
        pollTimerRef.current = null;
      }
      finishingRef.current = false;
    }

  }

  await tick();
  pollTimerRef.current = setInterval(tick, intervalMs);
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

  function onPickSmilesFile(e) {
    const file = e.target.files?.[0];
    if (!file) return;

    // --- Detect ligand case from filename ---
    // Examples:
    //   zm241385.smi → ligand_case = "zm241385"
    //   cpd28.smi    → ligand_case = "cpd28"
    //   mylig.smi    → ligand_case = "mylig"
    const name = file.name.toLowerCase();

    const caseKey = name
      .replace(/\.[^/.]+$/, "") // remove extension
      .trim()
      .replace(/\s/g, "_");     // replace spaces for safety

    setLigandCase(caseKey || "default");

    // --- Read SMILES ---
    const reader = new FileReader();
    reader.onload = (event) => {
      const text = event.target.result || "";
      const lines = text
        .split(/\r?\n/)
        .map((l) => l.trim())
        .filter(Boolean);

      if (lines.length === 0) return;

      const smiles = lines[0].split(/\s+/)[0];
      setOrthostericSmiles(smiles);
    };

    reader.readAsText(file);
  }


  const canSubmit =
    (!needsProtein || !!proteinFile) &&
    (!needsOrth || (!!orthostericFile && !!orthostericSmiles?.trim())) &&
    (!needsAllo || !!allostericPoseFile) &&
    !isSubmitting;

  const runningCount = useMemo(() => {
    return recentJobs.filter(j => j.status === "running").length;
  }, [recentJobs]);


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
      <Box sx={{ minWidth: { xs: "100%", md: 380 }, maxWidth: 460 }}>
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

              <Typography
                variant="body2"
                sx={{
                  ml: 1,
                  fontWeight: 600,
                  color: runningCount > 0 ? "#66bb6a" : "#888",
                  fontFamily: "monospace",
                }}
              >
                {runningCount} running
              </Typography>
            </Stack>
        <Box
                sx={{
                  border: "1px solid #333",
                  borderRadius: 1,
                  p: 1,
                  background: "#111",
                  maxHeight: "70vh",
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
                    <Box
                      onClick={() => openJob(j.job_id)}
                      sx={{
                        flex: 1,
                        minWidth: 0,
                        overflow: "hidden",
                      }}
                    >
                        <Tooltip
                          title={`${j.status || "?"} | ${j.scenario || "md"} | ${j.preset || "preset"}`}
                          arrow
                        >
                          <Typography
                            variant="body2"
                            sx={{
                              color: "#ddd",
                              whiteSpace: "nowrap",
                              overflow: "hidden",
                              textOverflow: "ellipsis",
                              cursor: "default",
                            }}
                          >
                            {j.status || "?"} | {j.scenario || "md"} | {j.preset || "preset"}
                          </Typography>
                        </Tooltip>

                        <Tooltip title={j.job_id} arrow>
                          <Typography
                            variant="body2"
                            sx={{
                              fontFamily: "monospace",
                              fontSize: "0.75rem",
                              fontWeight: 600,
                              color: "#66bb6a",
                              letterSpacing: 0.5,
                              whiteSpace: "nowrap",
                              overflow: "hidden",
                              textOverflow: "ellipsis",
                              maxWidth: 180,   // adjust if needed
                              cursor: "default",
                            }}
                          >
                            {j.job_id}
                          </Typography>
                        </Tooltip>


                    </Box>

                    {/* Right: lock + delete */}
                    <Stack direction="row" spacing={1} alignItems="center">
                      {/* Chimera launching button */}
                      <Button
                          size="small"
                          variant="outlined"
                          sx={{
                            minWidth: 32,
                            padding: "2px 6px",
                            borderColor: "#555",
                            color: "#66bb6a",
                            "&:hover": { borderColor: "#81c784" },
                          }}
                          onClick={(e) => {
                            e.stopPropagation();
                            openChimera(j.job_id);
                          }}
                        >
                          🧬
                        </Button>


                      <Button
                        size="small"
                        variant="outlined"
                        sx={{
                          minWidth: 32,
                          padding: "2px 6px",
                          color: lockedJobs.has(j.job_id) ? "#66bb6a" : "#aaa",
                          borderColor: lockedJobs.has(j.job_id)
                            ? "#66bb6a"
                            : "#555",
                          "&:hover": {
                            borderColor: lockedJobs.has(j.job_id)
                              ? "#81c784"
                              : "#888",
                          },
                        }}
                        onClick={(e) => {
                          e.stopPropagation();
                          setLockedJobs((prev) => {
                            const next = new Set(prev);
                            if (next.has(j.job_id)) next.delete(j.job_id);
                            else next.add(j.job_id);
                            return next;
                          });
                        }}
                      >
                        {lockedJobs.has(j.job_id) ? "🔒" : "🔓"}
                      </Button>

                      <Button
                        size="small"
                        color="error"
                        disabled={lockedJobs.has(j.job_id) === true}
                        onClick={async (e) => {
                          e.stopPropagation();
                          if (lockedJobs.has(j.job_id)) return;
                          await fetch(`/api/md/jobs/${j.job_id}`, {
                            method: "DELETE",
                          });
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

          {/* Buttons row */}
          <Stack
            direction="row"
            spacing={2}
            alignItems="center"
            sx={{ flexWrap: "wrap" }}
          >
        <Button
          variant="outlined"
          color="secondary"
          onClick={() => proteinInputRef.current?.click()}
          disabled={!needsProtein}
          sx={{ minWidth: 110 }}
        >
          Receptor
        </Button>


            <Button
              variant="outlined"
              color="secondary"
              onClick={() => orthostericInputRef.current?.click()}
              disabled={!needsOrth}
              sx={{ minWidth: 110 }}
            >
              {scenario === "ligand_water" ? "Ligand" : "Ortho"}
            </Button>


            <Button
                  variant="outlined"
                  color="secondary"
                  onClick={() => allostericInputRef.current?.click()}
                  disabled={!needsAllo}
                  sx={{ minWidth: 110 }}
                >
                  Allo
                </Button>
          {/* Selected files under buttons */}
          <Box>
            {proteinFile && (
              <Typography variant="caption" sx={{ color: "#aaa", display: "block" }}>
                PDB: {proteinFile.name}
              </Typography>
            )}
            {orthostericFile && (
              <Typography variant="caption" sx={{ color: "#aaa", display: "block" }}>
                Ortho: {orthostericFile.name}
              </Typography>
            )}
            {allostericPoseFile && (
              <Typography variant="caption" sx={{ color: "#aaa", display: "block" }}>
                Allo: {allostericPoseFile.name}
              </Typography>
            )}
          </Box>

            <input
              ref={proteinInputRef}
              type="file"
              accept=".pdb,.ent"
              style={{ display: "none" }}
              onChange={onPickProtein}
            />

            <input
              ref={orthostericInputRef}
              type="file"
              accept=".pdb,.sdf,.mol2,.pdbqt,.smi,.smiles,.txt"
              style={{ display: "none" }}
              onChange={onPickOrthosteric}
            />

            <input
              ref={allostericInputRef}
              type="file"
              accept=".pdb,.pdbqt,.sdf,.mol2,.smi,.smiles,.txt"
              style={{ display: "none" }}
              onChange={onPickAllostericPose}
            />

            <input
              ref={smilesInputRef}
              type="file"
              accept=".smi,.smiles,.txt"
              style={{ display: "none" }}
              onChange={onPickSmilesFile}
            />



          </Stack>

            {needsOrth && (
  <>
        <Button
      variant="outlined"
      color="secondary"
      size="small"
      onClick={() => smilesInputRef.current?.click()}
      sx={{ maxWidth: 200, mt: 1 }}
    >
      Load SMILES (.smi)
    </Button>

    <TextField
      label={
        scenario === "ligand_water"
          ? "Ligand SMILES"
          : "Orthosteric SMILES (required)"
      }
      value={orthostericSmiles}
      onChange={(e) => setOrthostericSmiles(e.target.value)}
      size="small"
      fullWidth
      sx={{ maxWidth: 520 }}
    />

  </>
)}







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
            sx={{ maxWidth: 200 }}
          >
            Build system
          </Button>

          <Divider sx={{ borderColor: "#333" }} />

          {/* Run MD */}
          <Typography variant="subtitle1">Run MD</Typography>

          <Stack direction="row" spacing={2} alignItems="center">
            <TextField
              label="Duration (ns)"
              type="number"
              value={mdDurationNs}
              onChange={(e) => setMdDurationNs(Number(e.target.value))}
              size="small"
              sx={{ width: 110 }}
              inputProps={{ min: 1, step: 10 }}
            />

            <TextField
              label="CPUs"
              type="number"
              value={numThreads}
              onChange={(e) => setNumThreads(Number(e.target.value))}
              size="small"
              sx={{ width: 90 }}
              inputProps={{ min: 1, step: 1 }}
            />

            <TextField
              select
              label="Preset"
              value={presetRun}
              onChange={(e) => setPresetRun(e.target.value)}
              size="small"
              sx={{ minWidth: 290 }}
            >
              {runPresets.map((p) => (
                <MenuItem key={p.id} value={p.id}>
                  {p.label}
                </MenuItem>
              ))}
            </TextField>
          </Stack>

          <Stack direction="row" spacing={2}>
            <Button
              variant="outlined"
              onClick={createRunJob}
              disabled={!canSubmit}
//               disabled={!canSubmit || status === "running"}
              sx={{ minWidth: 110 }}
            >
              Run MD
            </Button>
            <TextField
              label="Run Comment"
              fullWidth
              multiline
              rows={2}
              value={runComment}
              onChange={(e) => setRunComment(e.target.value)}
              sx={{ mt: 2 }}
            />


            <Button
              variant="outlined"
              color="error"
              onClick={stopRunJob}
              disabled={!jobId || status !== "running"}
              sx={{ minWidth: 110 }}
            >
              Stop
            </Button>

            {jobId && status === "done" && (
              <Button
                variant="contained"
                color="primary"
                href={`/api/md/jobs/${jobId}/download`}
              >
                Download
              </Button>
            )}
          </Stack>

          {error && (
            <Typography variant="body2" sx={{ color: "#ff8080" }}>
              {error}
            </Typography>
          )}

          <Divider sx={{ borderColor: "#333" }} />

          {jobComment && (
            <Typography variant="body2" sx={{ color: "#ccc", mb: 1 }}>
              Comment: {jobComment}
            </Typography>
          )}

          {/* Logs */}
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
                height: "260px",
                maxHeight: "260px",
                overflowY: "auto",
                color: "#0f0",
              }}
            />
          </Box>
        </Stack>
      </Box>
    </Stack>
  </Paper>
);
}

export default MolecularDynamics;
