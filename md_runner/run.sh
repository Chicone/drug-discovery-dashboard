#!/usr/bin/env bash
set -euo pipefail

JOB_DIR="/job"
IN_DIR="${JOB_DIR}/input"
OUT_DIR="${JOB_DIR}/out"
WORK_DIR="${JOB_DIR}/work"
LOG="${JOB_DIR}/log.txt"
STATUS="${JOB_DIR}/status.json"

mkdir -p "${OUT_DIR}" "${WORK_DIR}"
echo '{"status":"running"}' > "${STATUS}"
: > "${LOG}"

pdb_in="${IN_DIR}/protein.pdb"
params="${IN_DIR}/params.json"

if [ ! -f "${pdb_in}" ]; then
  echo "ERROR: Missing /job/input/protein.pdb" | tee -a "${LOG}"
  echo '{"status":"error","message":"missing protein.pdb"}' > "${STATUS}"
  exit 1
fi

preset="cg_popc_50ns"
if [ -f "${params}" ]; then
  preset=$(python3 -c "import json; print(json.load(open('${params}')).get('preset','cg_popc_50ns'))")
fi

echo "Preset: ${preset}" | tee -a "${LOG}"
echo "GROMACS: $(gmx --version | head -n 1)" | tee -a "${LOG}"
echo "INSANE:  $(command -v insane)" | tee -a "${LOG}"
echo "MARTINIZE2: $(command -v martinize2)" | tee -a "${LOG}"

cd "${WORK_DIR}"

# Copy Martini 3 bundle locally so includes resolve from current dir
# (We don't assume exact subfolder names inside /opt/martini)
echo "[setup] Copying Martini bundle into work dir..." | tee -a "${LOG}"
cp -r /opt/martini/* . 2>&1 | tee -a "${LOG}" || true
echo "[setup] Martini work dir contents:" | tee -a "${LOG}"
ls -la | head -n 50 | tee -a "${LOG}"

# 1) Coarse-grain protein
echo "[1/4] martinize2..." | tee -a "${LOG}"
martinize2 -f "${pdb_in}" -o topol.top -x cg.pdb -ff martini3001 \
  2>&1 | tee -a "${LOG}"

# 2) Build membrane + solvent + ions (POPC)
echo "[2/4] insane (membrane build)..." | tee -a "${LOG}"
insane \
  -f cg.pdb \
  -o system.gro \
  -p topol.top \
  -l POPC \
  -sol W \
  -salt 0.15 \
  -box 12,12,16 \
  2>&1 | tee -a "${LOG}"

# 3) GROMACS EM
cat > em.mdp << 'EOF'
integrator  = steep
nsteps      = 2000
emtol       = 500
emstep      = 0.01
cutoff-scheme = Verlet
nstlist     = 20
rlist       = 1.1
coulombtype = Reaction-Field
rcoulomb    = 1.1
vdwtype     = Cut-off
rvdw        = 1.1
pbc         = xyz
EOF

echo "[3/4] grompp+mdrun (EM)..." | tee -a "${LOG}"
gmx grompp -f em.mdp -c system.gro -p topol.top -o em.tpr -maxwarn 10 \
  2>&1 | tee -a "${LOG}"
gmx mdrun -deffnm em 2>&1 | tee -a "${LOG}"

# 4) Very short MD (just to prove it runs)
cat > md.mdp << 'EOF'
integrator    = md
dt            = 0.02
nsteps        = 25000
nstxout-compressed = 100
nstenergy     = 100
nstlog        = 100
cutoff-scheme = Verlet
nstlist       = 20
rlist         = 1.1
coulombtype   = Reaction-Field
rcoulomb      = 1.1
vdwtype       = Cut-off
rvdw          = 1.1
tcoupl        = v-rescale
tc-grps       = System
tau-t         = 1.0
ref-t         = 310
pcoupl        = berendsen
pcoupltype    = isotropic
tau-p         = 12.0
ref-p         = 1.0
compressibility = 3e-4
pbc           = xyz
EOF

echo "[4/4] grompp+mdrun (MD short)..." | tee -a "${LOG}"
gmx grompp -f md.mdp -c em.gro -p topol.top -o md.tpr -maxwarn 10 \
  2>&1 | tee -a "${LOG}"
gmx mdrun -deffnm md 2>&1 | tee -a "${LOG}"

# Export key outputs
cp -f md.xtc md.gro topol.top "${OUT_DIR}/" 2>/dev/null || true
cp -f "${LOG}" "${OUT_DIR}/md.log" 2>/dev/null || true
echo '{"status":"done"}' > "${STATUS}"
touch "${OUT_DIR}/done.marker"
echo "DONE" | tee -a "${LOG}"

