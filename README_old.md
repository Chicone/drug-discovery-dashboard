# Drug Discovery Dashboard

## Computational Drug Discovery Platform

This project is an interactive platform for computational drug discovery, combining molecular modeling, simulation, cheminformatics, and future AI-driven workflows within a unified environment.

It is designed to support practical research workflows rather than isolated tools, allowing users to move from molecular structures to simulation, analysis, and decision-making in a single interface.

---

## Current Capabilities

### Molecular Workspace

- 3D visualization of proteins and ligands
- render structures from SMILES
- inspect receptor-ligand systems
- browser-based interface

### Cheminformatics

- molecular weight
- logP
- tPSA
- rotatable bonds
- descriptor calculation
- compound loading and inspection

### Molecular Dynamics

- Martini coarse-grained workflows
- membrane and water systems
- protein-only and protein-ligand simulations
- automated build / equilibration / production runs
- job management and logs

### Multiscale Modeling

- selection of promising CG frames
- CG → AA backmapping
- ligand-aware atomistic reconstruction
- atomistic equilibration workflows

### Analysis

- RMSD
- center-of-mass metrics
- orientation metrics
- activation-related receptor metrics
- comparative run analysis

---

### Example Workflow

```text
Protein structure
→ ligand setup
→ Martini simulation
→ identify relevant frame
→ backmap to atomistic detail
→ equilibrate AA system
→ analyze metrics
```

## 🚀 Installation

1. Clone the repository
```bash  
git clone https://github.com/Chicone/drug-discovery-dashboard.git  
cd drug-discovery-dashboard
git checkout dev

2. Set up the backend (FastAPI)  

### 🔧 Environment Setup

This project uses different environment files depending on the operating system.

> Recommended environment manager: **conda**  
> Recommended Python version: **3.11**

---

#### 🟢 Ubuntu (Linux)

```bash
conda env create -f environment_ubuntu.yml
conda activate drugdash

```

#### 🟢 macOS Intel (x86_64)
```bash
conda create -n drugdash python=3.11
conda activate drugdash
conda env update -f environment_mac_intel.yml
```

#### 🟢 macOS ARM (M1 / M2 / M3 / M4) (not tested yet)
```bash
conda config --set subdir osx-arm64
conda create -n drugdash python=3.11
conda activate drugdash
conda env update -f environment_mac_arm.yml
```

#### 🟢 Windows (not tested yet)

[//]: # (Or using pip directly:  )
[//]: # (python -m venv venv  )
[//]: # (source venv/bin/activate  &#40;On Windows: venv\Scripts\activate&#41;  )
[//]: # (pip install -r requirements.txt  )

### 🚀 Running the backend and frontend
1. Move into the backend directory and start the FastAPI server: 
```bash
cd backend/
uvicorn main:app --reload
```
The backend will be available at:  
http://127.0.0.1:8000

2. Move into the frontend directory and set up the frontend (React)  
Open a new terminal and do:  
```bash
cd frontend  
conda activate drugdash
npm install  
npm run dev
```
The frontend will be available at:  
http://localhost:5173

3. Access the dashboard  
Once both servers are running, open your browser and visit:  
http://localhost:5173

You should see the interactive drug discovery dashboard interface.

---


## Generating Ligand Files for AA Simulations / Backmapping

Use this workflow when a ligand requires atomistic (AA) files for:

- CG → AA backmapping
- AA equilibration
- production MD
- ligand insertion into receptor systems

---

### Canonical Files Required by the AA Pipeline

Each ligand must end up with these standardized files:

```text
<LIGAND>.pdb
<LIGAND>.itp
<LIGAND>.prm
```

Example:

```text
ZM241385.pdb
ZM241385.itp
ZM241385.prm
```

These are the files used directly by the AA workflow.

---

### Step 1: Open CHARMM-GUI

Go to:

```text
https://www.charmm-gui.org/
```

Use the module:

```text
Ligand Reader & Modeler
```

---

### Step 2: Enter the Ligand SMILES

Paste the ligand SMILES string.

Example for `ZM241385`:

```text
C1=COC(=C1)C2=NN3C(=NC(=NC3=N2)NCCC4=CC=C(C=C4)O)N
```

Then continue.

---

### Step 3: Check Hydrogens / Protonation

If prompted, inspect the generated ligand.

Verify:

- no unrealistic formal charges
- OH groups present where expected
- sensible protonation state
- chemically reasonable structure

Then continue.

---

### Step 4: Generate CGenFF Topology

Select:

```text
Make CGenFF topology
```

Choose a short residue name (3 to 6 characters).

Examples:

```text
ZM24
NECA
CPD5
```

Use the same residue name consistently later.

---

### Step 5: Download the Results Archive

Download the generated `.tgz` archive.

Typical useful contents:

```text
ligandrm.pdb
ligandrm.str
toppar.str
gromacs/<RES>.itp
<res>/<res>.prm
```

Example:

```text
ligandrm.pdb
ligandrm.str
toppar.str
gromacs/ZM24.itp
zm24/zm24.prm
```

---

### Step 6: Create the Project Ligand Folder

Create:

```text
experiments/A2AR/ligands/ZM241385/
```

Recommended final structure (after renaming files):

```text
experiments/A2AR/ligands/ZM241385/
    ZM241385.pdb
    ZM241385.itp
    ZM241385.prm
    download.tgz
    ligandrm.pdb
    ligandrm.str
    toppar.str
```

All files remain in the ligand folder root.

---

### Step 7: Standardize File Names

Copy / rename:

```text
ligandrm.pdb      -> ZM241385.pdb
gromacs/ZM24.itp  -> ZM241385.itp
zm24/zm24.prm     -> ZM241385.prm
```

Keep the original downloaded files in the same folder for provenance.

---

### Step 8: Validate Before Use

Before running simulations, verify:

- atom names in `.pdb` match `.itp`
- atom count matches topology
- geometry is sensible
- expected total charge
- no corrupted coordinates

---

### Step 9: Use in the AA Pipeline

Use these canonical files:

```text
ZM241385.pdb
ZM241385.itp
ZM241385.prm
```

for:

- ligand insertion
- CG → AA backmapping
- AA equilibration
- production MD

---

### Notes

- CHARMM-GUI often already provides GROMACS-ready `.itp` files.
- Manual conversion scripts may not be necessary.
- Keep naming consistent across `.pdb`, `.itp`, `.prm`, and topology includes.
- Do not delete original downloaded files.
- Repeat this process for any new ligand.


## Tech Stack

- **Backend:** FastAPI (Python)  
- **Frontend:** React + 3Dmol.js  
- **Styling:** Bulma CSS / Tailwind  
- **Data:** Molecular descriptors, SMILES, SDF, CSV  

## Author

Developed by **Luis G. Camara, PhD**  
Computational Chemist & Data Scientist  
University of Geneva — Pharmaceutical Biochemistry Group


