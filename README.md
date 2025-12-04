# Drug Discovery Dashboard

This project is an interactive dashboard developed to support research in **computational drug discovery**, combining both **classical molecular modeling** and **AI-based methods**.  
It provides a unified environment for molecular visualization, descriptor calculation, and property prediction to facilitate compound exploration and hypothesis testing.

## 🚀 Installation

1. Clone the repository  
git clone https://github.com/Chicone/drug-discovery-dashboard.git  
cd drug-discovery-dashboard

2. Set up the backend (FastAPI)  

### 🔧 Environment Setup

This project uses different environment files depending on the operating system.

> Recommended environment manager: **conda**  
> Recommended Python version: **3.11**

---

#### 🟢 Ubuntu (Linux)

```bash
conda create -n drugdash python=3.11
conda activate drugdash
conda env create -f environment_ubuntu.yml
```

#### 🟢 macOS Intel (x86_64)
```bash
conda create -n drugdash python=3.11
conda activate drugdash
conda env create -f environment_mac_intel.yml
```

#### 🟢 macOS ARM (M1 / M2 / M3 / M4) (not tested yet)
```bash
conda config --set subdir osx-arm64
conda create -n drugdash python=3.11
conda activate drugdash
conda env create -f environment_mac_arm.yml
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

## Features

- 🧬 **3D Molecular Visualization** — render structures directly from SMILES using 3Dmol.js.  
- ⚗️ **Descriptor Panel** — compute classical molecular properties such as molecular weight, logP, tPSA, and rotatable bonds.  
- 🧠 **Predictive Modeling** — integrate both traditional QSAR and modern machine-learning approaches (e.g., regression models, GNNs, diffusion-based predictors).  
- 🧩 **Modular Architecture** — FastAPI backend and React frontend designed for easy extension and deployment.  
- 📊 **Dataset Management** — load and inspect molecular datasets (SMILES, SDF, CSV) for screening and comparison.

## Goals

This dashboard aims to:
- Provide an interactive workspace for exploring chemical libraries.  
- Combine classical cheminformatics with AI-driven property prediction.  
- Support visualization and interpretation of molecular descriptors and model outputs.  
- Serve as a foundation for future extensions (docking, molecular dynamics, ADMET prediction, etc.).

## Tech Stack

- **Backend:** FastAPI (Python)  
- **Frontend:** React + 3Dmol.js  
- **Styling:** Bulma CSS / Tailwind  
- **Data:** Molecular descriptors, SMILES, SDF, CSV  

## Author

Developed by **Luis G. Camara, PhD**  
Computational Chemist & Data Scientist  
University of Geneva — Pharmaceutical Biochemistry Group


