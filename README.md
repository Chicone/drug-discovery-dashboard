# Drug Discovery Dashboard

This project is an interactive dashboard developed to support research in **computational drug discovery**, combining both **classical molecular modeling** and **AI-based methods**.  
It provides a unified environment for molecular visualization, descriptor calculation, and property prediction to facilitate compound exploration and hypothesis testing.

## 🚀 Installation

1. Clone the repository  
git clone https://github.com/Chicone/drug-discovery-dashboard.git  
cd drug-discovery-dashboard

2. Set up the backend (FastAPI)  

It is recommended to use conda or venv for environment management.
Using conda:  
conda create -n drugdash python=3.10  
conda activate drugdash  
conda env create -f environment.yml

Or using pip directly:  
python -m venv venv  
source venv/bin/activate  (On Windows: venv\Scripts\activate)  
pip install -r requirements.txt  

Then start the FastAPI server:  
uvicorn app.main:app --reload  

The backend will be available at:  
http://127.0.0.1:8000

3. Set up the frontend (React)  
Open a new terminal and go to the frontend directory:  
cd frontend  
npm install  
npm start  

The frontend will be available at:  
http://localhost:3000

4. Access the dashboard  
Once both servers are running, open your browser and visit:  
http://localhost:3000

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


