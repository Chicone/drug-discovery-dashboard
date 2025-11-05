# Drug Discovery Dashboard

This project is an interactive dashboard developed to support research in **computational drug discovery**, combining both **classical molecular modeling** and **AI-based methods**.  
It provides a unified environment for molecular visualization, descriptor calculation, and property prediction to facilitate compound exploration and hypothesis testing.

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
