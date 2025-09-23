In Silico Discovery of NDM-1 Metallo-β-Lactamase Inhibitors for Klebsiella pneumoniae

 Overview

Antimicrobial resistance represents a critical global health challenge. Among resistant pathogens, Klebsiella pneumoniae has emerged as a significant concern due to its ability to produce New Delhi Metallo-β-Lactamase (NDM-1), which hydrolyzes a wide range of β-lactam antibiotics.
This project applies in silico approaches to identify novel inhibitors of NDM-1 through molecular docking, ADMET profiling, and molecular dynamics (MD) simulations.

 Objectives

Retrieve and prepare the NDM-1 target protein (PDB: 6NY7)

Screen standard ligands and sulfonamide derivatives for inhibitory potential

Perform molecular docking studies using AutoDock Vina

Predict pharmacokinetic and toxicity profiles through ADMET analysis

Evaluate the dynamic stability of ligand–protein complexes using 100 ns MD simulations in GROMACS

🛠️ Methodology
Protein and Ligand Preparation

Target protein retrieved from the Protein Data Bank (PDB ID: 6NY7)

Standard ligands: Penicillin, Ampicillin, Benzyl Penicillin, Sulfonamide, Mercaptopurine, etc.

Sulfonamide derivative dataset (~4500 compounds) obtained from PubChem and filtered using Lipinski’s Rule of Five and Veber’s criteria

Molecular Docking

Protein preparation and docking conducted using AutoDock Vina via PyRx

Virtual screening yielded ~60 promising sulfonamide derivatives with strong binding affinities

ADMET Analysis

Tools: ADMETLab 2.0 and Protox-II

Evaluated Absorption, Distribution, Metabolism, Excretion, and Toxicity properties

Two lead candidates were identified:

(3Z)-N-hydroxypenta-1,3-diene-2-sulfonamide (PubChem CID: 118156306)

N-hydroxyfuran-2-sulfonamide (PubChem CID: 46175386)

Molecular Dynamics Simulations

Conducted using GROMACS 5.0.4 for 100 ns

Structural analyses performed: Root Mean Square Deviation (RMSD), Root Mean Square Fluctuation (RMSF), Radius of Gyration (Rg), and Solvent Accessible Surface Area (SASA)

Results confirmed the stability of the selected ligand–protein complexes

 Results
Compound	PubChem CID	Docking Score (kcal/mol)	ADMET Profile	MD Stability
(3Z)-N-hydroxypenta-1,3-diene-2-sulfonamide	118156306	[insert value]	Favorable	Stable
N-hydroxyfuran-2-sulfonamide	46175386	[insert value]	Favorable	Stable
Sulfonamide (control)	–	[insert value]	Moderate	Stable

Figures to include:

Representative docking poses

Interaction diagrams

RMSD, RMSF, Rg, and SASA plots

(Store all figures in /figures and link them here)

 Workflow
graph TD
A[Protein Retrieval (PDB 6NY7)] --> B[Ligand Collection (PubChem)]
B --> C[Molecular Docking with AutoDock Vina]
C --> D[ADMET Screening]
D --> E[Lead Candidate Selection]
E --> F[MD Simulations (GROMACS)]
F --> G[Stable Inhibitor Identification]

 Conclusion

This in silico study identified two sulfonamide derivatives—(3Z)-N-hydroxypenta-1,3-diene-2-sulfonamide and N-hydroxyfuran-2-sulfonamide—as promising inhibitors of NDM-1. Both compounds demonstrated strong binding affinities, favorable ADMET profiles, and stable interactions during MD simulations. These findings highlight their potential as lead compounds for the development of novel therapeutics against NDM-1–mediated resistance, pending experimental validation.

Repository Structure
InSilico-NDM1-Inhibitors/
├── data/              # Protein and ligand datasets
├── docking/           # Docking configurations and results
├── admet/             # ADMET prediction reports
├── md_simulations/    # GROMACS input/output files
├── figures/           # Docking images and MD plots
├── scripts/           # Python/R/shell scripts
├── README.md          # Project documentation
└── requirements.txt   # Dependencies


Reproducibility Instructions

Clone this repository:

git clone https://github.com/your-username/InSilico-NDM1-Inhibitors.git
cd InSilico-NDM1-Inhibitors


Install required software:

AutoDock Vina / PyRx

GROMACS ≥ 5.0

Python ≥ 3.8 with RDKit, OpenBabel

Run docking:

vina --config config.txt


Run ADMET predictions (example scripts in /scripts)

Perform MD simulations:

gmx grompp -f md.mdp -c complex.gro -p topol.top -o md.tpr
gmx mdrun -deffnm md

 License

This project is distributed under the MIT License.

Contact

For questions, collaborations, or queries regarding this project, please feel free to reach out:

Name: Prashanth E

Email: prashantheprashanth584@gmail.com
