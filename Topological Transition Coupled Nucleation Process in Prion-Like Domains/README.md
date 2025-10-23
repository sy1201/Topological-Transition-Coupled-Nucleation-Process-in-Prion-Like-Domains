#############################################################
### Supporting Code for the Paper: FUS Protein Simulations ###
#############################################################

This repository contains all input files necessary to reproduce the simulations reported in the paper. The simulations were performed using the 
**LAMMPS molecular dynamics package (version 23 Jun 2023)**, running on Linux with standard libraries.

-------------------------------------------------------------
### 1. Overview
-------------------------------------------------------------

The simulations investigate the thermodynamic and dynamic behavior of FUS protein condensates using a coarse-grained molecular dynamics model.  
Interactions between coarse-grained beads are modeled using a soft Lennard-Jones potential and harmonic bonded interactions.

All simulations were performed in a cubic box under periodic boundary conditions,with a Langevin thermostat maintaining the temperature at 300 K.

-------------------------------------------------------------
### 2. Directory Structure
-------------------------------------------------------------

This archive contains several directories, each corresponding to a different simulation setup. Each folder includes its own **readme.txt** file describing
the purpose and contents of the simulation.

1-equilibration/ # Initial equilibration runs
2-production/ # Long production simulations
3-analysis/ # Post-processing and data analysis scripts
4-fus.data # Example initial data file

-------------------------------------------------------------
### 3. LAMMPS Compilation Instructions
-------------------------------------------------------------

Instructions on how to build LAMMPS are available at: https://docs.lammps.org/Build_make.html

A brief summary is provided below for convenience:

1. In a clean working directory, download the desired version of LAMMPS:  https://www.lammps.org/download.html

2. Compile the required packages:  
    make yes-MOLECULE
    make yes-USER-SOFT
    make yes-KSPACE
    make yes-RCB
These enable molecular topology, soft LJ potentials, and domain balancing.

3. Build the executable: make mpi or make serial

4. Copy the resulting executable file (lmp_mpi or lmp_serial) to the directory containing your simulation input files.

5. Run the simulation: ./lmp_mpi -in in.fus

-------------------------------------------------------------
### 4. Reproducibility
-------------------------------------------------------------
All simulations can be repeated by running: mpirun -np 16 lmp_mpi -in in.fus
Ensure that the fus.data file is placed in the same directory as in.fus. The results should reproduce the thermodynamic and structural data reported in the manuscript.