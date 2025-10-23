Place the files **relax.in** and **restart.fus.100000000** in the same directory, and run the simulation in LAMMPS using the input script: **mpirun -np 16 lmp_mpi -in relax.in**.

The main purpose of this input file is to allow the FUS protein solution to spontaneously form condensates within the cubic simulation box.

This run generates the files **dump.relax.lammpstrj** and **log.lammps** as outputs.
