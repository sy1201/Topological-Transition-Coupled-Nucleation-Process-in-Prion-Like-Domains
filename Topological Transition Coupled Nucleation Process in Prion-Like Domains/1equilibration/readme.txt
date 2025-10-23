Place the file **Anneal.in** and **fus.data** in the same directory, and run the simulation in LAMMPS using the input script: **mpirun -np 16 lmp_mpi -in Anneal.in**.

The main purpose of this input file is to ensure that the initial structure of the FUS solution does not form clusters at the beginning of the simulation.

It serves as a pre-equilibration stage to generate the file **restart.fus.100000000**, which will be used as the starting point for the subsequent equilibrium simulation.
