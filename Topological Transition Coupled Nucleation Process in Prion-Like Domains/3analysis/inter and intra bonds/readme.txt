## Purpose ##
This Python script analyzes **LAMMPS trajectory files (.lammpstrj)** to extract statistical information of **FUS protein systems**, including cluster identification, radius of gyration, and inter-chain bonding.  
Different command-line options allow the user to select specific analyses, such as:  
- Molecular radius of gyration (-rg)  
- Cluster composition and structure (-c)  
- Largest cluster information (-mc)  
- Inter-bond statistics of the largest cluster (-mcb)  
- Cluster size distribution (-cv)  
- Centering of the simulation box or the largest cluster (-bm, -cm)  

The program generates corresponding output files for further analysis and visualization.

## Usage ##

Run in terminal:

python analysis.py <lammpstrj_file> [options]

Example:

python analysis.py dump.relax.lammpstrj -sy -mc -mcb


This reads the file `dump.relax.lammpstrj` using the `sy` mode, outputs the largest cluster and its internal bonding information.


##  Command-Line Arguments ##

| Short | Long | Description |
|--------|------|-------------|
| lammpstrj| — | Input `.lammpstrj` file (required) |
| -bm | --box_mid     | Shift box center to (0,0,0) |
| -cm | --cluster_mid | Shift the center of mass of the largest cluster to (0,0,0)` |
| -yzs| --cluster_yzs | Use yzs-mode for cluster identification (choose one mode only) |
| -sy | --cluster_sy  | Use sy-mode for cluster identification (choose one mode only) |
| -rg | --mol_rg      | Output radius of gyration for each molecule |
| -c  | --cluster     | Output cluster composition information |
| -mc | --max_cluster | Output the largest cluster information |
| -mcb| --max_interbonds | Output the internal bond distribution of the largest cluster |
| -cv | --cluster_volum  | Output cluster size statistics |
| -o  | --output      | Write a new `.lammpstrj` file |
| -r  | --read        | Read frames from start to end with given step (e.g. `-r 0 1000 10`) |

---

##  Output Files ##

Depending on selected options, the program generates the following files:

| Filename | Description |
|-----------|--------------|
| cluster.txt | Composition of clusters at each timestep |
| drop.txt    | Average radius and count of clusters (generated with `-c`) |
| max_cluster.txt    | Chain indices of the largest cluster at each timestep |
| max_interbonds.txt | Bond counts and radial distribution of the largest cluster |
| cluster_volum.txt  | Cluster size distribution (1–9 chains and >10 chains) |
| <input>_rg.txt     | Radius of gyration for each molecule |
| bond_number.txt    | Intra- and inter-cluster bond counts (only in `-sy` mode) |
| Output `.lammpstrj` file | If `-o` is used, a shifted or modified trajectory is written |


The script also depends on the following local modules (must be in the same directory):
- dumpreader.py    — for reading LAMMPS `.lammpstrj` files  
- math_function.py — for computing centers of mass, radii, and cluster analysis  
