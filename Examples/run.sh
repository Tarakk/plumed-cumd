gmx_mpi_d grompp -f nvt.mdp -c confout.gro -p topol.top 
gmx_mpi_d mdrun -plumed plumed.dat
