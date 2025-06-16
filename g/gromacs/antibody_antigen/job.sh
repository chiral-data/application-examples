#!/bin/bash
#

GMX=/usr/local/gromacs/avx2_256/bin/gmx

# Create topology
echo 15 | $GMX pdb2gmx -f 5dk3-rs1-311_5k-100ns.pdb -ff charmm27 -water spc -o complex.gro 

## SET UP SIMULATION
# Define simulation box
gmx editconf -f complex.gro -o newbox.gro -bt dodecahedron -d 1.0

# Solvate simulation box
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro 

# Assemble atomic level description of system (.tpr file) 
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr 

# Add chloride atoms
echo 13 | $GMX genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

# Assemble .tpr file again (this time for input to mdrun)
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr


## RUN SIMULATION
# Run energy minimization
gmx mdrun -v -deffnm em
echo 10 0 | $GMX energy -f em.edr -o potential.xvg

# NVT equilibration 
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr 
gmx mdrun -deffnm nvt
echo 16 0 | $GMX energy -f nvt.edr -o temperature.xvg

# NPT equilibration 
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
gmx mdrun -deffnm npt
echo 18 0 | $GMX energy -f npt.edr -o pressure.xvg
echo 24 0 | $GMX energy -f npt.edr -o density.xvg

# Production MD
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_10.tpr
gmx mdrun -deffnm md_0_10

# Recentering and Rewrapping Coordinates
gmx trjconv -s md_0_10.tpr -f md_0_10.xtc -o dynamic-nopbc.xtc -pbc cluster


## ANALYSIS
# Root Mean Square Deviation (RMSD)
echo 4 4 | $GMX rms -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsd.xvg -tu ns
# RMSF
echo 4 4 | $GMX rmsf -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsf.xvg -res
# SASA
echo 4 4 | $GMX sasa -s md_0_10.tpr -f dynamic-nopbc.xtc -o sasa.xvg
# Radius of Gyration
echo 1 | $GMX gyrate -s md_0_10.tpr -f dynamic-nopbc.xtc -o gyrate.xvg
# Hydrogen Bonds
echo 4 4 | $GMX hbond -s md_0_10.tpr -f dynamic-nopbc.xtc -o hbond.xvg
