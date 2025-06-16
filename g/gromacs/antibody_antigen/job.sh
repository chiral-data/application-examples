#!/bin/bash
#
# Antibody-antigen simulation using GROMACS
# Reference: https://academic.oup.com/database/article/doi/10.1093/database/baae015/7631862?login=false#468860647
#            http://www.mdtutorials.com/gmx/lysozyme/index.html

GMX=/usr/local/gromacs/avx2_256/bin/gmx

# Create topology
echo 15 | $GMX pdb2gmx -f 5dk3-rs1-311_5k-100ns.pdb -ff charmm27 -water spc -o complex.gro -ignh 

## SET UP SIMULATION
# Define simulation box
$GMX editconf -f complex.gro -o newbox.gro -bt dodecahedron -d 1.0

# Solvate simulation box
$GMX solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro 

# Assemble atomic level description of system (.tpr file) 
$GMX grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr 

# Add chloride atoms
echo 13 | $GMX genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

# Assemble .tpr file again (this time for input to mdrun)
$GMX grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr


## RUN SIMULATION
# Run energy minimization
$GMX mdrun -v -deffnm em
echo 10 0 | $GMX energy -f em.edr -o potential.xvg

# Create index file for system groups
echo q  < /dev/null |  $GMX make_ndx -f em.gro -o index.ndx

# NVT equilibration 
$GMX grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr 
$GMX mdrun -deffnm nvt
echo 16 0 | $GMX energy -f nvt.edr -o temperature.xvg

# NPT equilibration 
$GMX grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
$GMX mdrun -deffnm npt
echo 18 0 | $GMX energy -f npt.edr -o pressure.xvg
echo 24 0 | $GMX energy -f npt.edr -o density.xvg

# Production MD
$GMX grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_10.tpr
$GMX mdrun -deffnm md_0_10

# Recentering and Rewrapping Coordinates
$GMX trjconv -s md_0_10.tpr -f md_0_10.xtc -o dynamic-nopbc.xtc -pbc cluster


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
