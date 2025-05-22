#!/bin/bash

# A Gromacs Tutorial from
# http://www.mdtutorials.com/gmx/lysozyme/01_pdb2gmx.html
# by Prof. Justin A. Lemkul

GMX=/usr/local/gromacs/avx2_256/bin/gmx

wget https://files.rcsb.org/download/1AKI.pdb
grep -v HOH 1AKI.pdb > 1AKI_clean.pdb

# Step One: Prepare the Topology
# Step Two: Examine the Topology
echo 15 | $GMX pdb2gmx -f 1AKI_clean.pdb -o 1AKI_processed.gro -water spce

# Step Three: Defining the Unit Cell & Adding Solvent
$GMX editconf -f 1AKI_processed.gro -o 1AKI_newbox.gro -c -d 1.0 -bt cubic
$GMX solvate -cp 1AKI_newbox.gro -cs spc216.gro -o 1AKI_solv.gro -p topol.top

# Step Four: Adding Ions
$GMX grompp -f ions.mdp -c 1AKI_solv.gro -p topol.top -o ions.tpr
echo 13 | $GMX genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral

# Step Five: Energy Minimization
$GMX grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr
$GMX mdrun -deffnm em
echo 10 0 | $GMX energy -f em.edr -o potential.xvg

# Step Six: Equilibration
$GMX grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
$GMX mdrun -deffnm nvt
echo 16 0 | $GMX energy -f nvt.edr -o temperature.xvg

# Step Seven: Equilibration, Part 2
$GMX grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
$GMX mdrun -deffnm npt
echo 18 0 | $GMX energy -f npt.edr -o pressure.xvg
echo 24 0 | $GMX energy -f npt.edr -o density.xvg

# Step Eight: Production MD
$GMX grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
$GMX mdrun -deffnm md_0_1
# check what have been produced
ls -lh

# Step Nine: Analysis
echo 1 0 | $GMX trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center
echo 4 4 | $GMX rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns
echo 4 4 | $GMX rms -s em.tpr -f md_0_1_noPBC.xtc -o rmsd_xtal.xvg -tu ns
echo 1 | $GMX gyrate -s md_0_1.tpr -f md_0_1_noPBC.xtc -o gyrate.xvg
