#!/bin/bash

# A Gromacs Tutorial from
# http://www.mdtutorials.com/gmx/lysozyme/01_pdb2gmx.html
# by Prof. Justin A. Lemkul

GMX=/usr/local/gromacs/avx2_256/bin/gmx

wget https://files.rcsb.org/download/1AKI.pdb
grep -v HOH 1AKI.pdb > 1AKI_clean.pdb
echo 15 | $GMX pdb2gmx -f 1AKI_clean.pdb -o 1AKI_processed.gro -water spce
$GMX editconf -f 1AKI_processed.gro -o 1AKI_newbox.gro -c -d 1.0 -bt cubic
$GMX solvate -cp 1AKI_newbox.gro -cs spc216.gro -o 1AKI_solv.gro -p topol.top
$GMX grompp -f ions.mdp -c 1AKI_solv.gro -p topol.top -o ions.tpr
echo 13 | $GMX genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral
$GMX grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr
$GMX mdrun -deffnm em
echo 10 0 | $GMX energy -f em.edr -o potential.xvg
