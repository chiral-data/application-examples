#!/bin/bash
#
# Antibody-antigen simulation using GROMACS
# Reference: https://academic.oup.com/database/article/doi/10.1093/database/baae015/7631862?login=false#468860647
#            http://www.mdtutorials.com/gmx/lysozyme/index.html

GMX=/usr/local/gromacs/avx2_256/bin/gmx

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
