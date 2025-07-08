#!/bin/bash
#
# Enhanced Antibody-antigen simulation using GROMACS (Docker version)
# Includes CDR-specific analysis based on 4G6K gevokizumab structure
# Reference: https://academic.oup.com/database/article/doi/10.1093/database/baae015/7631862?login=false#468860647
#            http://www.mdtutorials.com/gmx/lysozyme/index.html
#            https://www.bonvinlab.org/education/HADDOCK3/HADDOCK3-antibody-antigen/
#
# Note: This version uses printf instead of echo -e for Docker compatibility

GMX=/usr/local/gromacs/avx2_256/bin/gmx

# Define CDR residues for 4G6K gevokizumab (from HADDOCK tutorial)
# These correspond to all 6 CDR regions (L1, L2, L3, H1, H2, H3)
CDR_RESIDUES="31,32,33,34,35,52,54,55,56,100,101,102,103,104,105,106,151,152,169,170,173,211,212,213,214,216"

echo "=== GROMACS Antibody-Antigen Simulation with CDR Analysis ==="
echo "CDR residues (4G6K numbering): $CDR_RESIDUES"

## SYSTEM PREPARATION
echo "=== Setting up system ==="

# Create topology
echo 15 | $GMX pdb2gmx -f lightdock_0.pdb -ff charmm27 -water spc -o complex.gro -ignh 

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
echo "=== Running simulations ==="

# Run energy minimization
$GMX mdrun -v -deffnm em
echo 10 0 | $GMX energy -f em.edr -o potential.xvg

# Create enhanced index file with CDR groups
echo "=== Creating index groups including CDR regions ==="
cat > create_index.txt << EOF
r 31,32,33,34,35,52,54,55,56,100,101,102,103,104,105,106,151,152,169,170,173,211,212,213,214,216
name 20 CDR_All
q
EOF

$GMX make_ndx -f em.gro -o index.ndx < create_index.txt
rm create_index.txt

echo "CDR residues selected: 31,32,33,34,35,52,54,55,56,100,101,102,103,104,105,106,151,152,169,170,173,211,212,213,214,216"

echo "Index groups created. CDR_All group is number 20."

# NVT equilibration 
$GMX grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr 
$GMX mdrun -deffnm nvt
echo 16 0 | $GMX energy -f nvt.edr -o temperature.xvg

# NPT equilibration 
$GMX grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
$GMX mdrun -deffnm npt
echo 18 0 | $GMX energy -f npt.edr -o pressure.xvg
echo 24 0 | $GMX energy -f npt.edr -o density.xvg

# Production MD
$GMX grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_10.tpr
$GMX mdrun -deffnm md_0_10

# Recentering and Rewrapping Coordinates
echo "=== Processing trajectories ==="
echo "0" | $GMX trjconv -s md_0_10.tpr -f md_0_10.xtc -o dynamic-nopbc.xtc -pbc mol -center

## ANALYSIS
echo "=== Running analysis ==="

# Standard protein analysis
echo "--- Standard Protein Analysis ---"

# Root Mean Square Deviation (RMSD) - Group 4 (Backbone) for both reference and comparison
printf "4\n4\n" | $GMX rms -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsd_protein.xvg -tu ns

# RMSF - Group 4 (Backbone)
echo "4" | $GMX rmsf -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsf_protein.xvg -res

# SASA - Group 1 (Protein)
printf "1\n1\n" | $GMX sasa -s md_0_10.tpr -f dynamic-nopbc.xtc -o sasa_protein.xvg

# Radius of Gyration - Group 1 (Protein)
echo "1" | $GMX gyrate -s md_0_10.tpr -f dynamic-nopbc.xtc -o gyrate_protein.xvg

# Hydrogen Bonds - Group 1 (Protein) for both groups
printf "1\n1\n" | $GMX hbond -s md_0_10.tpr -f dynamic-nopbc.xtc -num hbond_protein.xvg

# CDR-specific analysis
echo "--- CDR-Specific Analysis ---"

# CDR RMSD analysis (Group 20 is CDR_All)
printf "20\n20\n" | $GMX rms -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsd_cdr.xvg -tu ns

# CDR RMSF analysis 
echo "20" | $GMX rmsf -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsf_cdr.xvg -res

# CDR SASA analysis
printf "20\n20\n" | $GMX sasa -s md_0_10.tpr -f dynamic-nopbc.xtc -o sasa_cdr.xvg

# CDR Radius of Gyration
echo "20" | $GMX gyrate -s md_0_10.tpr -f dynamic-nopbc.xtc -o gyrate_cdr.xvg

# Comparison analysis: Protein vs CDR RMSD in one plot
echo "--- Comparative Analysis ---"

# Create a combined RMSD comparison
printf "4\n20\n" | $GMX rms -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsd_comparison.xvg -tu ns

echo "=== Simulation and analysis complete ==="