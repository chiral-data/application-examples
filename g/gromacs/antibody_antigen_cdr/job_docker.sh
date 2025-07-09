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

# Create script to generate CDR index if it doesn't exist
cat > create_cdr_index_v2.sh << 'EOF'
#!/bin/bash
GMX=/usr/local/gromacs/avx2_256/bin/gmx

echo "=== Creating Enhanced Index File with Individual CDR Groups ==="

# Backup existing index file
if [ -f "index.ndx" ]; then
    cp index.ndx index_original.ndx
    echo "Backed up existing index.ndx to index_original.ndx"
fi

# Create commands for gmx make_ndx
cat > make_ndx_commands.txt << 'EOFINNER'
q
EOFINNER

# First, create a standard index file
echo "Creating standard index file..."
$GMX make_ndx -f em.gro -o index_temp.ndx < make_ndx_commands.txt

# Now add CDR groups to the standard index
cat > make_ndx_cdr_commands.txt << 'EOFINNER'
r 31 32 33 34 35
name 18 CDR_L1
r 52 54 55 56
name 19 CDR_L2
r 100 101 102 103 104 105 106
name 20 CDR_L3
r 151 152
name 21 CDR_H1
r 169 170 173
name 22 CDR_H2
r 211 212 213 214 216
name 23 CDR_H3
r 31 32 33 34 35 52 54 55 56 100 101 102 103 104 105 106 151 152 169 170 173 211 212 213 214 216
name 24 CDR_All
q
EOFINNER

# Add CDR groups to the index
echo "Adding CDR groups..."
$GMX make_ndx -f em.gro -n index_temp.ndx -o index_cdr.ndx < make_ndx_cdr_commands.txt

# Clean up
rm -f make_ndx_commands.txt make_ndx_cdr_commands.txt index_temp.ndx

if [ -f "index_cdr.ndx" ]; then
    echo "Successfully created index_cdr.ndx"
else
    echo "ERROR: Failed to create index_cdr.ndx"
    exit 1
fi
EOF

chmod +x create_cdr_index_v2.sh
./create_cdr_index_v2.sh
rm create_cdr_index_v2.sh

echo "=== Using index file with CDR groups ==="
echo "CDR groups in index_cdr.ndx:"
echo "  Group 17: CDR_L1 (residues 31-35)"
echo "  Group 18: CDR_L2 (residues 52,54-56)"
echo "  Group 19: CDR_L3 (residues 100-106)"
echo "  Group 20: CDR_H1 (residues 151-152)"
echo "  Group 21: CDR_H2 (residues 169-170,173)"
echo "  Group 22: CDR_H3 (residues 211-214,216)"
echo "  Group 23: CDR_All (all CDR residues)"

# NVT equilibration 
$GMX grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index_cdr.ndx -o nvt.tpr 
$GMX mdrun -deffnm nvt
echo 16 0 | $GMX energy -f nvt.edr -o temperature.xvg

# NPT equilibration 
$GMX grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index_cdr.ndx -o npt.tpr
$GMX mdrun -deffnm npt
echo 18 0 | $GMX energy -f npt.edr -o pressure.xvg
echo 24 0 | $GMX energy -f npt.edr -o density.xvg

# Production MD
$GMX grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index_cdr.ndx -o md_0_10.tpr
$GMX mdrun -deffnm md_0_10

# Recentering and Rewrapping Coordinates
echo "=== Processing trajectories ==="
printf "0\n0\n" | $GMX trjconv -s md_0_10.tpr -f md_0_10.xtc -o dynamic-nopbc.xtc -pbc mol -center

## ANALYSIS
echo "=== Running analysis ==="

# Standard protein analysis
echo "--- Standard Protein Analysis ---"

# Root Mean Square Deviation (RMSD) - Group 4 (Backbone) for both reference and comparison
printf "4\n4\n" | $GMX rms -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsd_protein.xvg -n index_cdr.ndx -tu ns

# RMSF - Group 4 (Backbone)
echo "4" | $GMX rmsf -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsf_protein.xvg -n index_cdr.ndx -res

# SASA - Group 1 (Protein)
printf "1\n1\n" | $GMX sasa -s md_0_10.tpr -f dynamic-nopbc.xtc -o sasa_protein.xvg -n index_cdr.ndx

# Radius of Gyration - Group 1 (Protein)
echo "1" | $GMX gyrate -s md_0_10.tpr -f dynamic-nopbc.xtc -o gyrate_protein.xvg -n index_cdr.ndx

# Hydrogen Bonds - Group 1 (Protein) for both groups
printf "1\n1\n" | $GMX hbond -s md_0_10.tpr -f dynamic-nopbc.xtc -num hbond_protein.xvg -n index_cdr.ndx

# CDR-specific analysis
echo "--- CDR-Specific Analysis ---"

# Analyze each CDR individually (using portable shell syntax)
# CDR_L1 - Group 17
echo "=== Analyzing CDR_L1 (Group 17) ==="
printf "17\n17\n" | $GMX rms -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsd_CDR_L1.xvg -n index_cdr.ndx -tu ns
echo "17" | $GMX rmsf -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsf_CDR_L1.xvg -n index_cdr.ndx -res
printf "17\n17\n" | $GMX sasa -s md_0_10.tpr -f dynamic-nopbc.xtc -o sasa_CDR_L1.xvg -n index_cdr.ndx
echo "17" | $GMX gyrate -s md_0_10.tpr -f dynamic-nopbc.xtc -o gyrate_CDR_L1.xvg -n index_cdr.ndx

# CDR_L2 - Group 18
echo "=== Analyzing CDR_L2 (Group 18) ==="
printf "18\n18\n" | $GMX rms -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsd_CDR_L2.xvg -n index_cdr.ndx -tu ns
echo "18" | $GMX rmsf -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsf_CDR_L2.xvg -n index_cdr.ndx -res
printf "18\n18\n" | $GMX sasa -s md_0_10.tpr -f dynamic-nopbc.xtc -o sasa_CDR_L2.xvg -n index_cdr.ndx
echo "18" | $GMX gyrate -s md_0_10.tpr -f dynamic-nopbc.xtc -o gyrate_CDR_L2.xvg -n index_cdr.ndx

# CDR_L3 - Group 19
echo "=== Analyzing CDR_L3 (Group 19) ==="
printf "19\n19\n" | $GMX rms -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsd_CDR_L3.xvg -n index_cdr.ndx -tu ns
echo "19" | $GMX rmsf -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsf_CDR_L3.xvg -n index_cdr.ndx -res
printf "19\n19\n" | $GMX sasa -s md_0_10.tpr -f dynamic-nopbc.xtc -o sasa_CDR_L3.xvg -n index_cdr.ndx
echo "19" | $GMX gyrate -s md_0_10.tpr -f dynamic-nopbc.xtc -o gyrate_CDR_L3.xvg -n index_cdr.ndx

# CDR_H1 - Group 20
echo "=== Analyzing CDR_H1 (Group 20) ==="
printf "20\n20\n" | $GMX rms -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsd_CDR_H1.xvg -n index_cdr.ndx -tu ns
echo "20" | $GMX rmsf -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsf_CDR_H1.xvg -n index_cdr.ndx -res
printf "20\n20\n" | $GMX sasa -s md_0_10.tpr -f dynamic-nopbc.xtc -o sasa_CDR_H1.xvg -n index_cdr.ndx
echo "20" | $GMX gyrate -s md_0_10.tpr -f dynamic-nopbc.xtc -o gyrate_CDR_H1.xvg -n index_cdr.ndx

# CDR_H2 - Group 21
echo "=== Analyzing CDR_H2 (Group 21) ==="
printf "21\n21\n" | $GMX rms -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsd_CDR_H2.xvg -n index_cdr.ndx -tu ns
echo "21" | $GMX rmsf -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsf_CDR_H2.xvg -n index_cdr.ndx -res
printf "21\n21\n" | $GMX sasa -s md_0_10.tpr -f dynamic-nopbc.xtc -o sasa_CDR_H2.xvg -n index_cdr.ndx
echo "21" | $GMX gyrate -s md_0_10.tpr -f dynamic-nopbc.xtc -o gyrate_CDR_H2.xvg -n index_cdr.ndx

# CDR_H3 - Group 22
echo "=== Analyzing CDR_H3 (Group 22) ==="
printf "22\n22\n" | $GMX rms -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsd_CDR_H3.xvg -n index_cdr.ndx -tu ns
echo "22" | $GMX rmsf -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsf_CDR_H3.xvg -n index_cdr.ndx -res
printf "22\n22\n" | $GMX sasa -s md_0_10.tpr -f dynamic-nopbc.xtc -o sasa_CDR_H3.xvg -n index_cdr.ndx
echo "22" | $GMX gyrate -s md_0_10.tpr -f dynamic-nopbc.xtc -o gyrate_CDR_H3.xvg -n index_cdr.ndx

# CDR_All - Group 23
echo "=== Analyzing CDR_All (Group 23) ==="
printf "23\n23\n" | $GMX rms -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsd_CDR_All.xvg -n index_cdr.ndx -tu ns
echo "23" | $GMX rmsf -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsf_CDR_All.xvg -n index_cdr.ndx -res
printf "23\n23\n" | $GMX sasa -s md_0_10.tpr -f dynamic-nopbc.xtc -o sasa_CDR_All.xvg -n index_cdr.ndx
echo "23" | $GMX gyrate -s md_0_10.tpr -f dynamic-nopbc.xtc -o gyrate_CDR_All.xvg -n index_cdr.ndx

# Comparison analysis: Protein backbone vs all CDRs
echo "--- Comparative Analysis ---"

# Create a combined RMSD comparison (Backbone vs All CDRs)
printf "4\n23\n" | $GMX rms -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsd_comparison.xvg -n index_cdr.ndx -tu ns

# Create a multi-CDR RMSD comparison (all individual CDRs in one plot)
echo "--- Multi-CDR RMSD Comparison ---"
printf "17\n18\n19\n20\n21\n22\n" | $GMX rms -s md_0_10.tpr -f dynamic-nopbc.xtc -o rmsd_all_cdrs.xvg -n index_cdr.ndx -tu ns -ng 6

echo "=== Simulation and analysis complete ==="

# List all generated files
echo "=== Generated analysis files ==="
echo "Trajectory file:"
ls -la dynamic-nopbc.xtc 2>/dev/null
echo ""
echo "Protein-wide analysis:"
ls -la *protein*.xvg 2>/dev/null | sort
echo ""
echo "Individual CDR analysis:"
ls -la *CDR*.xvg 2>/dev/null | sort
echo ""
echo "Comparison files:"
ls -la *comparison*.xvg *all_cdrs*.xvg 2>/dev/null | sort