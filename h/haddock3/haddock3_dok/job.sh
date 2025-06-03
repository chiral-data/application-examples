#!/bin/bash
#
# From example: https://www.bonvinlab.org/education/HADDOCK3/HADDOCK3-antibody-antigen/


echo "Downloading examples..."
wget wget https://surfdrive.surf.nl/files/index.php/s/MuC8YogNPj9Ac31/download -O HADDOCK3-antibody-antigen.zip
unzip HADDOCK3-antibody-antigen.zip
cd HADDOCK3-antibody-antigen/runs
rm -rf */
mkdir run1
cd ..

# Install PyMol
echo "Installing PyMol..."
#pip install pymol-open-source

echo "Preparing the antibody structure..."
pdb_fetch 4G6K | pdb_tidy -strict | pdb_selchain -H | pdb_delhetatm | pdb_fixinsert | pdb_selaltloc | pdb_keepcoord | pdb_selres -1:120 | pdb_tidy -strict > 4G6K_H.pdb
pdb_fetch 4G6K | pdb_tidy -strict | pdb_selchain -L | pdb_delhetatm | pdb_fixinsert | pdb_selaltloc | pdb_keepcoord | pdb_selres -1:107 | pdb_tidy -strict > 4G6K_L.pdb
pdb_merge 4G6K_H.pdb 4G6K_L.pdb | pdb_reres -1 | pdb_chain -A | pdb_chainxseg | pdb_tidy -strict > 4G6K_clean.pdb

echo "Preparing the antigen structure..."
pdb_fetch 4I1B | pdb_tidy -strict | pdb_delhetatm | pdb_selaltloc | pdb_keepcoord | pdb_chain -B | pdb_chainxseg | pdb_tidy -strict > 4I1B_clean.pdb

echo "Running HADDOCK3..."
# Change the number of cores
sed -i 's/ncores = 50/ncores = 2/g' workflows/docking-antibody-antigen.cfg
haddock3 workflows/docking-antibody-antigen.cfg 

# Create results directory
mkdir -p /opt/artifact/haddock3_results
cp -r . /opt/artifact/haddock3_results/

ls /opt/artifact/haddock3_results/
