#!/bin/bash
#
#

mkdir tmp_4HP0 
cd tmp_4HP0 
cp ../sample/4HP0.pdb . 

# check the pdb file
../bin/get_pdb_info.pl 4HP0.pdb

# create the input file for "pdbcheck"
echo 4HP0.pdb > inp_pdbcheck
echo 4HP0_1.pdb >> inp_pdbcheck
echo "–alt" >> inp_pdbcheck
echo "–ss" >> inp_pdbcheck
echo "–disableHet" >> inp_pdbcheck
# run "pdbcheck"
../bin/pdbcheck < inp_pdbcheck

# create the input file for "tplgene" 
echo 1 > inp_tplgeneX
echo FF_set1 >> inp_tplgeneX
echo 4HP0_1.pdb >> inp_tplgeneX
echo 1 >> inp_tplgeneX
echo Pro.pdb >> inp_tplgeneX
echo Pro.tpl >> inp_tplgeneX
# run "tplgeneX"
../bin/exec_tplgeneX.sh < inp_tplgeneX 

# create the target probe file 
grep "^HETATM" 4HP0.pdb | grep -v "HOH" > point.pdb

# copy the ligand file
cp ../sample/c001-1.mol2 ./ligand.mol2

# prepare the input file and run "sievgene"
cp ../input/inp_sievgene1 .
../bin/sievgene_d < inp_sievgene1
