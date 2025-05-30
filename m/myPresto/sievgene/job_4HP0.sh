#!/bin/bash
#
#

# check the pdb file
../bin/get_pdb_info.pl 4HP0.pdb

# create the input file for "pdbcheck"
echo 4HP0.pdb > inp_pdbcheck
echo 4HP0_1.pdb >> inp_pdbcheck
echo -alt >> inp_pdbcheck
echo -ss >> inp_pdbcheck
echo -disableHet >> inp_pdbcheck
# run "pdbcheck"
../bin/pdbcheck < inp_pdbcheck

# create the input file for "tplgeneX" 
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

# create the input file for "tplgeneL" 
echo 2 > inp_tplgeneL
echo ligand >> inp_tplgeneL
echo 3 >> inp_tplgeneL
echo gaff21.db >> inp_tplgeneL
echo no >> inp_tplgeneL
# run "tplgeneL"
../bin/exec_tplgeneL.sh < inp_tplgeneL

# prepare the input file and run "sievgene"
cp ../input/inp_sievgene1 .
../bin/sievgene_d < inp_sievgene1
