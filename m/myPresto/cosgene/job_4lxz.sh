#!/bin/bash
#

# 1. Prepare the PDB files
# Create a working directory for the job
mkdir work_MD_4lxz
cd work_MD_4lxz
# Cope sample file to the working directory
cp ../sample/pdb4lxz.ent .
# Check the PDB file
perl ../bin/get_pdb_info.pl pdb4lxz.ent

# 2. Process the PDB file
# Extract Chains A
perl ../bin/select_chain.pl A pdb4lxz.ent > tmp_selectChain
# Delete HOH
perl ../bin/del_res.pl HOH tmp_selectChain > tmp_delRes1
# Delete PG4
perl ../bin/del_res.pl PG4 tmp_delRes1 > tmp_delRes2
# Delete SHH
perl ../bin/del_res.pl SHH tmp_delRes2 > tmp_delRes3
# Confirm the result
perl ../bin/get_pdb_info.pl tmp_delRes3
# Extract restraint only
perl ../bin/select_res.pl SHH tmp_selectChain > tmp_ligand

# 3. Execute pdbcheck
# Run with protein
echo tmp_delRes3 > inp_pdbcheck
echo tmp_pdb_checked.pdb >> inp_pdbcheck
echo -alt >> inp_pdbcheck
cat inp_pdbcheck
../bin/pdbcheck < inp_pdbcheck

# Run with ligand
echo tmp_ligand > inp_pdbcheck_lig
echo lig.pdb >> inp_pdbcheck_lig
echo -alt >> inp_pdbcheck_lig
cat inp_pdbcheck_lig
../bin/pdbcheck < inp_pdbcheck_lig

# 4. Execute Hgene
../bin/Hgene -ipdb lig.pdb -p -mop AM1BCC -omol2 lig.mol2

# 5. Execute tplgeneL
echo 2 > inp_tplgeneL
echo lig.mol2 >> inp_tplgeneL
echo 1 >> inp_tplgeneL
echo gaff21.db >> inp_tplgeneL
echo no >> inp_tplgeneL
cat inp_tplgeneL
../bin/exec_tplgeneL.sh < inp_tplgeneL

# 6. Execute tylgeneX
echo pdb > inp_tplgeneX
echo C99 >> inp_tplgeneX
echo tmp_pdb_checked.pdb >> inp_tplgeneX
echo pdb >> inp_tplgeneX
echo Pro_0.pdb >> inp_tplgeneX
echo Pro_0.tpl >> inp_tplgeneX
../bin/exec_tplgeneX.sh < inp_tplgeneX

# 7. Merge the PDF files of the protein and ligand
cat Pro_0.pdb lig_tplL.pdb > Pro_1.pdb

# 8. Execute setwater
echo Pro_1.pdb > inp_setwater
echo N >> inp_setwater
echo wat.pdb >> inp_setwater
echo S >> inp_setwater
echo 8.0 >> inp_setwater
echo C >> inp_setwater
echo 1.0 >> inp_setwater
echo 1.0 >> inp_setwater
echo 3 >> inp_setwater
echo Y >> inp_setwater
echo setwater.log >> inp_setwater
echo N >> inp_setwater
cat inp_setwater
../bin/setwater < inp_setwater
cat Pro_1.pdb wat.pdb > Pro_2.pdb

# 9. Execute add_ion
echo Pro_2.pdb > inp_add_ion
echo Pro_3.pdb >> inp_add_ion
echo 4 >> inp_add_ion
echo Y >> inp_add_ion
echo 6.0 >> inp_add_ion
echo 0.2 >> inp_add_ion
cat inp_add_ion
../bin/add_ion < inp_add_ion

# 10. Execute tplgeneX again
echo pdb > inp_tplgeneX2
echo C99 >> inp_tplgeneX2
echo Pro_3.pdb >> inp_tplgeneX2
echo pdb >> inp_tplgeneX2
echo Pro_4.pdb >> inp_tplgeneX2
echo Pro_4.tpl >> inp_tplgeneX2
echo tip3p >> inp_tplgeneX2
cat inp_tplgeneX2
../bin/exec_tplgeneX.sh < inp_tplgeneX2
head -20 Pro_4.tpl

# 11. Execute tpl2capbc
echo Pro_4.tpl >> inp_tpl2capbc
echo wat.pdb >> inp_tpl2capbc
echo capbc >> inp_tpl2capbc
echo M_all.res >> inp_tpl2capbc
cat inp_tpl2capbc
../bin/tpl2capbc < inp_tpl2capbc

# 12. Prepare min.inp
cp ../sample/min.inp .

# 13. Execute cosgene (Energy minimization)
../bin/cosgene < min.inp > min.out &

# 14. Execute SHAKEinp
echo Pro_4.tpl > inp_shake
echo Pro_4_min.pdb >> inp_shake
echo shape_inp >> inp_shake
echo yes >> inp_shake
cat inp_shake
../bin/SHAKEinp < inp_shake

# 14. Prepare md.inp
cp ../sample/md.inp .

# 15. Execute cosgene (Molecular dynamics)
../bin/cosgene < md.inp > md.out &
grep LOOP md.out

# 16. Create and download results directory
mkdir -p /opt/artifact/cosgene_results
cp -r . /opt/artifact/cosgene_results/

ls /opt/artifact/cosgene_results/