#!/bin/bash
#
# the latest version of lightdock is released in 2023, so the python version has to be v3.8
# Otherwise there will be backward compatibility issues, such as
#   - since scipy 1.14: keyword argument 'turbo' removed from linalg.eigh() 
#   - since numpy 1.25: numpy.alltrue is deprecated
# 
# From example: https://lightdock.org/tutorials/0.9.3/simple_docking


#echo "Copying receptor and ligand protein results..."
#cp /home/ubuntu/chiral/application-examples/b/boltz/boltz_dok/boltz_results_4G6K/4G6K_rec_model_0.pdb ./4G6K_rec.pdb
#cp /home/ubuntu/chiral/application-examples/b/boltz/boltz_dok/boltz_results_4G6K/4G6K_lig_model_0.pdb ./4G6K_lig.pdb

echo "Enable the flags to remove OXT (--noxt) atoms, hydrogens (--noh) and waters (--now), and the ANM support"
lightdock3_setup.py 4G6K_rec_model_0.pdb 4G6K_lig_model_0.pdb --noxt --noh --now -anm

lightdock3.py setup.json 100 -c 1 -l 0

cd swarm_0
# Generate conformations
lgd_generate_conformations.py ../4G6K_rec_model_0.pdb ../4G6K_lig_model_0.pdb gso_100.out 20

# Clustering
lgd_cluster_bsas.py gso_100.out
cd ..

echo "Docking completed. Results are available in results folder."