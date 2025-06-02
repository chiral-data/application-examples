#!/bin/bash
#
# the latest version of lightdock is released in 2023, so the python version has to be v3.8
# Otherwise there will be backward compatibility issues, such as
#   - since scipy 1.14: keyword argument 'turbo' removed from linalg.eigh() 
#   - since numpy 1.25: numpy.alltrue is deprecated

echo "Downloading receptor and ligand protein examples..."A
wget https://raw.githubusercontent.com/lightdock/lightdock.github.io/master/tutorials/0.9.3/simple_docking/data/2UUY_rec.pdb
wget https://raw.githubusercontent.com/lightdock/lightdock.github.io/master/tutorials/0.9.3/simple_docking/data/2UUY_lig.pdb


echo "Enable the flags to remove OXT (--noxt) atoms, hydrogens (--noh) and waters (--now), and the ANM support"
lightdock3_setup.py 2UUY_rec.pdb 2UUY_lig.pdb --noxt --noh --now -anm

lightdock3.py setup.json 100 -c 1 -l 0

cd swarm_0
lgd_generate_conformations.py ../2UUY_rec.pdb ../2UUY_lig.pdb gso_100.out 200
cd ..

# Create results directory
mkdir -p /opt/artifact/lightdock_results
cp -r ./swarm_0 /opt/artifact/lightdock_results/

ls /opt/artifact/lightdock_results/
