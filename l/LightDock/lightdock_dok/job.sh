#!/bin/bash
#

echo "Downloading receptor and ligand protein examples..."
wget https://raw.githubusercontent.com/lightdock/lightdock.github.io/master/tutorials/0.9.3/simple_docking/data/2UUY_rec.pdb
wget https://raw.githubusercontent.com/lightdock/lightdock.github.io/master/tutorials/0.9.3/simple_docking/data/2UUY_lig.pdb


echo "Enable the flags to remove OXT (--noxt) atoms, hydrogens (--noh) and waters (--now), and the ANM support"
python3 -m lightdock3_setup.py 2UUY_rec.pdb 2UUY_lig.pdb --noxt --noh --now -anm

python3 -m lightdock3.py setup.json 100 -c 1 -l 0

cd swarm_0
python3 -m lgd_generate_conformations.py ../2UUY_rec.pdb ../2UUY_lig.pdb gso_100.out 200

# Create results directory
mkdir -p /opt/artifact/lightdock_results
cp -r ./swarm_0 /opt/artifact/lightdock_results/

ls /opt/artifact/lightdock_results/
