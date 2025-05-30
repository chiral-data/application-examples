#!/bin/bash
#
#

echo "Downloading new workflow file with use_gpu ..."
mkdir -p /app && cd /app
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.3-add-gromacs-py/g/gromacs/protein_ligand/app/workflow.py
mkdir -p /data && cd /data
wget https://github.com/chiral-data/application-examples/blob/v0.2.3-add-gromacs-py/g/gromacs/protein_ligand/data/ions.pdb
wget https://github.com/chiral-data/application-examples/blob/v0.2.3-add-gromacs-py/g/gromacs/protein_ligand/data/ligand.gro
wget https://github.com/chiral-data/application-examples/blob/v0.2.3-add-gromacs-py/g/gromacs/protein_ligand/data/ligand.itp
wget https://github.com/chiral-data/application-examples/blob/v0.2.3-add-gromacs-py/g/gromacs/protein_ligand/data/structure.pdb
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.4/g/gromacs/protein_ligand/workflow.yml
echo ""

cd /workspace
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.4/g/gromacs/protein_ligand/job.sh
echo "Running the job ..."
sh job.sh

echo "Copy the output files into the DOK artifact"
mkdir -p /opt/artifact/gromacs_result
cp -r /data /opt/artifact/gromacs_result/
ls /opt/artifact/gromacs_result/
