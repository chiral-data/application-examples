#!/bin/bash
#
#


echo "Downloading input files ..."

wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.3-add-gromacs-py/g/gromacs/protein_ligand/job.sh

mkdir app
mkdir data

cd app
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.3-add-gromacs-py/g/gromacs/protein_ligand/app/workflow.py

cd ../data
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.3-add-gromacs-py/g/gromacs/protein_ligand/data/ions.pdb
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.3-add-gromacs-py/g/gromacs/protein_ligand/data/ligand.gro
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.3-add-gromacs-py/g/gromacs/protein_ligand/data/ligand.itp
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.3-add-gromacs-py/g/gromacs/protein_ligand/data/structure.pdb
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.3-add-gromacs-py/g/gromacs/protein_ligand/data/workflow.yml
cd ..

echo ""

echo "Running the job ..."
sh job.sh

echo "Copy the output files into the DOK artifact"
mkdir -p /opt/artifact/gromacs_result
cp ./*.* /opt/artifact/gromacs_result/
ls /opt/artifact/gromacs_result/
