#!/bin/bash
#
#

echo "Downloading input files ..."
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.4_gromacsAb/g/gromacs/antibody_antigen/job.sh
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.4_gromacsAb/g/gromacs/antibody_antigen/ions.mdp
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.4_gromacsAb/g/gromacs/antibody_antigen/md.mdp
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.4_gromacsAb/g/gromacs/antibody_antigen/mdout.mdp
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.4_gromacsAb/g/gromacs/antibody_antigen/em.mdp
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.4_gromacsAb/g/gromacs/antibody_antigen/npt.mdp
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.4_gromacsAb/g/gromacs/antibody_antigen/nvt.mdp
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.4_gromacsAb/g/gromacs/antibody_antigen/5dk3-rs1-311_5k-100ns.pdb
echo ""

echo "Running the job ..."
sh job.sh

echo "Copy the output files into the DOK artifact"
mkdir -p /opt/artifact/gromacs_result
cp -r . /opt/artifact/gromacs_result/
ls /opt/artifact/gromacs_result/
