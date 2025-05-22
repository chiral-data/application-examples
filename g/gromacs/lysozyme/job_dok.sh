#!/bin/bash
#
#


echo "Downloading input files ..."
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.3-fix-gromacs/g/gromacs/lysozyme/job.sh
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.3-fix-gromacs/g/gromacs/lysozyme/ions.mdp
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.3-fix-gromacs/g/gromacs/lysozyme/md.mdp
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.3-fix-gromacs/g/gromacs/lysozyme/mdout.mdp
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.3-fix-gromacs/g/gromacs/lysozyme/minim.mdp
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.3-fix-gromacs/g/gromacs/lysozyme/npt.mdp
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.3-fix-gromacs/g/gromacs/lysozyme/nvt.mdp
echo ""

echo "Running the job ..."
sh job.sh

echo "Copy the output files into the DOK artifact"
mkdir -p /opt/artifact/gromacs_result
cp ./potential.xvg /opt/artifact/gromacs_result/
cp ./*.trr /opt/artifact/
ls /opt/artifact/gromacs_result/
