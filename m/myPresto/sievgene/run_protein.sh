#!/bin/bash
#
#

$PROTEIN=4HP0

echo "Downloading input files ..."
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.3-add-myPresto/m/myPresto/sievgene/job_$PROTEIN.sh
echo ""

mkdir tmp_$PROTEIN
cd tmp_$PROTEIN
cp ../sample/$PROTEIN.pdb . 

echo "Running the sievgene job for $PROTEIN..."
sh job_$PROTEIN.sh

echo "Copy the output files into the DOK artifact"
mkdir -p /opt/artifact/results
cp -r ./tmp_$PROTEIN /opt/artifact/results/
ls /opt/artifact/results/

