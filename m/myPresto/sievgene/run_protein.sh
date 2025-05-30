#!/bin/bash
#
#

PROTEIN="4HP0"

mkdir tmp_$PROTEIN
cd tmp_$PROTEIN
cp ../sample/$PROTEIN.pdb . 

echo "Downloading the launching script ..."
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.3-add-myPresto/m/myPresto/sievgene/job_$PROTEIN.sh
echo "Running the sievgene job for $PROTEIN..."
sh job_$PROTEIN.sh
echo ""

echo "Copy the output files into the DOK artifact"
mkdir -p /opt/artifact/results
cd ..
cp -r ./tmp_$PROTEIN /opt/artifact/results/
ls /opt/artifact/results/

