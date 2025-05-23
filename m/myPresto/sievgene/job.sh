#!/bin/bash
#
#


echo "Downloading input files ..."
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.3-add-myPresto/m/myPresto/sievgene/job.sh
echo ""

echo "Running the job ..."
sh job.sh

echo "Copy the output files into the DOK artifact"
mkdir -p /opt/artifact/results
cp -r ./tmp_4HP0 /opt/artifact/results/
ls /opt/artifact/results/

