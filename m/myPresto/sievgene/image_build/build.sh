#!/bin/bash
#

# Put "sievgene_pack_240118.tar.gz" under same directory as this script
cp ../../../../@common/run_dok.sh ./

docker build -t mypresto_sievgene_dok .
rm sievgene_pack_240118.tar.gz
rm run_dok.sh

