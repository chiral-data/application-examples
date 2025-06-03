#!/bin/bash
#

# Put "cosgene_pack_240724.tar.gz" under same directory as this script
cp ../../../../@common/run_dok.sh ./

docker build -t mypresto_cosgene_dok .
rm cosgene_pack_240724.tar.gz
rm run_dok.sh

