#!/bin/bash
#

cp ../../../@common/run_dok.sh ./
docker build -t mypresto_sievgene_dok .
rm run_dok.sh
