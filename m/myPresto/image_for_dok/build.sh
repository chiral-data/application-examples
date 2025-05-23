#!/bin/bash
#

cp ../../../@common/run.sh ./
docker build -t mypresto_sievgene_dok .
rm run.sh
