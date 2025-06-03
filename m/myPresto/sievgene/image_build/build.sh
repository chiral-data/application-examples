#!/bin/bash
#

cp ${MYPRESTO_SRC}/sievgene_pack_240118.tar.gz ./
cp ../../../@common/run_dok.sh ./

docker build -t mypresto_sievgene_dok .
rm sievgene_pack_240118.tar.gz
rm run_dok.sh

