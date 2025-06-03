#!/bin/bash
#

cp ${MYPRESTO_SRC}/sievgene_pack_240118.tar.gz ./
sudo docker build -t mypresto_sievgene .
rm sievgene_pack_240118.tar.gz
