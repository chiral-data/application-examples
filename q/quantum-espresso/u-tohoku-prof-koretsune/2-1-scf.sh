#!/bin/bash
#

mkdir 2-1-scf
cd 2-1-scf
wget http://www.cmpt.phys.tohoku.ac.jp/~koretsune/SATL_qe_tutorial/files/graphene_band/graphene.scf.in
wget http://www.cmpt.phys.tohoku.ac.jp/~koretsune/SATL_qe_tutorial/files/C.pz-van_ak.UPF
docker run -it --rm --gpus all --runtime=nvidia --ipc=host -v ${PWD}:/host_pwd -w /host_pwd nvcr.io/hpc/quantum_espresso:qe-7.3.1 pw.x -input graphene.scf.in > graphene.scf.out
# singularity run --nv -B ${PWD}:/host_pwd --pwd /host_pwd docker://nvcr.io/hpc/quantum_espresso:qe-7.3.1 pw.x -input graphene.scf.in > graphene.scf.out
grep -a "^\!" graphene.scf.out
