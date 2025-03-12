#!/bin/bash

cd ausurf
docker run -it --rm --gpus all --runtime=nvidia --ipc=host -v ${PWD}:/host_pwd -w /host_pwd nvcr.io/hpc/quantum_espresso:qe-7.3.1 ./run_qe.sh
# singularity run --nv -B${PWD}:/host_pwd --pwd /host_pwd docker://nvcr.io/hpc/quantum_espresso:qe-7.3.1 ./run_qe.sh
