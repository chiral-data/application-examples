#!/bin/bash
# This is a script to be executed before job starts


echo "preprocessing"
mkdir ausurf
cd ausurf
wget https://repository.prace-ri.eu/git/UEABS/ueabs/-/raw/master/quantum_espresso/test_cases/small/Au.pbe-nd-van.UPF
wget https://repository.prace-ri.eu/git/UEABS/ueabs/-/raw/master/quantum_espresso/test_cases/small/ausurf.in
wget https://gitlab.com/NVHPC/ngc-examples/-/raw/master/qe/single-node/run_qe.sh
chmod +x run_qe.sh
