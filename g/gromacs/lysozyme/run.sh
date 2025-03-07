#!/bin/bash

# A Gromacs Tutorial from
# http://www.mdtutorials.com/gmx/lysozyme/01_pdb2gmx.html
# by Prof. Justin A. Lemkul

grep -v HOH 1aki.pdb > 1AKI_clean.pdb
echo 15 | gmx pdb2gmx -f 1AKI_clean.pdb -o 1AKI_processed.gro -water spce
