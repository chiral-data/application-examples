#!/bin/bash
#
# This is a script to be executed after the job completion
#

echo "post processing"
rm -f *.pdb
rm -f *.gro
rm -f *.itp
rm -f *.top
rm -f *.edr
rm -f *.log
rm -f *.tpr
rm -f *.trr
rm -f *.xvg
rm -f *.cpt
rm -f *.xtc
