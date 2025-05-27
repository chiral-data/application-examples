#!/bin/bash
#
#

# conda run --no-capture-output -n complex_setup_env python /app/workflow.py --config /data/workflow.yml
conda run --no-capture-output -n biobb_wf_protein-complex_md_setup python /app/workflow.py --config /data/workflow.yml
