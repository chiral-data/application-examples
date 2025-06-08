#!/bin/bash
#
# Boltz-1 Democratizing biomolecular interaction modeling
#

echo "Downloading prot_no_msa.yaml ..."
wget https://raw.githubusercontent.com/jwohlwend/boltz/refs/heads/main/examples/prot_no_msa.yaml

# Run 
echo "Run boltz calculation ..."
python3 -m boltz.main predict prot_no_msa.yaml

# Create results directory
mkdir -p /opt/artifact/boltz_results_prot_no_msa
cp ./boltz_results_prot_no_msa/predictions/prot_no_msa/prot_no_msa_model_0.cif /opt/artifact/boltz_results_prot_no_msa/
cp ./boltz_results_prot_no_msa/predictions/prot_no_msa/confidence_prot_no_msa_model_0.json /opt/artifact/boltz_results_prot_no_msa/
cp ./boltz_results_prot_no_msa/predictions/prot_no_msa/plddt_prot_no_msa_model_0.npz /opt/artifact/boltz_results_prot_no_msa/

echo "Prediction completed. Results are available in results folder."
ls /opt/artifact/boltz_results_prot_no_msa/