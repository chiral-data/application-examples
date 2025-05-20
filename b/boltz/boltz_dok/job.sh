#!/bin/bash
#
# Boltz-1 Democratizing biomolecular interaction modeling
#

# Install required packages
pip install boltz -U

# Run 
python3 -m boltz.main predict prot_no_msa.yaml

# Create results directory
mkdir -p /opt/artifact/boltz_results/prot_no_msa/
cp prot_no_msa_model_0.cif /opt/artifact/boltz_results/prot_no_msa/
cp confidence_prot_no_msa_model_0.json /opt/artifact/boltz_results/prot_no_msa/
#cp pae_prot_no_msa_model_0.npz /opt/artifact/boltz_results/prot_no_msa/
#cp pde_prot_no_msa_model_0.npz /opt/artifact/boltz_results/prot_no_msa/
cp plddt_prot_no_msa_model_0.npz /opt/artifact/boltz_results/prot_no_msa/
cp prot_no_msa.npz /opt/artifact/boltz_results/prot_no_msa/
#cp prot_no_msa_model_diffusion_samples-1.cif /opt/artifact/boltz_results/prot_no_msa/

echo "Prediction completed. Results are available in results folder."
ls /opt/artifact/boltz_results/prot_no_msa/
