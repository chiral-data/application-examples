#!/bin/bash
#
# Boltz-1 Democratizing biomolecular interaction modeling
#

# Run 
echo "Run boltz calculation ..."
python3 -m boltz.main predict Antibody-1-heavy_chain.fasta --use_msa_server --output_format pdb
python3 -m boltz.main predict Antibody-1-light_chain.fasta --use_msa_server --output_format pdb
python3 -m boltz.main predict Antibody-2-heavy_chain.fasta --use_msa_server --output_format pdb
python3 -m boltz.main predict Antibody-2-light_chain.fasta --use_msa_server --output_format pdb
python3 -m boltz.main predict Antibody-3-heavy_chain.fasta --use_msa_server --output_format pdb
python3 -m boltz.main predict Antibody-3-light_chain.fasta --use_msa_server --output_format pdb

# Create results directory
mkdir -p /opt/artifact/boltz_results
echo "Results for Antibody-1:"
ls ./boltz_results_Antibody-1-heavy_chain/predictions/Antibody-1-heavy_chain/
ls ./boltz_results_Antibody-1-light_chain/predictions/Antibody-1-light_chain/
echo "Results for Antibody-2:"
ls ./boltz_results_Antibody-2-heavy_chain/predictions/Antibody-2-heavy_chain/
ls ./boltz_results_Antibody-2-light_chain/predictions/Antibody-2-light_chain/
echo "Results for Antibody-3:"
ls ./boltz_results_Antibody-3-heavy_chain/predictions/Antibody-3-heavy_chain/
ls ./boltz_results_Antibody-3-light_chain/predictions/Antibody-3-light_chain/

cp ./boltz_results_Antibody-1-heavy_chain/predictions/Antibody-1-heavy_chain/Antibody-1-heavy_chain_model_0.pdb /opt/artifact/boltz_results/
cp ./boltz_results_Antibody-1-light_chain/predictions/Antibody-1-light_chain/Antibody-1-light_chain_model_0.pdb /opt/artifact/boltz_results/
cp ./boltz_results_Antibody-2-heavy_chain/predictions/Antibody-2-heavy_chain/Antibody-2-heavy_chain_model_0.pdb /opt/artifact/boltz_results/
cp ./boltz_results_Antibody-2-light_chain/predictions/Antibody-2-light_chain/Antibody-2-light_chain_model_0.pdb /opt/artifact/boltz_results/
cp ./boltz_results_Antibody-3-heavy_chain/predictions/Antibody-3-heavy_chain/Antibody-3-heavy_chain_model_0.pdb /opt/artifact/boltz_results/
cp ./boltz_results_Antibody-3-light_chain/predictions/Antibody-3-light_chain/Antibody-3-light_chain_model_0.pdb /opt/artifact/boltz_results/

echo "Prediction completed. Results are available in results folder."
ls /opt/artifact/boltz_results/
