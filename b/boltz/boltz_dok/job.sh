#!/bin/bash
#
# Boltz-1 Democratizing biomolecular interaction modeling
#

echo "Downloading 4G6K.fasta ..."
wget https://www.rcsb.org/fasta/entry/4G6K -O 4G6K.fasta
sed -n '1,2p' 4G6K.fasta > 4G6K_rec.fasta
sed -n '3,4p' 4G6K.fasta > 4G6K_lig.fasta

# Run 
echo "Run boltz calculation ..."
python3 -m boltz.main predict 4G6K_rec.fasta
python3 -m boltz.main predict 4G6K_lig.fasta

# Create results directory
mkdir -p /opt/artifact/boltz_results_4G6K
cp ./boltz_results_4G6K/predictions/4G6K_rec/4G6K_rec_model_0.pdb /opt/artifact/boltz_results_4G6K/
cp ./boltz_results_4G6K/predictions/4G6K_rec/confidence_4G6K_rec_model_0.json /opt/artifact/boltz_results_4G6K/
cp ./boltz_results_4G6K/predictions/4G6K_rec/plddt_4G6K_rec_model_0.npz /opt/artifact/boltz_results_4G6K/
cp ./boltz_results_4G6K/predictions/4G6K_lig/4G6K_lig_model_0.pdb /opt/artifact/boltz_results_4G6K/
cp ./boltz_results_4G6K/predictions/4G6K_lig/confidence_4G6K_lig_model_0.json /opt/artifact/boltz_results_4G6K/
cp ./boltz_results_4G6K/predictions/4G6K_lig/plddt_4G6K_lig_model_0.npz /opt/artifact/boltz_results_4G6K/

echo "Prediction completed. Results are available in results folder."
ls /opt/artifact/boltz_results_4G6K/
