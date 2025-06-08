#!/bin/bash
#
# Boltz-1 Democratizing biomolecular interaction modeling
#

echo "Downloading 4G6K.fasta ..."
#wget https://www.rcsb.org/fasta/entry/4G6K -O 4G6K.fasta
#sed -n '1,2p' 4G6K.fasta > 4G6K_rec.fasta
#sed -n '3,4p' 4G6K.fasta > 4G6K_lig.fasta
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.4_boltz2lightdock/b/boltz/boltz_dok/4G6K_rec.fasta -O 4G6K_rec.fasta
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.4_boltz2lightdock/b/boltz/boltz_dok/4G6K_lig.fasta -O 4G6K_lig.fasta

# Run 
echo "Run boltz calculation ..."
python3 -m boltz.main predict 4G6K_rec.fasta --use_msa_server
python3 -m boltz.main predict 4G6K_lig.fasta --use_msa_server

# Create results directory
mkdir -p /opt/artifact/boltz_results_4G6K
echo "Results for 4G6K_rec:"
ls ./boltz_results_4G6K_rec/predictions/4G6K_rec/
echo "Results for 4G6K_lig:"
ls ./boltz_results_4G6K_lig/predictions/4G6K_lig/
cp ./boltz_results_4G6K_rec/predictions/4G6K_rec/4G6K_rec_model_0.pdb /opt/artifact/boltz_results_4G6K/
cp ./boltz_results_4G6K_rec/predictions/4G6K_rec/confidence_4G6K_rec_model_0.json /opt/artifact/boltz_results_4G6K/
cp ./boltz_results_4G6K_rec/predictions/4G6K_rec/plddt_4G6K_rec_model_0.npz /opt/artifact/boltz_results_4G6K/
cp ./boltz_results_4G6K_lig/predictions/4G6K_lig/4G6K_lig_model_0.pdb /opt/artifact/boltz_results_4G6K/
cp ./boltz_results_4G6K_lig/predictions/4G6K_lig/confidence_4G6K_lig_model_0.json /opt/artifact/boltz_results_4G6K/
cp ./boltz_results_4G6K_lig/predictions/4G6K_lig/plddt_4G6K_lig_model_0.npz /opt/artifact/boltz_results_4G6K/

echo "Prediction completed. Results are available in results folder."
ls /opt/artifact/boltz_results_4G6K/
