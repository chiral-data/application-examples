#!/bin/bash
#
# https://huggingface.co/collections/google/txgemma-release-67dd92e931c857d15e4d1e87

# Output file for exported text
OUTPUT_FILE="txgemma_output.txt"

echo "Downloading the Python script ..." | tee -a "$OUTPUT_FILE"
wget https://raw.githubusercontent.com/chiral-data/application-examples/refs/heads/v0.2.4-TxGemma/t/TxGemma/TxGemma_dok/example.py -O example.py 2>&1 | tee -a "$OUTPUT_FILE"

# Run 
echo "Run TxGemma ..." | tee -a "$OUTPUT_FILE"
python3 -m example 2>&1 | tee -a "$OUTPUT_FILE"

echo "" | tee -a "$OUTPUT_FILE"
echo "Output saved to: $OUTPUT_FILE" | tee -a "$OUTPUT_FILE"