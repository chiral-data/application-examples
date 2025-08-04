#!/bin/bash
#
#

echo "--- Individual Arguments ---"
echo "Running script: $0 ..." 
echo "Job script to be downloaded: $1"
# Extract filename without query parameters
filename=$(basename "$1" | cut -d'?' -f1)
echo "Script to run: $filename"
echo ""

# --- Checking the number of arguments ---
echo "--- Number of Arguments ---"
echo "Total number of arguments received: $#"
if [ "$#" -eq 0 ]; then
  echo "No arguments provided. Please provide some arguments when running the script."
  echo "Usage: ./run.sh <url_for_script>  ..."
  exit
elif [ "$#" -gt 1 ]; then
  echo "You provided 2 or more arguments."
  exit
fi
echo ""

echo "Downloading the job script ..."
wget -O "$filename" "$1"
echo ""

echo "Running the job script ..."
sh $filename
echo ""



