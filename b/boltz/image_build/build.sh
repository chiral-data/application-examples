#!/bin/bash
#

# --- Checking the number of arguments ---
echo "--- Number of Arguments ---"
echo "Total number of arguments received: $#"
if [ "$#" -eq 0 ]; then
  echo "No arguments provided. Please provide some arguments when running the script."
  echo "Usage: ./build.sh nvidia  ..."
  exit
elif [ "$#" -gt 1 ]; then
  echo "You provided 2 or more arguments."
  exit
fi

echo "Building container image with base image from $1"
echo ""

mkdir -p $1 && cd $1
cp ../../../../@common/run_dok.sh ./
cp ../Dockerfile_$1 ./Dockerfile
docker build -t boltz_dok_$1 --platform linux/amd64 .
cd ..
rm -rf $1
