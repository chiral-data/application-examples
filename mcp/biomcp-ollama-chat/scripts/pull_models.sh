#!/bin/bash

# Pull multiple models for Ollama
echo "Pulling Ollama models..."

# List of models to pull
MODELS=(
    "gemma3:4b"
    "llama3.2:latest"
    "deepseek-r1:14b"
    "qwen3:8b"
)

for model in "${MODELS[@]}"; do
    echo "Pulling model: $model"
    ollama pull "$model"
    if [ $? -eq 0 ]; then
        echo "Successfully pulled $model"
    else
        echo "Failed to pull $model"
    fi
done

echo "Model pulling complete!"

# Start Ollama server
exec ollama serve