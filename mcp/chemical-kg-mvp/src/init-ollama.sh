#!/bin/bash

# Ollama Model Initialization Script
# This script runs inside the ollama-setup container to download required models

echo "Pulling Ollama models..."

# Get model names from environment variables with defaults
LLM_MODEL=${OLLAMA_MODEL:-llama3.2:latest}
EMBEDDING_MODEL=${EMBEDDING_MODEL:-nomic-embed-text}

# List of models to pull
MODELS=(
    "$LLM_MODEL"
    "$EMBEDDING_MODEL"
)

# Set Ollama host for client commands
export OLLAMA_HOST=ollama:11434

# Wait a bit for Ollama service to be ready
echo "Waiting for Ollama service to start..."
sleep 10

for model in "${MODELS[@]}"; do
    echo "Pulling model: $model"
    # Capture output while filtering verbose messages
    OUTPUT=$(ollama pull "$model" 2>&1)
    PULL_STATUS=$?
    echo "$OUTPUT" | grep -v "pulling manifest" | grep -v "pulling [0-9a-f]" || true
    
    if [ $PULL_STATUS -eq 0 ]; then
        echo "✓ Successfully pulled $model"
    else
        echo "✗ Failed to pull $model"
        # Try some common alternatives
        if [ "$model" = "llama3.2:latest" ]; then
            echo "Trying alternative: mistral"
            OUTPUT=$(ollama pull "mistral" 2>&1)
            echo "$OUTPUT" | grep -v "pulling manifest" | grep -v "pulling [0-9a-f]" || true
        elif [ "$model" = "nomic-embed-text" ]; then
            echo "Trying alternative: all-minilm"
            OUTPUT=$(ollama pull "all-minilm" 2>&1)
            echo "$OUTPUT" | grep -v "pulling manifest" | grep -v "pulling [0-9a-f]" || true
        fi
    fi
done

echo "Model pulling complete!"

# List available models
echo "Available models:"
ollama list