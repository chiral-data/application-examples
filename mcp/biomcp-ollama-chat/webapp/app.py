import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from flask import Flask, render_template, request, jsonify, Response, stream_with_context
import asyncio
import json
import logging
from dotenv import load_dotenv
import aiohttp

from mcp_client.ollama_integration import OllamaBioMCPIntegration

# Load environment variables
load_dotenv()

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = Flask(__name__)

# Configuration
OLLAMA_HOST = os.getenv('OLLAMA_HOST', 'http://localhost:11434')
DEFAULT_MODEL = os.getenv('OLLAMA_MODEL', 'gemma3:4b')
BIOMCP_VERSION = os.getenv('BIOMCP_VERSION', '0.1.1')

# Available models (can be overridden by actual Ollama models)
AVAILABLE_MODELS = [
    "gemma3:4b",
    "llama3.2:latest",
    "deepseek-r1:14b",
    "qwen3:8b"
]

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/models', methods=['GET'])
def get_models():
    """Get available Ollama models"""
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    
    async def fetch_models():
        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(f"{OLLAMA_HOST}/api/tags") as response:
                    if response.status == 200:
                        data = await response.json()
                        models = [model['name'] for model in data.get('models', [])]
                        # If no models found, return our default list
                        if not models:
                            models = AVAILABLE_MODELS
                        return {"models": models, "default": DEFAULT_MODEL}
                    else:
                        # Return default models if Ollama is not ready
                        return {"models": AVAILABLE_MODELS, "default": DEFAULT_MODEL}
        except Exception as e:
            logger.error(f"Error fetching models: {e}")
            return {"models": AVAILABLE_MODELS, "default": DEFAULT_MODEL}
    
    result = loop.run_until_complete(fetch_models())
    loop.close()
    
    return jsonify(result)

@app.route('/chat', methods=['POST'])
def chat():
    data = request.json
    message = data.get('message', '')
    model = data.get('model', DEFAULT_MODEL)
    
    def generate():
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        
        async def process():
            async with OllamaBioMCPIntegration(
                OLLAMA_HOST, 
                model, 
                BIOMCP_VERSION
            ) as integration:
                async for chunk in integration.process_message(message):
                    yield f"data: {json.dumps(chunk)}\n\n"
        
        # Run the async generator
        for item in loop.run_until_complete(collect_async_gen(process())):
            yield item
        
        loop.close()
    
    return Response(stream_with_context(generate()), content_type='text/event-stream')

async def collect_async_gen(async_gen):
    """Helper to collect async generator results"""
    results = []
    async for item in async_gen:
        results.append(item)
    return results

@app.route('/tools', methods=['GET'])
def get_tools():
    """Get available BioMCP tools"""
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    
    async def fetch_tools():
        async with OllamaBioMCPIntegration(
            OLLAMA_HOST, 
            DEFAULT_MODEL, 
            BIOMCP_VERSION
        ) as integration:
            return await integration.biomcp_client.get_available_tools()
    
    tools = loop.run_until_complete(fetch_tools())
    loop.close()
    
    return jsonify(tools)

@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint"""
    return jsonify({
        "status": "healthy",
        "ollama_host": OLLAMA_HOST,
        "default_model": DEFAULT_MODEL,
        "biomcp_version": BIOMCP_VERSION
    })


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=False)