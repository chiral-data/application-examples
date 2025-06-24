import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from flask import Flask, render_template, request, jsonify, Response, stream_with_context
import asyncio
import json
import logging
from dotenv import load_dotenv

from mcp_client.ollama_integration import OllamaBioMCPIntegration

# Load environment variables
load_dotenv()

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = Flask(__name__)

# Configuration
OLLAMA_HOST = os.getenv('OLLAMA_HOST', 'http://localhost:11434')
OLLAMA_MODEL = os.getenv('OLLAMA_MODEL', 'llama3.2:latest')
BIOMCP_VERSION = os.getenv('BIOMCP_VERSION', '0.1.1')

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/chat', methods=['POST'])
def chat():
    data = request.json
    message = data.get('message', '')
    
    def generate():
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        
        async def process():
            async with OllamaBioMCPIntegration(
                OLLAMA_HOST, 
                OLLAMA_MODEL, 
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
            OLLAMA_MODEL, 
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
        "model": OLLAMA_MODEL,
        "biomcp_version": BIOMCP_VERSION
    })


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=False)