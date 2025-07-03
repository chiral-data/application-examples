"""
LLM utility for formatting RAG responses using Ollama
"""
import os
import requests
import json
from typing import Dict, Any

class ResponseFormatter:
    def __init__(self):
        self.ollama_host = os.getenv("OLLAMA_HOST", "localhost")
        self.ollama_port = os.getenv("OLLAMA_PORT", "11434")
        self.model_name = os.getenv("OLLAMA_MODEL", "llama3.2:latest")
        self.base_url = f"http://{self.ollama_host}:{self.ollama_port}"
        
        # Load the formatting prompt
        self.prompt_template = self._load_prompt_template()
    
    def _load_prompt_template(self) -> str:
        """Load the formatting prompt template"""
        try:
            prompt_path = os.path.join(os.path.dirname(__file__), "format_prompt.txt")
            with open(prompt_path, 'r', encoding='utf-8') as f:
                return f.read().strip()
        except FileNotFoundError:
            # Fallback prompt if file not found
            return """You are a scientific document assistant. Format this RAG response for better readability:
- Remove technical SMILES strings and molecular formulas from main text
- Use simple structure references like "Structure 1"
- Include chemical properties like toxicity and activity
- Keep it scientific but readable

RAW RESPONSE TO FORMAT:
{raw_response}

Please reformat this response:"""
    
    def format_response(self, raw_response: Dict[str, Any]) -> str:
        """Format a RAG response using Ollama LLM"""
        try:
            # Convert response to string if it's a dict
            if isinstance(raw_response, dict):
                response_text = json.dumps(raw_response, indent=2)
            else:
                response_text = str(raw_response)
            
            # Create the formatting prompt
            formatting_prompt = self.prompt_template.format(raw_response=response_text)
            
            # Make request to Ollama
            response = self._call_ollama(formatting_prompt)
            
            if response:
                return response
            else:
                # Fallback to basic formatting if LLM fails
                return self._basic_format(raw_response)
                
        except Exception as e:
            print(f"Error formatting response with LLM: {e}")
            return self._basic_format(raw_response)
    
    def _call_ollama(self, prompt: str) -> str:
        """Make API call to Ollama"""
        try:
            payload = {
                "model": self.model_name,
                "prompt": prompt,
                "stream": False,
                "options": {
                    "temperature": 0.3,  # Low temperature for consistent formatting
                    "top_p": 0.9,
                    "max_tokens": 1000
                }
            }
            
            response = requests.post(
                f"{self.base_url}/api/generate",
                json=payload,
                timeout=30
            )
            
            if response.status_code == 200:
                result = response.json()
                return result.get('response', '').strip()
            else:
                print(f"Ollama API error: {response.status_code}")
                return None
                
        except requests.exceptions.RequestException as e:
            print(f"Request error calling Ollama: {e}")
            return None
    
    def _basic_format(self, raw_response: Dict[str, Any]) -> str:
        """Basic fallback formatting without LLM"""
        if isinstance(raw_response, str):
            return raw_response
        
        if isinstance(raw_response, dict) and 'answer' in raw_response:
            answer = raw_response.get('answer', '')
            sources = raw_response.get('sources', [])
            
            # Basic cleanup - remove long SMILES strings
            import re
            cleaned = re.sub(r'SMILES:\s*[A-Za-z0-9@\[\]()=+\-#\.:,\\\/]+', 'Structure ID', answer)
            cleaned = re.sub(r'Formula:\s*[A-Za-z0-9]+\s*MW:\s*[\d.]+', '', cleaned)
            cleaned = re.sub(r'\n\s*\n+', '\n\n', cleaned).strip()
            
            # Add basic source info
            if sources:
                cleaned += f"\n\n*Sources: {len(sources)} document sections*"
            
            return cleaned
        
        return str(raw_response)

# Global formatter instance
formatter = ResponseFormatter()

def format_rag_response(response: Dict[str, Any]) -> str:
    """Main function to format RAG responses"""
    return formatter.format_response(response)