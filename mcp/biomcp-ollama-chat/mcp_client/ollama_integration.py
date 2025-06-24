import re
import json
import httpx
import asyncio
from typing import Dict, Any, AsyncGenerator, Tuple, Optional, List
import logging

from .biomcp_client import BioMCPClient, ToolResult

logger = logging.getLogger(__name__)

class OllamaBioMCPIntegration:
    def __init__(self, ollama_host: str, model_name: str, biomcp_version: str = "0.1.1"):
        self.ollama_host = ollama_host
        self.model_name = model_name
        self.biomcp_client = BioMCPClient(biomcp_version)
        self.http_client = httpx.AsyncClient(timeout=60.0)
        
    async def __aenter__(self):
        await self.biomcp_client.connect()
        return self
        
    async def __aexit__(self, exc_type, exc_val, exc_tb):
        await self.biomcp_client.disconnect()
        await self.http_client.aclose()
    
    def _create_system_prompt(self, tools: List[Dict[str, Any]]) -> str:
        """Create system prompt with BioMCP tools context"""
        tools_context = self.biomcp_client.format_tools_for_llm(tools)
        
        # Load instructions and researcher prompts
        import os
        prompts_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "prompts")
        
        instructions_prompt = ""
        researcher_prompt = ""
        
        try:
            with open(os.path.join(prompts_dir, "instructions.md"), 'r') as f:
                instructions_prompt = f.read()
        except FileNotFoundError:
            pass
            
        try:
            with open(os.path.join(prompts_dir, "researcher.md"), 'r') as f:
                researcher_prompt = f.read()
        except FileNotFoundError:
            pass
        
        # Combine the prompts
        combined_prompt = f"""{instructions_prompt}

{researcher_prompt}

## Available BioMCP Tools:
{tools_context}

## Tool Call Format
Use this exact format for tool calls:
tool_name: {{"parameter_name": "value"}}

Examples:
- article_searcher: {{"query": "BRAF V600E clinical implications"}}
- variant_searcher: {{"query": "breast cancer"}}
- trial_searcher: {{"query": "melanoma treatment"}}"""
        
        return combined_prompt
    
    def _extract_tool_calls(self, text: str) -> List[Tuple[str, Dict[str, Any]]]:
        """Extract tool calls from LLM response"""
        tool_calls = []
        
        valid_tools = ["variant_details", "variant_searcher", "article_details", "article_searcher", "trial_searcher"]
        
        # Pattern to match tool_name: {json} (primary format)
        pattern1 = r'(\w+):\s*(\{[^}]*\})'
        for match in re.finditer(pattern1, text):
            try:
                tool_name = match.group(1).strip()
                tool_json = json.loads(match.group(2))
                
                if tool_name in valid_tools:
                    tool_calls.append((tool_name, tool_json))
            except json.JSONDecodeError:
                logger.error(f"Failed to parse tool call JSON: {match.group(2)}")
                continue
        
        # Pattern to match tool_name(param="value") (fallback format)
        pattern2 = r'(\w+)\s*\(\s*(\w+)\s*=\s*["\']([^"\']*)["\']'
        for match in re.finditer(pattern2, text):
            try:
                tool_name = match.group(1).strip()
                param_name = match.group(2).strip()
                param_value = match.group(3).strip()
                
                if tool_name in valid_tools:
                    tool_calls.append((tool_name, {param_name: param_value}))
            except Exception as e:
                logger.error(f"Failed to parse function call format: {e}")
                continue
        
        return tool_calls
    
    async def process_message(self, user_message: str) -> AsyncGenerator[Dict[str, Any], None]:
        """Process user message with Ollama and BioMCP integration"""
        
        # Get available tools
        tools = await self.biomcp_client.get_available_tools()
        system_prompt = self._create_system_prompt(tools)
        
        # Prepare the conversation
        conversation = []
        conversation.append({"role": "system", "content": system_prompt})
        conversation.append({"role": "user", "content": user_message})
        
        # Create a single prompt from the conversation
        full_prompt = ""
        for msg in conversation:
            if msg["role"] == "system":
                full_prompt += f"System: {msg['content']}\n\n"
            elif msg["role"] == "user":
                full_prompt += f"User: {msg['content']}\n\nAssistant: "
        
        # Stream response from Ollama
        async with self.http_client.stream(
            "POST",
            f"{self.ollama_host}/api/generate",
            json={
                "model": self.model_name,
                "prompt": full_prompt,
                "stream": True
            }
        ) as response:
            accumulated_response = ""
            
            async for line in response.aiter_lines():
                if line:
                    try:
                        chunk = json.loads(line)
                        if "response" in chunk:
                            content = chunk["response"]
                            accumulated_response += content
                            
                            # Check for tool calls in accumulated response
                            tool_calls = self._extract_tool_calls(accumulated_response)
                            
                            # Always yield content, but filter out tool calls
                            if not any(tool_name in content for tool_name in ["variant_details:", "variant_searcher:", "article_details:", "article_searcher:", "trial_searcher:"]):
                                yield {"type": "content", "content": content}
                            
                            # If message is complete and contains tool calls, execute them
                            if chunk.get("done", False) and tool_calls:
                                for tool_name, parameters in tool_calls:
                                    yield {"type": "tool_call", "tool": tool_name, "parameters": parameters}
                                    
                                    # Execute the tool
                                    result = await self.biomcp_client.execute_tool(tool_name, parameters)
                                    
                                    if result.success:
                                        yield {
                                            "type": "tool_result",
                                            "tool": tool_name,
                                            "result": result.content
                                        }
                                        
                                        # Add tool result to conversation and get LLM's interpretation
                                        conversation.append({
                                            "role": "assistant",
                                            "content": accumulated_response
                                        })
                                        conversation.append({
                                            "role": "user",
                                            "content": f"Tool '{tool_name}' returned:\n\n{result.content}\n\nPlease interpret and explain these results."
                                        })
                                        
                                        # Create interpretation prompt
                                        interp_prompt = f"Tool '{tool_name}' returned this data:\n\n{result.content}\n\nPlease interpret and explain these results in a clear, helpful way:"
                                        
                                        # Get interpretation from LLM
                                        async with self.http_client.stream(
                                            "POST",
                                            f"{self.ollama_host}/api/generate",
                                            json={
                                                "model": self.model_name,
                                                "prompt": interp_prompt,
                                                "stream": True
                                            }
                                        ) as interp_response:
                                            async for interp_line in interp_response.aiter_lines():
                                                if interp_line:
                                                    try:
                                                        interp_chunk = json.loads(interp_line)
                                                        if "response" in interp_chunk:
                                                            yield {
                                                                "type": "interpretation",
                                                                "content": interp_chunk["response"]
                                                            }
                                                    except json.JSONDecodeError:
                                                        continue
                                    else:
                                        yield {
                                            "type": "tool_error",
                                            "tool": tool_name,
                                            "error": result.error
                                        }
                        
                    except json.JSONDecodeError:
                        continue