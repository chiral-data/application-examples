import re
import json
import httpx
import asyncio
from typing import Dict, Any, AsyncGenerator, Tuple, Optional
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
        
        return f"""You are a helpful bioinformatics assistant with access to specialized tools through BioMCP.

Available BioMCP Tools:
{tools_context}

When users ask questions that could benefit from these tools:
1. Identify which tool would be most helpful
2. Gather any required parameters from the user if not provided
3. Use the tool by responding with: TOOL_CALL: {{"tool": "tool_name", "parameters": {{"param1": "value1"}}}}
4. Wait for the tool result and then interpret and explain it clearly

Important:
- Always use tools when they would provide specific, accurate information
- The TOOL_CALL must be valid JSON on a single line
- Explain what you're doing before calling a tool
- After receiving results, provide clear interpretation

Be helpful and explain bioinformatics concepts clearly."""
    
    def _extract_tool_calls(self, text: str) -> List[Tuple[str, Dict[str, Any]]]:
        """Extract tool calls from LLM response"""
        tool_calls = []
        
        # Pattern to match TOOL_CALL: {json}
        pattern = r'TOOL_CALL:\s*(\{[^}]+\})'
        matches = re.finditer(pattern, text)
        
        for match in matches:
            try:
                tool_json = json.loads(match.group(1))
                tool_name = tool_json.get("tool")
                parameters = tool_json.get("parameters", {})
                
                if tool_name:
                    tool_calls.append((tool_name, parameters))
            except json.JSONDecodeError:
                logger.error(f"Failed to parse tool call JSON: {match.group(1)}")
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
        
        # Stream response from Ollama
        async with self.http_client.stream(
            "POST",
            f"{self.ollama_host}/api/chat",
            json={
                "model": self.model_name,
                "messages": conversation,
                "stream": True
            }
        ) as response:
            accumulated_response = ""
            
            async for line in response.aiter_lines():
                if line:
                    try:
                        chunk = json.loads(line)
                        if "message" in chunk and "content" in chunk["message"]:
                            content = chunk["message"]["content"]
                            accumulated_response += content
                            
                            # Check for tool calls in accumulated response
                            tool_calls = self._extract_tool_calls(accumulated_response)
                            
                            if tool_calls and not chunk.get("done", False):
                                # Don't yield the tool call syntax to user
                                if "TOOL_CALL:" not in content:
                                    yield {"type": "content", "content": content}
                            else:
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
                                        
                                        # Get interpretation from LLM
                                        async with self.http_client.stream(
                                            "POST",
                                            f"{self.ollama_host}/api/chat",
                                            json={
                                                "model": self.model_name,
                                                "messages": conversation,
                                                "stream": True
                                            }
                                        ) as interp_response:
                                            async for interp_line in interp_response.aiter_lines():
                                                if interp_line:
                                                    try:
                                                        interp_chunk = json.loads(interp_line)
                                                        if "message" in interp_chunk and "content" in interp_chunk["message"]:
                                                            yield {
                                                                "type": "interpretation",
                                                                "content": interp_chunk["message"]["content"]
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