import asyncio
import json
from typing import Dict, Any, List, Optional
from dataclasses import dataclass
import logging

from mcp.client.session import ClientSession
from mcp.client.stdio import StdioServerParameters, stdio_client
from mcp.types import TextContent, Tool, Resource

logger = logging.getLogger(__name__)

@dataclass
class ToolResult:
    success: bool
    content: str
    error: Optional[str] = None

class BioMCPClient:
    def __init__(self, biomcp_version: str = "0.1.1"):
        self.biomcp_version = biomcp_version
        self.session: Optional[ClientSession] = None
        self._tools_cache: Optional[List[Tool]] = None
        self._resources_cache: Optional[List[Resource]] = None
        
    async def __aenter__(self):
        await self.connect()
        return self
        
    async def __aexit__(self, exc_type, exc_val, exc_tb):
        await self.disconnect()
        
    async def connect(self):
        """Establish connection to BioMCP server"""
        try:
            # Use uv to run biomcp-python package
            server_params = StdioServerParameters(
                command="uv",
                args=["run", "--with", f"biomcp-python=={self.biomcp_version}", "biomcp", "run"],
            )
            
            # Create stdio client connection
            self.read_stream, self.write_stream = await stdio_client(server_params).__aenter__()
            
            # Create and initialize session
            self.session = ClientSession(self.read_stream, self.write_stream)
            await self.session.__aenter__()
            await self.session.initialize()
            
            logger.info("Successfully connected to BioMCP server")
            
            # Cache tools and resources on connection
            await self._cache_tools_and_resources()
            
        except Exception as e:
            logger.error(f"Failed to connect to BioMCP: {e}")
            raise
    
    async def disconnect(self):
        """Close connection to BioMCP server"""
        if self.session:
            await self.session.__aexit__(None, None, None)
        if hasattr(self, 'read_stream'):
            # Close stdio streams
            pass
            
    async def _cache_tools_and_resources(self):
        """Cache available tools and resources"""
        if self.session:
            tool_result = await self.session.list_tools()
            self._tools_cache = tool_result.tools
            
            resource_result = await self.session.list_resources()
            self._resources_cache = resource_result.resources
    
    async def get_available_tools(self) -> List[Dict[str, Any]]:
        """Get list of available BioMCP tools"""
        if not self._tools_cache:
            if self.session:
                tool_result = await self.session.list_tools()
                self._tools_cache = tool_result.tools
            else:
                return []
        
        # Convert tools to dict format
        tools_list = []
        for tool in self._tools_cache:
            tool_dict = {
                "name": tool.name,
                "description": tool.description,
                "parameters": []
            }
            
            if tool.inputSchema:
                # Parse parameters from input schema
                properties = tool.inputSchema.get("properties", {})
                required = tool.inputSchema.get("required", [])
                
                for param_name, param_schema in properties.items():
                    param_info = {
                        "name": param_name,
                        "type": param_schema.get("type", "string"),
                        "description": param_schema.get("description", ""),
                        "required": param_name in required
                    }
                    tool_dict["parameters"].append(param_info)
            
            tools_list.append(tool_dict)
        
        return tools_list
    
    async def execute_tool(self, tool_name: str, parameters: Dict[str, Any]) -> ToolResult:
        """Execute a BioMCP tool with given parameters"""
        if not self.session:
            return ToolResult(
                success=False,
                content="",
                error="Not connected to BioMCP server"
            )
        
        try:
            # Call the tool through MCP protocol
            result = await self.session.call_tool(tool_name, parameters)
            
            if result.isError:
                return ToolResult(
                    success=False,
                    content="",
                    error=f"Tool execution error: {result.content}"
                )
            
            # Extract text content from result
            content_text = ""
            if result.content:
                for content_block in result.content:
                    if isinstance(content_block, TextContent):
                        content_text += content_block.text + "\n"
            
            return ToolResult(
                success=True,
                content=content_text.strip()
            )
            
        except Exception as e:
            logger.error(f"Error executing tool {tool_name}: {e}")
            return ToolResult(
                success=False,
                content="",
                error=str(e)
            )
    
    def format_tools_for_llm(self, tools: List[Dict[str, Any]]) -> str:
        """Format tool descriptions for LLM context"""
        tool_descriptions = []
        
        for tool in tools:
            desc = f"Tool: {tool['name']}\n"
            desc += f"Description: {tool['description']}\n"
            
            if tool['parameters']:
                desc += "Parameters:\n"
                for param in tool['parameters']:
                    req = " (required)" if param.get('required', False) else " (optional)"
                    desc += f"  - {param['name']} ({param['type']}){req}: {param.get('description', '')}\n"
            
            tool_descriptions.append(desc)
        
        return "\n".join(tool_descriptions)