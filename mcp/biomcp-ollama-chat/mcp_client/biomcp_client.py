import asyncio
import json
import subprocess
import logging
from typing import Dict, Any, List, Optional
from dataclasses import dataclass

logger = logging.getLogger(__name__)

@dataclass
class ToolResult:
    success: bool
    content: str
    error: Optional[str] = None

class BioMCPClient:
    """BioMCP client that uses direct subprocess calls for reliable tool execution"""
    
    def __init__(self, biomcp_version: str = "0.1.1"):
        self.biomcp_version = biomcp_version
        
    async def __aenter__(self):
        return self
        
    async def __aexit__(self, exc_type, exc_val, exc_tb):
        pass
        
    async def connect(self):
        """Connect to BioMCP (no-op for subprocess-based client)"""
        logger.info("BioMCPClient connected")
        pass
    
    async def disconnect(self):
        """Disconnect from BioMCP (no-op for subprocess-based client)"""
        pass
        
    async def get_available_tools(self) -> List[Dict[str, Any]]:
        """Get list of available BioMCP tools"""
        # Define the tools based on the actual biomcp CLI capabilities
        tools = [
            {
                "name": "variant_searcher",
                "description": "Search for genetic variants using various query terms",
                "parameters": [
                    {"name": "query", "required": True, "description": "Search query for variants"}
                ]
            },
            {
                "name": "variant_details", 
                "description": "Get detailed information about a specific variant",
                "parameters": [
                    {"name": "variant_id", "required": True, "description": "Variant ID (e.g., rs113488022)"}
                ]
            },
            {
                "name": "article_searcher",
                "description": "Search for biomedical articles and publications",
                "parameters": [
                    {"name": "query", "required": True, "description": "Search query for articles"}
                ]
            },
            {
                "name": "article_details",
                "description": "Get detailed information about a specific article",
                "parameters": [
                    {"name": "pmid", "required": True, "description": "PubMed ID of the article"}
                ]
            },
            {
                "name": "trial_searcher",
                "description": "Search for clinical trials",
                "parameters": [
                    {"name": "query", "required": True, "description": "Search query for clinical trials"}
                ]
            }
        ]
        return tools
    
    async def execute_tool(self, tool_name: str, parameters: Dict[str, Any]) -> ToolResult:
        """Execute a BioMCP tool using subprocess"""
        try:
            # Build the command based on tool name and actual biomcp CLI structure
            if tool_name == "variant_details":
                variant_id = parameters.get("variant_id", "")
                cmd = ["biomcp", "variant", "get", variant_id]
            elif tool_name == "variant_searcher":
                query = parameters.get("query", "")
                cmd = ["biomcp", "variant", "search"]
                
                # Use general search with pathogenic variants if no specific parameters found
                cmd.extend(["--significance", "pathogenic"])
                cmd.extend(["--size", "20"])
            elif tool_name == "article_details":
                pmid = parameters.get("pmid", "")
                cmd = ["biomcp", "article", "get", pmid]
            elif tool_name == "article_searcher":
                query = parameters.get("query", "")
                cmd = ["biomcp", "article", "search"]
                # Parse query for relevant search terms
                if "BRAF" in query.upper():
                    cmd.extend(["--gene", "BRAF"])
                if "V600E" in query.upper():
                    cmd.extend(["--variant", "V600E"])
                elif "mutation" in query.lower():
                    cmd.extend(["--keyword", "mutation"])
                # Extract other keywords
                keywords = [word for word in query.split() if len(word) > 3 and word.lower() not in ["what", "are", "the", "and", "mutation", "implications"]]
                for keyword in keywords[:3]:  # Limit to 3 keywords
                    cmd.extend(["--keyword", keyword])
            elif tool_name == "trial_searcher":
                query = parameters.get("query", "")
                cmd = ["biomcp", "trial", "search"]
                
                # Use general search term
                cmd.extend(["--term", query])
                cmd.extend(["--status", "open"])
            else:
                return ToolResult(
                    success=False,
                    content="",
                    error=f"Unknown tool: {tool_name}"
                )
            
            logger.info(f"Executing biomcp command: {' '.join(cmd)}")
            
            # Execute the command with timeout
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )
            
            try:
                stdout, stderr = await asyncio.wait_for(
                    process.communicate(), 
                    timeout=30.0
                )
            except asyncio.TimeoutError:
                process.kill()
                await process.wait()
                return ToolResult(
                    success=False,
                    content="",
                    error="Tool execution timed out"
                )
            
            if process.returncode == 0:
                output = stdout.decode('utf-8').strip()
                return ToolResult(
                    success=True,
                    content=output,
                    error=None
                )
            else:
                error_msg = stderr.decode('utf-8').strip()
                return ToolResult(
                    success=False,
                    content="",
                    error=f"Tool failed with exit code {process.returncode}: {error_msg}"
                )
                
        except Exception as e:
            logger.error(f"Error executing tool {tool_name}: {e}")
            return ToolResult(
                success=False,
                content="",
                error=str(e)
            )
    
    def format_tools_for_llm(self, tools: List[Dict[str, Any]]) -> str:
        """Format tools for LLM context"""
        tool_descriptions = []
        for tool in tools:
            params = []
            if tool.get("parameters"):
                for param in tool["parameters"]:
                    param_desc = f"{param['name']}"
                    if param.get("required"):
                        param_desc += " (required)"
                    if param.get("description"):
                        param_desc += f": {param['description']}"
                    params.append(param_desc)
            
            tool_desc = f"- {tool['name']}: {tool.get('description', 'No description')}"
            if params:
                tool_desc += f"\n  Parameters: {', '.join(params)}"
            tool_descriptions.append(tool_desc)
        
        return "Available BioMCP tools:\n" + "\n".join(tool_descriptions)