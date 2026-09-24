"""Test the MCP server end-to-end using the MCP client SDK."""
import asyncio
import sys
from pathlib import Path

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PYTHON = PROJECT_ROOT / ".venv" / "Scripts" / "python.exe"
SERVER = PROJECT_ROOT / "tools" / "vault_mcp_server.py"


async def main():
    params = StdioServerParameters(
        command=str(PYTHON),
        args=[str(SERVER)],
    )

    async with stdio_client(params) as (read, write):
        async with ClientSession(read, write) as session:
            await session.initialize()

            # List available tools
            tools = await session.list_tools()
            print("Available tools:")
            for t in tools.tools:
                print(f"  - {t.name}: {t.description[:80]}")
            print()

            # Call vault_stats
            result = await session.call_tool("vault_stats", {})
            print("vault_stats result:")
            for content in result.content:
                if hasattr(content, "text"):
                    print(f"  {content.text}")
            print()

            # Call search_vault
            result = await session.call_tool(
                "search_vault",
                {"query": "E_MAX_RATIO decision", "k": 2},
            )
            print("search_vault('E_MAX_RATIO decision', k=2):")
            for content in result.content:
                if hasattr(content, "text"):
                    text = content.text
                    print(text[:500])
                    if len(text) > 500:
                        print(f"  ... [{len(text) - 500} more chars]")


if __name__ == "__main__":
    asyncio.run(main())