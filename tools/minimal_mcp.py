"""Minimal MCP server — no dependencies beyond mcp package."""
from mcp.server.fastmcp import FastMCP

mcp = FastMCP("test-server")


@mcp.tool()
def hello(name: str) -> str:
    """Say hello to someone."""
    return f"Hello, {name}!"


if __name__ == "__main__":
    print("Server starting...", file=__import__("sys").stderr)
    mcp.run()
    print("Server exited.", file=__import__("sys").stderr)