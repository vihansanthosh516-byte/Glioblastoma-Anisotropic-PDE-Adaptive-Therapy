"""Send a raw MCP initialize request to the server and print everything."""
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PYTHON = PROJECT_ROOT / ".venv" / "Scripts" / "python.exe"
SERVER = PROJECT_ROOT / "tools" / "vault_mcp_server.py"

# Minimal MCP initialize request (JSON-RPC 2.0)
request = (
    '{"jsonrpc":"2.0","id":1,"method":"initialize",'
    '"params":{"protocolVersion":"2024-11-05","capabilities":{},'
    '"clientInfo":{"name":"debug","version":"0.1"}}}\n'
)

proc = subprocess.Popen(
    [str(PYTHON), str(SERVER)],
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True,
    encoding="utf-8",
)

try:
    stdout, stderr = proc.communicate(input=request, timeout=15)
except subprocess.TimeoutExpired:
    proc.kill()
    stdout, stderr = proc.communicate()

print("=== STDOUT ===")
print(stdout if stdout else "(empty)")
print()
print("=== STDERR ===")
print(stderr if stderr else "(empty)")
print()
print(f"=== RETURN CODE === {proc.returncode}")