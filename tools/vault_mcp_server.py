"""
MCP server exposing semantic search over the Obsidian vault.

OpenCode calls the search_vault tool to retrieve relevant notes.

Requires mcp>=2 (uses MCPServer, formerly FastMCP).
"""
from __future__ import annotations

from pathlib import Path

import chromadb
from chromadb.config import Settings
from mcp.server.mcpserver import MCPServer
from sentence_transformers import SentenceTransformer

PROJECT_ROOT = Path(__file__).resolve().parent.parent
INDEX_DIR = PROJECT_ROOT / "tools" / ".vault_index"
MODEL_NAME = "all-MiniLM-L6-v2"

# Initialize once at startup
_model = None
_collection = None


def _get_model():
    global _model
    if _model is None:
        _model = SentenceTransformer(MODEL_NAME)
    return _model


def _get_collection():
    global _collection
    if _collection is None:
        if not INDEX_DIR.exists():
            raise RuntimeError(
                f"Vault index not found at {INDEX_DIR}. "
                f"Run `python tools/vault_index.py` first."
            )
        client = chromadb.PersistentClient(
            path=str(INDEX_DIR),
            settings=Settings(anonymized_telemetry=False),
        )
        _collection = client.get_collection("vault")
    return _collection


mcp = MCPServer("obsidian-vault-search")


@mcp.tool()
def search_vault(query: str, k: int = 5) -> str:
    """Search the Obsidian vault for notes relevant to the query.

    Returns the top-k most relevant chunks with source file paths.

    Args:
        query: Natural language search query
        k: Number of results to return (default 5, max 20)
    """
    k = max(1, min(k, 20))
    model = _get_model()
    collection = _get_collection()

    query_embedding = model.encode([query])[0].tolist()

    results = collection.query(
        query_embeddings=[query_embedding],
        n_results=k,
    )

    if not results["documents"] or not results["documents"][0]:
        return f"No results for query: {query}"

    output = [f"Top {len(results['documents'][0])} results for: {query}\n"]
    for i, (doc, meta, dist) in enumerate(zip(
        results["documents"][0],
        results["metadatas"][0],
        results["distances"][0],
    )):
        similarity = 1 - dist
        output.append(
            f"\n--- Result {i + 1} (similarity {similarity:.3f}) ---\n"
            f"Source: {meta['source']} (chunk {meta['chunk_idx']})\n"
            f"{doc}\n"
        )

    return "\n".join(output)


@mcp.tool()
def vault_stats() -> str:
    """Return statistics about the vault index (number of files, chunks)."""
    try:
        collection = _get_collection()
    except RuntimeError as e:
        return str(e)
    n = collection.count()
    return f"Vault index: {n} chunks indexed from {INDEX_DIR}"


if __name__ == "__main__":
    mcp.run()