"""
Index the Obsidian vault into a Chroma vector database.

Run this whenever vault content changes significantly:
    python tools/vault_index.py

It reads every .md file under vault/, chunks it, embeds the chunks,
and stores them in tools/.vault_index/.
"""
from __future__ import annotations

import hashlib
import shutil
from pathlib import Path

import chromadb
from chromadb.config import Settings
from langchain_text_splitters import RecursiveCharacterTextSplitter
from sentence_transformers import SentenceTransformer

# Paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent
VAULT_DIR = PROJECT_ROOT / "vault"
INDEX_DIR = PROJECT_ROOT / "tools" / ".vault_index"

# Embedding model (runs locally, no API key)
MODEL_NAME = "all-MiniLM-L6-v2"
CHUNK_SIZE = 800
CHUNK_OVERLAP = 150


def find_markdown_files(root: Path) -> list[Path]:
    """Find all .md files, excluding .obsidian/ and other noise."""
    ignore_dirs = {".obsidian", ".git", ".trash", "node_modules", ".vault_index"}
    files = []
    for path in root.rglob("*.md"):
        if any(part in ignore_dirs for part in path.parts):
            continue
        files.append(path)
    return sorted(files)


def chunk_text(text: str, source: str) -> list[dict]:
    """Split text into chunks with metadata."""
    splitter = RecursiveCharacterTextSplitter(
        chunk_size=CHUNK_SIZE,
        chunk_overlap=CHUNK_OVERLAP,
        separators=["\n## ", "\n### ", "\n\n", "\n", " "],
    )
    chunks = splitter.split_text(text)
    return [
        {"text": c, "source": source, "chunk_idx": i, "total_chunks": len(chunks)}
        for i, c in enumerate(chunks)
    ]


def main():
    print(f"Indexing vault: {VAULT_DIR}")
    if not VAULT_DIR.exists():
        raise SystemExit(f"Vault not found: {VAULT_DIR}")

    # Fresh index each time (delete old)
    if INDEX_DIR.exists():
        shutil.rmtree(INDEX_DIR)
    INDEX_DIR.mkdir(parents=True)

    # Load embedding model
    print(f"Loading embedding model: {MODEL_NAME}")
    model = SentenceTransformer(MODEL_NAME)

    # Chroma client
    client = chromadb.PersistentClient(
        path=str(INDEX_DIR),
        settings=Settings(anonymized_telemetry=False),
    )
    collection = client.create_collection(
        name="vault",
        metadata={"hnsw:space": "cosine"},
    )

    # Find and chunk files
    md_files = find_markdown_files(VAULT_DIR)
    print(f"Found {len(md_files)} markdown files")

    all_chunks = []
    for f in md_files:
        try:
            text = f.read_text(encoding="utf-8")
        except Exception as e:
            print(f"  SKIP {f.name}: {e}")
            continue
        rel = f.relative_to(PROJECT_ROOT)
        chunks = chunk_text(text, str(rel))
        all_chunks.extend(chunks)

    print(f"Total chunks: {len(all_chunks)}")

    if not all_chunks:
        print("No chunks to index. Add notes to vault/ and re-run.")
        return

    # Embed and store
    print("Embedding chunks...")
    texts = [c["text"] for c in all_chunks]
    embeddings = model.encode(texts, show_progress_bar=True, batch_size=32)

    ids = [
        hashlib.sha1(f"{c['source']}:{c['chunk_idx']}".encode()).hexdigest()[:16]
        for c in all_chunks
    ]
    metadatas = [
        {"source": c["source"], "chunk_idx": c["chunk_idx"]} for c in all_chunks
    ]

    collection.add(
        ids=ids,
        embeddings=embeddings.tolist(),
        documents=texts,
        metadatas=metadatas,
    )

    print(f"Indexed {len(ids)} chunks into {INDEX_DIR}")
    print("Done.")


if __name__ == "__main__":
    main()