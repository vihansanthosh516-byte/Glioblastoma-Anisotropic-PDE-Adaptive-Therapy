"""Quick test: does the vault RAG return sensible results?"""
import chromadb
from chromadb.config import Settings
from sentence_transformers import SentenceTransformer

PROJECT_ROOT = None
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
INDEX_DIR = PROJECT_ROOT / "tools" / ".vault_index"

client = chromadb.PersistentClient(
    path=str(INDEX_DIR),
    settings=Settings(anonymized_telemetry=False),
)
collection = client.get_collection("vault")
print(f"Indexed chunks: {collection.count()}")
print()

model = SentenceTransformer("all-MiniLM-L6-v2")

queries = [
    "E_MAX_RATIO decision",
    "adaptive therapy resistance",
    "threshold setpoint high rho",
    "script 44 results",
    "negative rho responders",
]

for query in queries:
    print(f"=== Query: {query} ===")
    emb = model.encode([query])[0].tolist()
    results = collection.query(query_embeddings=[emb], n_results=3)
    for doc, meta, dist in zip(
        results["documents"][0],
        results["metadatas"][0],
        results["distances"][0],
    ):
        sim = 1 - dist
        source = meta["source"]
        chunk = meta["chunk_idx"]
        preview = doc[:140].replace("\n", " ")
        print(f"  [{sim:.3f}] {source} chunk {chunk}")
        print(f"         {preview}...")
    print()