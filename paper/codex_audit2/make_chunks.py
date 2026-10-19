"""Builds paper/codex_audit2/chunks.json and progress.md: a fixed, complete work list for the old-track audit.
Deterministic. Re-running it must give the same chunks (it overwrites progress.md only if --reset is passed).
"""
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = Path(__file__).resolve().parent
SKIP_MD = ("graphify-out/", ".claude/", ".kilo/", "paper/codex_audit")


def ledger_entries():
    ids = []
    for line in (ROOT / "vault" / "PAPER_FINDINGS_LEDGER.md").read_text(encoding="utf-8", errors="replace").split("\n"):
        m = re.match(r"^(\d+[a-z]{0,2})\. ", line)
        if m:
            ids.append(m.group(1))
    return ids


def main():
    chunks = []
    ids = ledger_entries()
    cut = ids.index("2br") if "2br" in ids else len(ids)
    old_ids = ids[:cut]
    for i in range(0, len(old_ids), 6):
        chunks.append(dict(kind="ledger", items=old_ids[i:i + 6]))
    md = subprocess.check_output(["git", "ls-files", "*.md"], cwd=ROOT, text=True).split("\n")
    md = sorted(f for f in md if f and not f.startswith(SKIP_MD))
    cur, size = [], 0
    for f in md:
        s = (ROOT / f).stat().st_size if (ROOT / f).exists() else 0
        if cur and (size + s > 30000 or len(cur) >= 5):
            chunks.append(dict(kind="markdown", items=cur)); cur, size = [], 0
        cur.append(f); size += s
    if cur:
        chunks.append(dict(kind="markdown", items=cur))
    scripts = sorted(str(p.relative_to(ROOT)).replace("\\", "/") for p in (ROOT / "src").rglob("*.py") if "__pycache__" not in p.parts)
    for i in range(0, len(scripts), 8):
        chunks.append(dict(kind="script", items=scripts[i:i + 8]))
    tests = sorted(str(p.relative_to(ROOT)).replace("\\", "/") for p in (ROOT / "tests").glob("*.py"))
    chunks.append(dict(kind="tests", items=tests))
    for n, c in enumerate(chunks, 1):
        c["id"] = f"{n:03d}"
    (OUT / "chunks.json").write_text(json.dumps(chunks, indent=1))
    prog = OUT / "progress.md"
    if "--reset" in sys.argv or not prog.exists():
        lines = ["# Progress (update after EVERY chunk)", "", "| chunk | kind | items | status |", "|---|---|---|---|"]
        for c in chunks:
            lines.append(f"| {c['id']} | {c['kind']} | {len(c['items'])} | TODO |")
        prog.write_text("\n".join(lines) + "\n")
    print(len(chunks), "chunks;", {k: sum(1 for c in chunks if c["kind"] == k) for k in ("ledger", "markdown", "script", "tests")},
          "| ledger entries before 2br:", len(old_ids), "| md files:", len(md), "| scripts:", len(scripts))


if __name__ == "__main__":
    main()
