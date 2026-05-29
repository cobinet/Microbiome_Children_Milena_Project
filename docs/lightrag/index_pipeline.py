#!/usr/bin/env python3
# =============================================================================
# Author      : Dr. rer. nat. Guillermo G. Torres <guigotoe@gmail.com>
# Project     : MILENA — Microbiome of Children (Pastoral Kids Cohort)
# Repository  : https://github.com/cobinet/Microbiome_Children_Milena_Project
# Script      : docs/lightrag/index_pipeline.py
# Description :
#   Indexes all MILENA analysis scripts into a local LightRAG knowledge base
#   (naive mode — no API key required). After indexing, runs a set of
#   structured queries to extract:
#     (a) full pipeline outline
#     (b) per-script purpose summaries
#     (c) complete dependency lists (R packages / Python libraries)
#   Output is written to docs/lightrag/pipeline_outline.md for inclusion
#   in the repository README.
#
# Usage :
#   cd <repo_root>
#   python docs/lightrag/index_pipeline.py
#
# Dependencies : lightrag-hku (pip install lightrag-hku)
#   LightRAG runs fully locally in "naive" mode — no LLM API key is needed
#   for the indexing step. Query results use the local nano-vectordb backend.
# =============================================================================

import asyncio
import json
from pathlib import Path

# LightRAG local imports
from lightrag import LightRAG, QueryParam
from lightrag.llm.openai import openai_complete_if_cache, openai_embed
from lightrag.utils import EmbeddingFunc

# ── Configuration ─────────────────────────────────────────────────────────────
REPO_ROOT   = Path(__file__).parent.parent.parent   # <repo_root>/
WORKING_DIR = Path(__file__).parent / "rag_storage" # local vector store
WORKING_DIR.mkdir(parents=True, exist_ok=True)

# Scripts to index
SCRIPT_PATHS = [
    REPO_ROOT / "scripts" / "AssocAnalysis.R",
    # Notebook source extracted to plain text for indexing
    REPO_ROOT / "scripts" / "funAnalysis_metagPastoKids_source.txt",
]

# Structured queries to run after indexing
QUERIES = {
    "pipeline_outline": (
        "Describe the complete analysis pipeline step by step, "
        "from raw data inputs through preprocessing, statistical analysis, "
        "and final outputs."
    ),
    "script_purposes": (
        "For each script (AssocAnalysis.R and funAnalysis_metagPastoKids), "
        "provide a one-paragraph description of its purpose, inputs, outputs, "
        "and statistical methods used."
    ),
    "dependencies": (
        "List all external R packages and Python libraries used across "
        "all scripts, grouped by language."
    ),
    "replication_guide": (
        "Write step-by-step instructions for a computational biologist "
        "to replicate this analysis from scratch, including data acquisition, "
        "environment setup, and script execution order."
    ),
}


def extract_notebook_source(nb_path: Path, out_path: Path):
    """Extract code and markdown cells from .ipynb as plain text for indexing."""
    if not nb_path.exists():
        print(f"  Notebook not found: {nb_path} — skipping extraction.")
        return

    with open(nb_path) as f:
        nb = json.load(f)

    lines = [f"# Notebook: {nb_path.name}\n"]
    for i, cell in enumerate(nb.get("cells", [])):
        ctype = cell.get("cell_type", "code")
        src   = "".join(cell.get("source", []))
        if not src.strip():
            continue
        lines.append(f"\n## Cell {i} [{ctype}]\n{src}\n")

    out_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"  Notebook source extracted → {out_path}")


def collect_documents() -> list[str]:
    """Read all indexed scripts into a list of strings."""
    docs = []
    for path in SCRIPT_PATHS:
        if path.exists():
            docs.append(path.read_text(encoding="utf-8", errors="replace"))
            print(f"  Loaded: {path.name} ({len(docs[-1])} chars)")
        else:
            print(f"  MISSING: {path} — skipping.")
    return docs


async def run_lightrag_naive(docs: list[str], queries: dict) -> dict:
    """
    Index documents and run queries in LightRAG naive mode.
    Naive mode uses direct embedding search without graph construction —
    suitable for local use without an LLM API key for the index step.
    """
    # Minimal embedding setup using a local sentence-level hash
    # (replace with a real embedding function if an API key is available)
    async def mock_embed(texts):
        import hashlib, numpy as np
        vecs = []
        for t in texts:
            seed = int(hashlib.md5(t[:200].encode()).hexdigest(), 16) % (2**31)
            rng  = np.random.default_rng(seed)
            vecs.append(rng.standard_normal(384).astype(np.float32).tolist())
        return vecs

    rag = LightRAG(
        working_dir   = str(WORKING_DIR),
        embedding_func = EmbeddingFunc(
            embedding_dim     = 384,
            max_token_size    = 512,
            func              = mock_embed,
        ),
    )

    print("\nIndexing documents...")
    for i, doc in enumerate(docs):
        await rag.ainsert(doc)
        print(f"  Indexed document {i+1}/{len(docs)}")

    results = {}
    print("\nRunning queries...")
    for qname, qtext in queries.items():
        print(f"  Query: {qname}")
        try:
            result = await rag.aquery(qtext, param=QueryParam(mode="naive"))
            results[qname] = result
        except Exception as e:
            results[qname] = f"[Query failed: {e}]"

    return results


def write_outline(results: dict, out_path: Path):
    """Save LightRAG query results as a Markdown document."""
    lines = [
        "# MILENA Pipeline — LightRAG Knowledge Extraction\n",
        "> Auto-generated by `docs/lightrag/index_pipeline.py`  \n",
        "> Author: Dr. rer. nat. Guillermo G. Torres\n\n",
    ]
    section_titles = {
        "pipeline_outline":  "## 1. Full Pipeline Outline",
        "script_purposes":   "## 2. Script Purposes",
        "dependencies":      "## 3. Dependencies",
        "replication_guide": "## 4. Replication Guide",
    }
    for key, title in section_titles.items():
        lines.append(f"{title}\n\n")
        lines.append(results.get(key, "_No result generated._"))
        lines.append("\n\n---\n\n")

    out_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nOutline written → {out_path}")


def main():
    print("=== MILENA LightRAG Pipeline Indexer ===\n")

    # Extract notebook source text
    nb_src = REPO_ROOT / "scripts" / "funAnalysis_metagPastoKids_source.txt"
    nb_orig = REPO_ROOT / "scripts" / "funAnalysis_metagPastoKids.ipynb"
    extract_notebook_source(nb_orig, nb_src)

    docs = collect_documents()
    if not docs:
        print("ERROR: No documents found to index. Check SCRIPT_PATHS.")
        return

    results = asyncio.run(run_lightrag_naive(docs, QUERIES))

    out_path = Path(__file__).parent / "pipeline_outline.md"
    write_outline(results, out_path)
    print("\nDone.")


if __name__ == "__main__":
    main()
