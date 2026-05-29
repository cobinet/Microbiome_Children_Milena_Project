#!/usr/bin/env python3
# =============================================================================
# Author      : Dr. rer. nat. Guillermo G. Torres <guigotoe@gmail.com>
# Project     : MILENA — Microbiome of Children (Pastoral Kids Cohort)
# Repository  : https://github.com/cobinet/Microbiome_Children_Milena_Project
# Script      : docs/graph/build_pipeline_graph.py
# Description :
#   Builds a directed knowledge graph of the MILENA analysis pipeline using
#   networkx. Nodes represent data inputs, scripts, and output artefacts.
#   Edges represent data-flow relationships (read / transform / produce).
#   Exports the graph as:
#     - pipeline_graph.dot  (Graphviz DOT source — for use in publications)
#     - pipeline_graph.png  (rendered PNG via matplotlib — included in README)
# Usage :
#   python docs/graph/build_pipeline_graph.py
# Dependencies : networkx, matplotlib (both in requirements.txt)
# =============================================================================

import networkx as nx
import matplotlib
matplotlib.use("Agg")          # headless rendering — no display required
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# ── Node catalog ─────────────────────────────────────────────────────────────
# Each tuple: (node_id, label, node_type)
#   node_type ∈ {input, script, intermediate, output}
NODES = [
    # Raw / external data inputs
    ("ASV_table",      "ASV_table.tsv\n(relative abundance)",           "input"),
    ("ASV_tax",        "ASV_tax_species.user.tsv\n(taxonomy)",           "input"),
    ("EC_raw",         "EC_pred_metagenome_unstrat.tsv\n(PICRUSt2 EC)", "input"),
    ("KO_raw",         "KO_pred_metagenome_unstrat.tsv\n(PICRUSt2 KO)", "input"),
    ("MetaCyc_raw",    "METACYC_path_abun_unstrat.tsv\n(PICRUSt2 MC)",  "input"),
    ("metadata_xlsx",  "data_formicrobiome.xlsx\n(clinical metadata)",  "input"),

    # Scripts
    ("nb_preproc",     "funAnalysis_metagPastoKids.ipynb\n[Python/Jupyter]", "script"),
    ("r_assoc",        "AssocAnalysis.R\n[R]",                               "script"),
    ("r_graph",        "build_pipeline_graph.py\n[Python]",                  "script"),

    # Intermediate data (preprocessing outputs — not committed)
    ("asv_clr",        "asv_table_clr_transformed.tsv",              "intermediate"),
    ("KO_clr",         "KO_Functions_clr_transformed.tsv",           "intermediate"),
    ("EC_clr",         "EC_Functions_clr_transformed.tsv",           "intermediate"),
    ("MC_clr",         "MetaCyc_Pathways_clr_transformed.tsv",       "intermediate"),
    ("meta_clean",     "metadata_clean.csv",                         "intermediate"),
    ("tax_clean",      "asv_taxonomy_clean.csv",                     "intermediate"),

    # Final outputs (not committed)
    ("html_table",     "{TARGET}_publication_table.html",            "output"),
    ("lollipop_png",   "{TARGET}_lollipop_*.png",                    "output"),
    ("volcano_png",    "volcano_{BLOCK}_{contrast}.png",             "output"),
    ("csv_results",    "{TARGET}_ALL_features_all_contrasts.csv",    "output"),
    ("graph_dot",      "pipeline_graph.dot",                         "output"),
    ("graph_png",      "pipeline_graph.png",                         "output"),
]

# ── Edge catalog ─────────────────────────────────────────────────────────────
# Each tuple: (source_id, target_id, edge_label)
EDGES = [
    # Notebook reads raw inputs
    ("ASV_table",   "nb_preproc",  "reads"),
    ("ASV_tax",     "nb_preproc",  "reads"),
    ("EC_raw",      "nb_preproc",  "reads"),
    ("KO_raw",      "nb_preproc",  "reads"),
    ("MetaCyc_raw", "nb_preproc",  "reads"),
    ("metadata_xlsx","nb_preproc", "reads"),

    # Notebook produces CLR-transformed intermediates + cleaned metadata
    ("nb_preproc",  "asv_clr",    "CLR transform"),
    ("nb_preproc",  "KO_clr",     "CLR transform"),
    ("nb_preproc",  "EC_clr",     "CLR transform"),
    ("nb_preproc",  "MC_clr",     "CLR transform"),
    ("nb_preproc",  "meta_clean", "QC + clean"),
    ("nb_preproc",  "tax_clean",  "QC + clean"),

    # R script reads CLR intermediates and metadata
    ("asv_clr",    "r_assoc",  "reads"),
    ("KO_clr",     "r_assoc",  "reads"),
    ("EC_clr",     "r_assoc",  "reads"),
    ("MC_clr",     "r_assoc",  "reads"),
    ("meta_clean", "r_assoc",  "reads"),
    ("tax_clean",  "r_assoc",  "reads"),

    # R script produces final outputs
    ("r_assoc",  "html_table",  "gt export"),
    ("r_assoc",  "lollipop_png","ggplot2"),
    ("r_assoc",  "volcano_png", "ggplot2"),
    ("r_assoc",  "csv_results", "write_csv"),

    # Graph builder produces DOT + PNG
    ("r_graph",  "graph_dot",  "networkx export"),
    ("r_graph",  "graph_png",  "matplotlib render"),
]

# ── Colour scheme ─────────────────────────────────────────────────────────────
COLOR_MAP = {
    "input":        "#AED6F1",   # light blue — raw data
    "script":       "#A9DFBF",   # light green — code
    "intermediate": "#FAD7A0",   # light orange — processed data
    "output":       "#D2B4DE",   # light purple — final results
}

SHAPE_MAP = {
    "input":        "s",    # square
    "script":       "D",    # diamond
    "intermediate": "o",    # circle
    "output":       "^",    # triangle
}


def build_graph():
    G = nx.DiGraph()

    # Add nodes with attributes
    for node_id, label, ntype in NODES:
        G.add_node(node_id, label=label, ntype=ntype)

    # Add edges
    for src, dst, elabel in EDGES:
        G.add_edge(src, dst, label=elabel)

    return G


def export_dot(G: nx.DiGraph, dot_path: Path):
    """Write a Graphviz DOT file from the networkx graph."""
    lines = ["digraph MILENA_Pipeline {",
             '    graph [rankdir=LR, fontname="Helvetica", splines=ortho];',
             '    node  [fontname="Helvetica", fontsize=10];',
             '    edge  [fontname="Helvetica", fontsize=8];',
             ""]

    ntype_shape = {"input": "box", "script": "diamond",
                   "intermediate": "ellipse", "output": "parallelogram"}
    ntype_color = {"input": "#AED6F1", "script": "#A9DFBF",
                   "intermediate": "#FAD7A0", "output": "#D2B4DE"}

    for node_id, data in G.nodes(data=True):
        label = data["label"].replace("\n", "\\n")
        shape = ntype_shape.get(data["ntype"], "box")
        color = ntype_color.get(data["ntype"], "white")
        lines.append(
            f'    "{node_id}" [label="{label}", shape={shape}, '
            f'style=filled, fillcolor="{color}"];'
        )

    lines.append("")
    for src, dst, data in G.edges(data=True):
        elabel = data.get("label", "")
        lines.append(f'    "{src}" -> "{dst}" [label="{elabel}"];')

    lines.append("}")
    dot_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"DOT file written: {dot_path}")


def export_png(G: nx.DiGraph, png_path: Path):
    """Render the graph as a PNG using matplotlib."""
    fig, ax = plt.subplots(figsize=(20, 12))
    ax.set_title("MILENA Analysis Pipeline — Data Flow Graph",
                 fontsize=14, fontweight="bold", pad=15)

    # Hierarchical layout via dot-like positioning
    try:
        pos = nx.nx_agraph.graphviz_layout(G, prog="dot", args="-Grankdir=LR")
    except Exception:
        # Fall back to spring layout if Graphviz is not installed
        pos = nx.spring_layout(G, seed=42, k=3.0)

    node_types  = [G.nodes[n]["ntype"] for n in G.nodes()]
    node_colors = [COLOR_MAP[t] for t in node_types]

    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                           node_size=2200, alpha=0.95)
    nx.draw_networkx_labels(
        G, pos, ax=ax,
        labels={n: G.nodes[n]["label"] for n in G.nodes()},
        font_size=7, font_family="monospace"
    )
    nx.draw_networkx_edges(G, pos, ax=ax, arrows=True,
                           arrowsize=18, width=1.2,
                           edge_color="#555555", connectionstyle="arc3,rad=0.08")
    edge_labels = {(u, v): d["label"] for u, v, d in G.edges(data=True)}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels,
                                 ax=ax, font_size=6, font_color="#333333")

    # Legend
    legend_patches = [
        mpatches.Patch(color=COLOR_MAP["input"],        label="Raw input data"),
        mpatches.Patch(color=COLOR_MAP["script"],       label="Analysis script"),
        mpatches.Patch(color=COLOR_MAP["intermediate"], label="Intermediate data (CLR)"),
        mpatches.Patch(color=COLOR_MAP["output"],       label="Final output artefact"),
    ]
    ax.legend(handles=legend_patches, loc="lower left", fontsize=9,
              framealpha=0.9, title="Node type")
    ax.axis("off")

    plt.tight_layout()
    plt.savefig(png_path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"PNG written: {png_path}")


if __name__ == "__main__":
    out_dir = Path(__file__).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    G = build_graph()
    print(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    export_dot(G, out_dir / "pipeline_graph.dot")
    export_png(G, out_dir / "pipeline_graph.png")
    print("Done.")
