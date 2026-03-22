#!/usr/bin/env python3
"""
visualize_model.py
==================
Comprehensive visualization module for a COBRApy genome-scale metabolic model.

Produces (all outputs written to ``out_dir``):
  1. summary_stats.png          – bar chart of reactions / metabolites / genes /
                                   blocked reactions
  2. compartment_breakdown.png  – pie chart of metabolite compartments
  3. degree_distribution.png    – histogram of reaction connectivity
  4. subsystem_barchart.png     – horizontal bar chart of reactions per subsystem
  5. flux_distribution.png      – FBA flux bar chart (top-N active reactions)
  6. reaction_network.html      – interactive Plotly network graph
  7. dashboard.html             – interactive Plotly multi-panel dashboard
  8. mutation_analysis.html     – interactive knockout + knockin gene analysis
  9. escher_map.html            – interactive Escher metabolic map (if escher
                                   is installed and a base-map is available)
  10. model_summary.csv         – machine-readable summary statistics table
  11. index.html                – unified single-page report hub with all
                                   figures, interactive iframes, and environment
                                   editor

Usage (standalone)
------------------
    python src/visualize_model.py --model output/model.xml --out-dir output/viz

Usage (from generate_model.py)
-------------------------------
    from visualize_model import run_all_visualizations
    run_all_visualizations(model, out_dir=Path("output/viz"), log=log)
"""

from __future__ import annotations

import argparse
import logging
import warnings
from pathlib import Path
from typing import TYPE_CHECKING

# ---------------------------------------------------------------------------
# Soft imports – warn if missing rather than crashing at import time
# ---------------------------------------------------------------------------
try:
    import matplotlib
    matplotlib.use("Agg")          # headless / non-interactive backend
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    _HAS_MPL = True
except ImportError:
    _HAS_MPL = False

try:
    import seaborn as sns
    _HAS_SNS = True
except ImportError:
    _HAS_SNS = False

try:
    import networkx as nx
    _HAS_NX = True
except ImportError:
    _HAS_NX = False

try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    _HAS_PLOTLY = True
except ImportError:
    _HAS_PLOTLY = False

try:
    import pandas as pd
    from collections import Counter
    _HAS_PD = True
except ImportError:
    _HAS_PD = False

try:
    import escher
    _HAS_ESCHER = True
except ImportError:
    _HAS_ESCHER = False

if TYPE_CHECKING:
    import cobra

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _check_mpl(log: logging.Logger) -> bool:
    if not _HAS_MPL:
        log.warning("matplotlib not installed – skipping static plots.")
    return _HAS_MPL


def _apply_style() -> None:
    """Apply a clean, publication-quality style."""
    if _HAS_SNS:
        sns.set_theme(style="whitegrid", palette="muted", font_scale=1.1)
    else:
        plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid"
                      in plt.style.available else "ggplot")


def _save(fig: "plt.Figure", path: Path, log: logging.Logger, dpi: int = 150) -> None:
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    log.info("  Saved: %s", path)


# ---------------------------------------------------------------------------
# 1. Summary statistics bar chart
# ---------------------------------------------------------------------------

def plot_summary_stats(
    model: "cobra.Model",
    blocked: list[str],
    fba_value: float | None,
    out_dir: Path,
    log: logging.Logger,
) -> None:
    """Bar chart: reactions, metabolites, genes, blocked reactions."""
    if not _check_mpl(log):
        return

    _apply_style()

    n_rxn   = len(model.reactions)
    n_met   = len(model.metabolites)
    n_gene  = len(model.genes)
    n_blk   = len(blocked)
    n_active = n_rxn - n_blk

    labels  = ["Reactions", "Metabolites", "Genes", "Active\nReactions", "Blocked\nReactions"]
    values  = [n_rxn,       n_met,         n_gene,  n_active,            n_blk]
    colors  = ["#4C72B0",   "#55A868",     "#C44E52","#8172B2",           "#CCB974"]

    fig, ax = plt.subplots(figsize=(9, 5))
    bars = ax.bar(labels, values, color=colors, width=0.55, edgecolor="white", linewidth=0.8)

    for bar, val in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(values) * 0.015,
            f"{val:,}",
            ha="center", va="bottom", fontsize=10, fontweight="bold",
        )

    title = f"Model: {model.id or 'GSMM'}"
    if fba_value is not None:
        title += f"   |   FBA growth = {fba_value:.4f}"
    ax.set_title(title, fontsize=13, fontweight="bold", pad=14)
    ax.set_ylabel("Count", fontsize=11)
    ax.set_ylim(0, max(values) * 1.15)
    ax.spines[["top", "right"]].set_visible(False)

    _save(fig, out_dir / "summary_stats.png", log)


# ---------------------------------------------------------------------------
# 2. Compartment breakdown pie chart
# ---------------------------------------------------------------------------

def plot_compartment_breakdown(
    model: "cobra.Model",
    out_dir: Path,
    log: logging.Logger,
) -> None:
    """Pie chart of metabolite counts per compartment."""
    if not _check_mpl(log):
        return

    _apply_style()

    from collections import Counter
    comp_counts = Counter(m.compartment for m in model.metabolites)

    if not comp_counts:
        log.warning("  No compartment data – skipping compartment plot.")
        return

    # Friendly names where available
    comp_names = model.compartments  # dict {id: name}
    labels = [comp_names.get(c, c) or c for c in comp_counts.keys()]
    sizes  = list(comp_counts.values())

    palette = (sns.color_palette("pastel", len(sizes))
               if _HAS_SNS else None)

    fig, ax = plt.subplots(figsize=(7, 7))
    pie_result = ax.pie(
        sizes,
        labels=labels,
        autopct="%1.1f%%",
        startangle=140,
        colors=palette,
        wedgeprops=dict(linewidth=0.7, edgecolor="white"),
        pctdistance=0.82,
    )
    # pie() returns (wedges, texts) or (wedges, texts, autotexts) depending on autopct
    if len(pie_result) == 3:
        for at in pie_result[2]:
            at.set_fontsize(9)

    ax.set_title("Metabolite distribution by compartment",
                 fontsize=13, fontweight="bold", pad=16)

    _save(fig, out_dir / "compartment_breakdown.png", log)


# ---------------------------------------------------------------------------
# 3. Reaction-degree distribution histogram
# ---------------------------------------------------------------------------

def plot_degree_distribution(
    model: "cobra.Model",
    out_dir: Path,
    log: logging.Logger,
) -> None:
    """Histogram of how many metabolites participate in each reaction."""
    if not _check_mpl(log):
        return

    _apply_style()

    degrees = [len(rxn.metabolites) for rxn in model.reactions]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(degrees, bins=30, color="#4C72B0", edgecolor="white", linewidth=0.6)
    ax.axvline(sum(degrees) / max(len(degrees), 1), color="#C44E52",
               linestyle="--", linewidth=1.5, label=f"Mean = {sum(degrees)/max(len(degrees),1):.1f}")
    ax.set_xlabel("Number of metabolites per reaction", fontsize=11)
    ax.set_ylabel("Number of reactions", fontsize=11)
    ax.set_title("Reaction connectivity (degree) distribution",
                 fontsize=13, fontweight="bold", pad=14)
    ax.legend(fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)

    _save(fig, out_dir / "degree_distribution.png", log)


# ---------------------------------------------------------------------------
# 4. Subsystem bar chart
# ---------------------------------------------------------------------------

def plot_subsystem_barchart(
    model: "cobra.Model",
    out_dir: Path,
    log: logging.Logger,
    top_n: int = 25,
) -> None:
    """Horizontal bar chart of the top-N subsystems by reaction count."""
    if not _check_mpl(log):
        return

    _apply_style()

    from collections import Counter
    subsystems = Counter(
        rxn.subsystem for rxn in model.reactions if rxn.subsystem
    )
    if not subsystems:
        log.warning("  No subsystem annotations – skipping subsystem plot.")
        return

    top = subsystems.most_common(top_n)
    labels = [s for s, _ in top][::-1]
    counts = [c for _, c in top][::-1]

    palette = (sns.color_palette("Blues_d", len(counts))
               if _HAS_SNS else ["#4C72B0"] * len(counts))

    fig, ax = plt.subplots(figsize=(10, max(5, len(labels) * 0.35)))
    bars = ax.barh(labels, counts, color=palette, edgecolor="white", linewidth=0.6)

    for bar, val in zip(bars, counts):
        ax.text(val + max(counts) * 0.01, bar.get_y() + bar.get_height() / 2,
                str(val), va="center", fontsize=8)

    ax.set_xlabel("Number of reactions", fontsize=11)
    ax.set_title(f"Top {min(top_n, len(subsystems))} subsystems by reaction count",
                 fontsize=13, fontweight="bold", pad=14)
    ax.spines[["top", "right"]].set_visible(False)

    _save(fig, out_dir / "subsystem_barchart.png", log)


# ---------------------------------------------------------------------------
# 5. FBA flux distribution
# ---------------------------------------------------------------------------

def plot_flux_distribution(
    model: "cobra.Model",
    out_dir: Path,
    log: logging.Logger,
    top_n: int = 30,
) -> None:
    """Bar chart of the top-N reactions by absolute FBA flux."""
    if not _check_mpl(log):
        return

    try:
        solution = model.optimize()
        if solution.status != "optimal":
            log.warning("  FBA status '%s' – skipping flux distribution plot.", solution.status)
            return
    except Exception as exc:
        log.warning("  FBA failed (%s) – skipping flux distribution plot.", exc)
        return

    _apply_style()

    fluxes = solution.fluxes
    top = fluxes.abs().nlargest(top_n)
    rxn_ids = top.index.tolist()
    flux_vals = fluxes[rxn_ids].tolist()

    colors = ["#4C72B0" if v >= 0 else "#C44E52" for v in flux_vals]

    # Use short reaction names if available
    def _short(rid: str) -> str:
        rxn = model.reactions.get_by_id(rid)
        name = rxn.name or rid
        return name[:35] + "…" if len(name) > 35 else name

    labels = [_short(r) for r in rxn_ids]

    fig, ax = plt.subplots(figsize=(10, max(5, len(labels) * 0.35)))
    bars = ax.barh(labels[::-1], flux_vals[::-1], color=colors[::-1],
                   edgecolor="white", linewidth=0.5)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Flux (mmol gDW⁻¹ h⁻¹)", fontsize=11)
    ax.set_title(f"Top {len(labels)} reactions by absolute FBA flux",
                 fontsize=13, fontweight="bold", pad=14)

    pos_patch = mpatches.Patch(color="#4C72B0", label="Forward flux")
    neg_patch = mpatches.Patch(color="#C44E52", label="Reverse flux")
    ax.legend(handles=[pos_patch, neg_patch], fontsize=9, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)

    _save(fig, out_dir / "flux_distribution.png", log)


# ---------------------------------------------------------------------------
# 6. Interactive Plotly network graph
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 6b. Interactive Plotly dashboard (overview of all key statistics)
# ---------------------------------------------------------------------------

def plot_reaction_network(
    model: "cobra.Model",
    out_dir: Path,
    log: logging.Logger,
    max_nodes: int = 300,
) -> None:
    """
    Build a semantically-clustered reaction-metabolite network and export it as
    interactive HTML.

    Enhancements:
      - Reactions grouped by subsystem/pathway with consistent colors.
      - Force-directed layout with subsystem hub attraction.
      - Dual view toggle (detailed nodes vs subsystem modules) as semantic zoom.
    """
    if not _HAS_PLOTLY:
        log.warning("plotly not installed – skipping interactive network graph.")
        return
    if not _HAS_NX:
        log.warning("networkx not installed – skipping interactive network graph.")
        return

    from collections import Counter, defaultdict
    import json as _json
    import math
    import random

    log.info("  Building reaction network graph (semantic clustering) …")

    reactions = list(model.reactions)
    if len(reactions) > max_nodes:
        random.seed(42)
        reactions = random.sample(reactions, max_nodes)
        log.info("  (sampled %d / %d reactions for legibility)", max_nodes, len(model.reactions))

    G = nx.DiGraph()
    for rxn in reactions:
        subsystem = (rxn.subsystem or "Unknown subsystem").strip() or "Unknown subsystem"
        G.add_node(
            rxn.id,
            node_type="reaction",
            label=rxn.name or rxn.id,
            subsystem=subsystem,
        )
        for met in rxn.reactants:
            G.add_node(
                met.id,
                node_type="metabolite",
                label=met.name or met.id,
                compartment=met.compartment or "?",
            )
            G.add_edge(met.id, rxn.id, direction="import")
        for met in rxn.products:
            G.add_node(
                met.id,
                node_type="metabolite",
                label=met.name or met.id,
                compartment=met.compartment or "?",
            )
            G.add_edge(rxn.id, met.id, direction="export")

    rxn_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "reaction"]
    met_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "metabolite"]

    if not rxn_nodes:
        log.warning("  No reaction nodes available – skipping reaction network plot.")
        return

    subsystem_counts = Counter(G.nodes[n].get("subsystem", "Unknown subsystem") for n in rxn_nodes)
    subsystem_order = [s for s, _ in subsystem_counts.most_common()]
    n_subsystems = max(len(subsystem_order), 1)

    palette = [
        "#2563eb", "#14b8a6", "#f97316", "#a855f7", "#ef4444", "#22c55e", "#f59e0b",
        "#0ea5e9", "#d946ef", "#84cc16", "#4f46e5", "#e11d48", "#06b6d4", "#7c3aed",
        "#16a34a", "#ea580c", "#3b82f6", "#ca8a04", "#0f766e", "#6d28d9",
    ]
    subsystem_color = {sub: palette[i % len(palette)] for i, sub in enumerate(subsystem_order)}

    # Build an auxiliary layout graph with subsystem hubs to create semantic clusters.
    layout_graph = nx.Graph()
    layout_graph.add_nodes_from(G.nodes())
    for u, v in G.edges():
        layout_graph.add_edge(u, v, weight=1.0)

    hub_radius = max(2.0, 0.62 * math.sqrt(n_subsystems) + 0.8)
    hub_nodes: dict[str, str] = {}
    fixed_hubs: list[str] = []
    initial_pos: dict[str, tuple[float, float]] = {}
    rnd = random.Random(42)

    for i, sub in enumerate(subsystem_order):
        hub = f"__subsystem_hub__{i}"
        hub_nodes[sub] = hub
        layout_graph.add_node(hub)
        angle = (2.0 * math.pi * i) / n_subsystems
        hx = hub_radius * math.cos(angle)
        hy = hub_radius * math.sin(angle)
        initial_pos[hub] = (hx, hy)
        fixed_hubs.append(hub)

    for rxn in rxn_nodes:
        sub = G.nodes[rxn].get("subsystem", "Unknown subsystem")
        hub = hub_nodes.get(sub)
        if hub:
            layout_graph.add_edge(rxn, hub, weight=6.0)
            hx, hy = initial_pos[hub]
            initial_pos[rxn] = (hx + rnd.uniform(-0.25, 0.25), hy + rnd.uniform(-0.25, 0.25))

    for met in met_nodes:
        neighbors = [n for n in (list(G.predecessors(met)) + list(G.successors(met))) if n in initial_pos]
        if neighbors:
            mx = sum(initial_pos[n][0] for n in neighbors) / len(neighbors)
            my = sum(initial_pos[n][1] for n in neighbors) / len(neighbors)
            initial_pos[met] = (mx + rnd.uniform(-0.18, 0.18), my + rnd.uniform(-0.18, 0.18))
        else:
            initial_pos[met] = (rnd.uniform(-0.4, 0.4), rnd.uniform(-0.4, 0.4))

    try:
        pos_all = nx.spring_layout(
            layout_graph,
            seed=42,
            pos=initial_pos,
            fixed=fixed_hubs,
            iterations=280,
            k=1.35 / math.sqrt(max(layout_graph.number_of_nodes(), 4)),
            weight="weight",
        )
    except Exception:
        try:
            pos_all = nx.kamada_kawai_layout(G)
        except Exception:
            pos_all = nx.spring_layout(G, seed=42)

    pos = {n: pos_all[n] for n in G.nodes() if n in pos_all}

    # Assign metabolites to dominant neighboring subsystem (for semantic coloring).
    met_dominant_subsystem: dict[str, str] = {}
    for met in met_nodes:
        neighbor_reactions = [
            n for n in (list(G.predecessors(met)) + list(G.successors(met)))
            if G.nodes[n].get("node_type") == "reaction"
        ]
        if not neighbor_reactions:
            met_dominant_subsystem[met] = "Unknown subsystem"
            continue
        vote = Counter(G.nodes[r].get("subsystem", "Unknown subsystem") for r in neighbor_reactions)
        met_dominant_subsystem[met] = vote.most_common(1)[0][0]

    # Detailed view traces
    edge_x, edge_y = [], []
    rxn_to_metabolites: dict[str, set[str]] = defaultdict(set)
    rxn_to_edge_segments: dict[str, list[tuple[float, float, float, float]]] = defaultdict(list)
    for u, v in G.edges():
        if u not in pos or v not in pos:
            continue
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        edge_x += [x0, x1, None]
        edge_y += [y0, y1, None]

        u_type = G.nodes[u].get("node_type")
        v_type = G.nodes[v].get("node_type")
        if u_type == "reaction" and v_type == "metabolite":
            rxn_to_metabolites[u].add(v)
            rxn_to_edge_segments[u].append((x0, y0, x1, y1))
        elif u_type == "metabolite" and v_type == "reaction":
            rxn_to_metabolites[v].add(u)
            rxn_to_edge_segments[v].append((x0, y0, x1, y1))

    edge_trace = go.Scatter(
        x=edge_x,
        y=edge_y,
        mode="lines",
        line=dict(width=0.7, color="#9ca3af"),
        hoverinfo="none",
        name="Reaction links",
        showlegend=False,
        visible=True,
    )

    highlight_edge_trace = go.Scatter(
        x=[],
        y=[],
        mode="lines",
        line=dict(width=2.8, color="#ef4444"),
        hoverinfo="none",
        name="Selected pathways",
        showlegend=False,
        visible=True,
        opacity=0.0,
    )

    rxn_x, rxn_y, rxn_hover, rxn_color, rxn_ids_for_plot = [], [], [], [], []
    for n in rxn_nodes:
        if n not in pos:
            continue
        x, y = pos[n]
        data = G.nodes[n]
        subsystem = data.get("subsystem", "Unknown subsystem")
        rxn_x.append(x)
        rxn_y.append(y)
        rxn_ids_for_plot.append(n)
        rxn_color.append(subsystem_color.get(subsystem, "#334155"))
        rxn_hover.append(
            f"<b>{data.get('label', n)}</b><br>"
            f"Reaction ID: {n}<br>"
            f"Subsystem: {subsystem}<br>"
            f"In-degree: {G.in_degree(n)}<br>"
            f"Out-degree: {G.out_degree(n)}"
        )

    rxn_trace = go.Scatter(
        x=rxn_x,
        y=rxn_y,
        mode="markers",
        marker=dict(size=11, color=rxn_color, symbol="square", line=dict(width=0.8, color="#ffffff")),
        hovertext=rxn_hover,
        hoverinfo="text",
        customdata=rxn_ids_for_plot,
        name="Reactions",
        visible=True,
    )

    met_x, met_y, met_hover, met_color, met_size, met_ids_for_plot = [], [], [], [], [], []
    for n in met_nodes:
        if n not in pos:
            continue
        x, y = pos[n]
        data = G.nodes[n]
        dom_sub = met_dominant_subsystem.get(n, "Unknown subsystem")
        met_x.append(x)
        met_y.append(y)
        met_ids_for_plot.append(n)
        met_color.append(subsystem_color.get(dom_sub, "#94a3b8"))
        met_size.append(min(16, 6 + 1.7 * (G.in_degree(n) + G.out_degree(n))))
        met_hover.append(
            f"<b>{data.get('label', n)}</b><br>"
            f"Metabolite ID: {n}<br>"
            f"Compartment: {data.get('compartment', '?')}<br>"
            f"Dominant subsystem: {dom_sub}<br>"
            f"Connections: {G.in_degree(n) + G.out_degree(n)}"
        )

    met_trace = go.Scatter(
        x=met_x,
        y=met_y,
        mode="markers",
        marker=dict(size=met_size, color=met_color, symbol="circle", opacity=0.86,
                    line=dict(width=0.8, color="#f8fafc")),
        hovertext=met_hover,
        hoverinfo="text",
        customdata=met_ids_for_plot,
        name="Metabolites",
        visible=True,
    )

    # Subsystem module (semantic zoom) view
    module_links: dict[tuple[str, str], int] = defaultdict(int)
    for met in met_nodes:
        neighbors = [
            n for n in (list(G.predecessors(met)) + list(G.successors(met)))
            if G.nodes[n].get("node_type") == "reaction"
        ]
        subs = sorted({G.nodes[r].get("subsystem", "Unknown subsystem") for r in neighbors})
        for i in range(len(subs)):
            for j in range(i + 1, len(subs)):
                module_links[(subs[i], subs[j])] += 1

    module_pos: dict[str, tuple[float, float]] = {}
    for i, sub in enumerate(subsystem_order):
        angle = (2.0 * math.pi * i) / n_subsystems
        module_pos[sub] = ((hub_radius + 0.3) * math.cos(angle), (hub_radius + 0.3) * math.sin(angle))

    mod_edge_x, mod_edge_y = [], []
    for (s1, s2), _w in module_links.items():
        x0, y0 = module_pos[s1]
        x1, y1 = module_pos[s2]
        mod_edge_x += [x0, x1, None]
        mod_edge_y += [y0, y1, None]

    module_edge_trace = go.Scatter(
        x=mod_edge_x,
        y=mod_edge_y,
        mode="lines",
        line=dict(width=2.0, color="#94a3b8"),
        hoverinfo="none",
        showlegend=False,
        name="Module links",
        visible=False,
    )

    module_x, module_y, module_sizes, module_colors, module_hover = [], [], [], [], []
    for sub in subsystem_order:
        x, y = module_pos[sub]
        cnt = subsystem_counts[sub]
        module_x.append(x)
        module_y.append(y)
        module_sizes.append(min(95, 24 + 6 * math.sqrt(max(cnt, 1))))
        module_colors.append(subsystem_color.get(sub, "#64748b"))
        module_hover.append(
            f"<b>{sub}</b><br>"
            f"Reactions in module: {cnt}<br>"
            f"Inter-module links: {sum(1 for (a, b) in module_links if a == sub or b == sub)}"
        )

    module_trace = go.Scatter(
        x=module_x,
        y=module_y,
        mode="markers+text",
        text=[sub if len(sub) <= 28 else sub[:25] + "…" for sub in subsystem_order],
        textposition="middle center",
        textfont=dict(size=10, color="#0f172a"),
        marker=dict(size=module_sizes, color=module_colors, symbol="circle", opacity=0.92,
                    line=dict(width=2.5, color="#ffffff")),
        hovertext=module_hover,
        hoverinfo="text",
        name="Subsystem modules",
        visible=False,
    )

    fig = go.Figure(data=[edge_trace, highlight_edge_trace, met_trace, rxn_trace, module_edge_trace, module_trace])
    fig.update_layout(
        title=dict(
            text=(
                "Reaction–Metabolite Network (Semantic Clustered)<br>"
                f"<sup>{len(G.nodes())} nodes • {len(G.edges())} edges • {n_subsystems} subsystems</sup>"
            ),
            font=dict(size=16),
            x=0.5,
        ),
        showlegend=True,
        hovermode="closest",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor="#f8fafc",
        paper_bgcolor="#f8fafc",
        margin=dict(l=20, r=20, t=92, b=20),
        updatemenus=[
            dict(
                type="buttons",
                direction="right",
                x=0.01,
                y=1.15,
                xanchor="left",
                yanchor="top",
                buttons=[
                    dict(
                        label="Detailed view",
                        method="update",
                        args=[
                            {"visible": [True, True, True, True, False, False]},
                            {
                                "title.text": (
                                    "Reaction–Metabolite Network (Detailed semantic layout)<br>"
                                    f"<sup>{len(G.nodes())} nodes • {len(G.edges())} edges • {n_subsystems} subsystems</sup>"
                                )
                            },
                        ],
                    ),
                    dict(
                        label="Subsystem modules",
                        method="update",
                        args=[
                            {"visible": [False, False, False, False, True, True]},
                            {
                                "title.text": (
                                    "Reaction–Metabolite Network (Subsystem module view)<br>"
                                    f"<sup>{n_subsystems} pathway modules • switch back to expand details</sup>"
                                )
                            },
                        ],
                    ),
                ],
            )
        ],
        annotations=[
            dict(
                x=0.01,
                y=1.08,
                xref="paper",
                yref="paper",
                text=(
                    "<b>Semantic clustering:</b> reaction nodes are grouped by subsystem; "
                    "metabolites inherit dominant neighboring subsystem color."
                ),
                showarrow=False,
                align="left",
                font=dict(size=11, color="#334155"),
            )
        ],
    )

    out_html = out_dir / "reaction_network.html"
    reaction_post_script = f"""
    (function() {{
      const gd = document.getElementsByClassName('plotly-graph-div')[0];
      if (!gd || !window.Plotly) return;

      const rxnToMets = {_json.dumps({k: sorted(v) for k, v in rxn_to_metabolites.items()})};
      const rxnToSegments = {_json.dumps({k: v for k, v in rxn_to_edge_segments.items()})};

      const EDGE_TRACE_IDX = 0;
      const HILITE_EDGE_TRACE_IDX = 1;
      const MET_TRACE_IDX = 2;
      const RXN_TRACE_IDX = 3;

      function clearHighlight() {{
        const rxnCount = (gd.data[RXN_TRACE_IDX].x || []).length;
        const metCount = (gd.data[MET_TRACE_IDX].x || []).length;
        Plotly.restyle(gd, {{
          'marker.opacity': [Array(rxnCount).fill(1.0)],
          'marker.size': [Array(rxnCount).fill(11)]
        }}, [RXN_TRACE_IDX]);
        Plotly.restyle(gd, {{
          'marker.opacity': [Array(metCount).fill(0.86)]
        }}, [MET_TRACE_IDX]);
        Plotly.restyle(gd, {{'opacity': [1.0]}}, [EDGE_TRACE_IDX]);
        Plotly.restyle(gd, {{'x': [[]], 'y': [[]], 'opacity': [0.0]}}, [HILITE_EDGE_TRACE_IDX]);
        Plotly.relayout(gd, {{'xaxis.autorange': true, 'yaxis.autorange': true}});
      }}

      function applyHighlight(reactionIds) {{
        const selected = new Set((reactionIds || []).map(String));
        if (!selected.size) {{
          clearHighlight();
          return;
        }}

        const rxnCustom = gd.data[RXN_TRACE_IDX].customdata || [];
        const rxnOpacity = [];
        const rxnSizes = [];
        const selectedRxnXY = [];
        for (let i = 0; i < rxnCustom.length; i += 1) {{
          const rid = String(rxnCustom[i]);
          const on = selected.has(rid);
          rxnOpacity.push(on ? 1.0 : 0.10);
          rxnSizes.push(on ? 16 : 9);
          if (on) selectedRxnXY.push([gd.data[RXN_TRACE_IDX].x[i], gd.data[RXN_TRACE_IDX].y[i]]);
        }}
        Plotly.restyle(gd, {{
          'marker.opacity': [rxnOpacity],
          'marker.size': [rxnSizes]
        }}, [RXN_TRACE_IDX]);

        const linkedMets = new Set();
        selected.forEach(rid => (rxnToMets[rid] || []).forEach(mid => linkedMets.add(String(mid))));
        const metCustom = gd.data[MET_TRACE_IDX].customdata || [];
        const metOpacity = [];
        const selectedMetXY = [];
        for (let i = 0; i < metCustom.length; i += 1) {{
          const mid = String(metCustom[i]);
          const on = linkedMets.has(mid);
          metOpacity.push(on ? 0.96 : 0.08);
          if (on) selectedMetXY.push([gd.data[MET_TRACE_IDX].x[i], gd.data[MET_TRACE_IDX].y[i]]);
        }}
        Plotly.restyle(gd, {{'marker.opacity': [metOpacity]}}, [MET_TRACE_IDX]);

        const hx = [];
        const hy = [];
        selected.forEach(rid => {{
          (rxnToSegments[rid] || []).forEach(seg => {{
            hx.push(seg[0], seg[2], null);
            hy.push(seg[1], seg[3], null);
          }});
        }});
        Plotly.restyle(gd, {{'x': [hx], 'y': [hy], 'opacity': [hx.length ? 1.0 : 0.0]}}, [HILITE_EDGE_TRACE_IDX]);
        Plotly.restyle(gd, {{'opacity': [0.15]}}, [EDGE_TRACE_IDX]);

        const pts = selectedRxnXY.concat(selectedMetXY);
        if (pts.length) {{
          const xs = pts.map(p => p[0]);
          const ys = pts.map(p => p[1]);
          const minX = Math.min.apply(null, xs);
          const maxX = Math.max.apply(null, xs);
          const minY = Math.min.apply(null, ys);
          const maxY = Math.max.apply(null, ys);
          const padX = Math.max(0.2, (maxX - minX) * 0.35 + 0.12);
          const padY = Math.max(0.2, (maxY - minY) * 0.35 + 0.12);
          Plotly.relayout(gd, {{
            'xaxis.range': [minX - padX, maxX + padX],
            'yaxis.range': [minY - padY, maxY + padY]
          }});
        }}
      }}

      window.addEventListener('message', function(event) {{
        const msg = event.data || {{}};
        if (msg.type === 'gsmm-highlight-reactions') {{
          applyHighlight(msg.reactionIds || []);
        }}
        if (msg.type === 'gsmm-clear-reaction-highlight') {{
          clearHighlight();
        }}
      }});
    }})();
    """
    fig.write_html(
        str(out_html),
        include_plotlyjs="cdn",
        config={"displayModeBar": True, "scrollZoom": True, "displaylogo": False},
        post_script=reaction_post_script,
    )
    log.info("  Saved: %s", out_html)


def plot_dashboard(
    model: "cobra.Model",
    blocked: list[str],
    fba_value: float | None,
    out_dir: Path,
    log: logging.Logger,
) -> None:
    """
    Single self-contained HTML dashboard combining:
      • Summary stats bar chart
      • Compartment pie chart
      • Top-20 subsystems bar chart
      • FBA flux scatter
    """
    if not _HAS_PLOTLY:
        log.warning("plotly not installed – skipping dashboard.")
        return
    if not _HAS_PD:
        log.warning("pandas not installed – skipping dashboard.")
        return

    log.info("  Building interactive dashboard …")

    from collections import Counter

    # --- data ---
    n_rxn    = len(model.reactions)
    n_met    = len(model.metabolites)
    n_gene   = len(model.genes)
    n_blk    = len(blocked)
    n_active = n_rxn - n_blk

    comp_counts = Counter(m.compartment for m in model.metabolites)
    comp_names  = model.compartments
    comp_labels = [comp_names.get(c, c) or c for c in comp_counts]
    comp_sizes  = list(comp_counts.values())

    subsystems = Counter(r.subsystem for r in model.reactions if r.subsystem)
    top_sub = subsystems.most_common(20)
    sub_labels = [s for s, _ in top_sub]
    sub_counts = [c for _, c in top_sub]

    # FBA fluxes
    fba_df = None
    try:
        solution = model.optimize()
        if solution.status == "optimal":
            top_flux = solution.fluxes.abs().nlargest(30)
            fba_df = pd.DataFrame({
                "reaction": top_flux.index,
                "flux":     solution.fluxes[top_flux.index].values,
            })
    except Exception:
        pass

    # --- layout ---
    n_rows = 3 if fba_df is not None else 2
    row_heights = [0.28, 0.28, 0.44] if fba_df is not None else [0.45, 0.55]
    specs: list[list[dict]] = (
        [[{"type": "bar"}, {"type": "pie"}],
         [{"type": "bar", "colspan": 2}, None],
         [{"type": "bar", "colspan": 2}, None]]
        if fba_df is not None else
        [[{"type": "bar"}, {"type": "pie"}],
         [{"type": "bar", "colspan": 2}, None]]
    )

    fig = make_subplots(
        rows=n_rows, cols=2,
        subplot_titles=(
            "Model statistics",
            "Metabolites by compartment",
            f"Top {len(sub_labels)} subsystems",
            *(["Top 30 FBA fluxes"] if fba_df is not None else []),
        ),
        specs=specs,
        row_heights=row_heights,
        vertical_spacing=0.10,
        horizontal_spacing=0.08,
    )

    # Row 1 left – summary bar
    fig.add_trace(go.Bar(
        x=["Reactions", "Metabolites", "Genes", "Active Rxns", "Blocked Rxns"],
        y=[n_rxn, n_met, n_gene, n_active, n_blk],
        marker_color=["#4C72B0","#55A868","#C44E52","#8172B2","#CCB974"],
        text=[f"{v:,}" for v in [n_rxn, n_met, n_gene, n_active, n_blk]],
        textposition="outside",
        showlegend=False,
    ), row=1, col=1)

    # Row 1 right – compartment pie
    fig.add_trace(go.Pie(
        labels=comp_labels, values=comp_sizes,
        hole=0.35, showlegend=True,
        textinfo="percent+label",
    ), row=1, col=2)

    # Row 2 – subsystem bar
    fig.add_trace(go.Bar(
        x=sub_counts[::-1], y=sub_labels[::-1],
        orientation="h",
        marker_color="#4C72B0",
        showlegend=False,
        text=[str(c) for c in sub_counts[::-1]],
        textposition="outside",
    ), row=2, col=1)

    # Row 3 – flux bar (if available)
    if fba_df is not None:
        colors_flux = ["#4C72B0" if v >= 0 else "#C44E52" for v in fba_df["flux"]]
        fig.add_trace(go.Bar(
            x=fba_df["flux"].tolist()[::-1],
            y=fba_df["reaction"].tolist()[::-1],
            orientation="h",
            marker_color=colors_flux[::-1],
            showlegend=False,
        ), row=3, col=1)

    title_str = f"GSMM Dashboard — {model.id or 'model'}"
    if fba_value is not None:
        title_str += f"  |  FBA growth = {fba_value:.4f}"

    fig.update_layout(
        title_text=title_str,
        title_font_size=16,
        height=300 * n_rows + 100,
        paper_bgcolor="white",
        plot_bgcolor="white",
    )

    out_html = out_dir / "dashboard.html"
    fig.write_html(str(out_html))
    log.info("  Saved: %s", out_html)


# ---------------------------------------------------------------------------
# 7. Escher interactive map
# ---------------------------------------------------------------------------

def export_escher_map(
    model: "cobra.Model",
    out_dir: Path,
    log: logging.Logger,
) -> None:
    """
    Export an Escher interactive HTML map.

    Escher needs a JSON map to overlay reactions on.  We attempt to use the
    built-in "iJO1366.Central metabolism" map as a demo base; for a real
    organism you would supply your own map JSON.
    """
    if not _HAS_ESCHER:
        log.warning("escher not installed – skipping Escher map.")
        return

    log.info("  Building Escher map …")

    try:
        solution = model.optimize()
        reaction_data = (
            solution.fluxes.to_dict() if solution.status == "optimal" else None
        )
    except Exception:
        reaction_data = None

    try:
        builder = escher.Builder(
            map_name="iJO1366.Central metabolism",
            model=model,
            reaction_data=reaction_data,
            reaction_scale=[
                {"type": "min",  "color": "#C44E52", "size": 12},
                {"type": "zero", "color": "#eeeeee", "size":  5},
                {"type": "max",  "color": "#4C72B0", "size": 12},
            ],
        )
        out_html = out_dir / "escher_map.html"
        builder.save_html(str(out_html))
        log.info("  Saved: %s", out_html)
    except Exception as exc:
        log.warning("  Escher map generation failed: %s", exc)


# ---------------------------------------------------------------------------
# 7b. Escher maps for gene mutations (knockout / knockin overlays)
# ---------------------------------------------------------------------------

def export_escher_map_mutations(
    model: "cobra.Model",
    out_dir: Path,
    log: logging.Logger,
    gene_results: "list[dict] | None" = None,
    map_name: str = "iJO1366.Central metabolism",
) -> None:
    """
    Generate one Escher HTML map per mutation entry in *gene_results*,
    each showing the FBA flux distribution under that gene perturbation.

    If *gene_results* is None or empty the function is a no-op.
    The maps are written to ``<out_dir>/escher_mutations/`` as
    ``<type>_<gene_id>.html`` (e.g. ``knockout_b2215.html``).
    A companion ``escher_mutations_index.html`` lists and iframes all maps.
    """
    if not _HAS_ESCHER:
        log.warning("escher not installed – skipping mutation Escher maps.")
        return
    if not gene_results:
        return

    mut_dir = out_dir / "escher_mutations"
    mut_dir.mkdir(parents=True, exist_ok=True)

    log.info("  Building Escher maps for %d gene mutations …", len(gene_results))

    # ── wild-type fluxes (base reference) ────────────────────────────────────
    try:
        wt_sol = model.optimize()
        wt_fluxes = wt_sol.fluxes.to_dict() if wt_sol.status == "optimal" else {}
    except Exception:
        wt_fluxes = {}

    generated: list[tuple[str, str, str]] = []   # (filename, label, type)

    for entry in gene_results:
        gid   = entry.get("gene_id", "unknown")
        gname = entry.get("gene_name", gid)
        etype = entry.get("type", "knockout")
        note  = entry.get("note", "")
        rxn_affected = entry.get("reactions_affected", [])

        try:
            with model:
                if etype == "knockout":
                    gene_obj = model.genes.query(lambda g: g.id == gid)
                    if gene_obj:
                        gene_obj[0].knock_out()
                elif etype == "knockin":
                    # Open all reactions associated with this gene
                    gene_obj = model.genes.query(lambda g: g.id == gid)
                    if gene_obj:
                        for rxn in gene_obj[0].reactions:
                            if rxn.lower_bound >= 0:
                                rxn.lower_bound = 0
                                rxn.upper_bound = max(rxn.upper_bound, 1000)
                            else:
                                rxn.lower_bound = -1000
                                rxn.upper_bound = max(rxn.upper_bound, 1000)
                sol = model.optimize()
                mut_fluxes = sol.fluxes.to_dict() if sol.status == "optimal" else {}
                growth = sol.objective_value if sol.status == "optimal" else 0.0
        except Exception as exc:
            log.warning("    Skipping Escher map for %s: %s", gid, exc)
            continue

        # Compute flux *difference* relative to WT for colour scale
        diff_fluxes = {rxn_id: mut_fluxes.get(rxn_id, 0.0) - wt_fluxes.get(rxn_id, 0.0)
                       for rxn_id in set(mut_fluxes) | set(wt_fluxes)}

        # Highlight affected reactions in orange
        highlight = {r: 99999 for r in rxn_affected}

        try:
            builder = escher.Builder(
                map_name=map_name,
                model=model,
                reaction_data=mut_fluxes,
                reaction_scale=[
                    {"type": "min",  "color": "#C44E52", "size": 8},
                    {"type": "zero", "color": "#eeeeee", "size": 3},
                    {"type": "max",  "color": "#4C72B0", "size": 8},
                ],
            )
            safe_id  = gid.replace("/", "_").replace(" ", "_")
            filename = f"{etype}_{safe_id}.html"
            out_html = mut_dir / filename
            builder.save_html(str(out_html))
            label = f"{etype.capitalize()}: {gname} ({gid})  |  growth={growth:.4f}"
            if note:
                label += f"  [{note}]"
            generated.append((filename, label, etype))
            log.info("    Saved: %s", out_html)
        except Exception as exc:
            log.warning("    Escher map failed for %s: %s", gid, exc)
            continue

    if not generated:
        return

    # ── index page ────────────────────────────────────────────────────────────
    tab_btns  = ""
    tab_iframes = ""
    for i, (fn, lbl, etype) in enumerate(generated):
        active  = "active" if i == 0 else ""
        display = "block"  if i == 0 else "none"
        color   = "#C44E52" if etype == "knockout" else "#4C72B0"
        tab_btns += (
            f'<button class="mtab {active}" '
            f'style="border-left:4px solid {color};" '
            f'onclick="switchMut(this,\'mf-{i}\')">{lbl}</button>\n'
        )
        tab_iframes += (
            f'<iframe id="mf-{i}" src="{fn}" '
            f'style="display:{display};width:100%;height:85vh;border:none;"></iframe>\n'
        )

    index_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8"/>
  <title>Escher – Gene Mutation Maps</title>
  <style>
    body {{ font-family: sans-serif; margin: 0; background: #f4f6f9; }}
    header {{ background: #2c5282; color: #fff; padding: 1rem 1.5rem; }}
    header h1 {{ font-size: 1.3rem; margin: 0; }}
    .sidebar {{ width: 300px; float: left; height: calc(100vh - 60px);
                overflow-y: auto; background: #fff;
                border-right: 1px solid #e2e8f0; padding: .5rem; }}
    .main {{ margin-left: 300px; }}
    .mtab {{ display: block; width: 100%; text-align: left; padding: .55rem .75rem;
              border: none; background: none; cursor: pointer; font-size: .82rem;
              border-bottom: 1px solid #f0f0f0; border-radius: 4px; margin-bottom: 2px; }}
    .mtab:hover {{ background: #f7fafc; }}
    .mtab.active {{ background: #ebf4ff; font-weight: 600; }}
  </style>
</head>
<body>
<header><h1>Escher Gene Mutation Maps  ({len(generated)} mutations)</h1></header>
<div class="sidebar">
{tab_btns}
</div>
<div class="main">
{tab_iframes}
</div>
<script>
function switchMut(btn, id) {{
  document.querySelectorAll('.mtab').forEach(b => b.classList.remove('active'));
  document.querySelectorAll('iframe').forEach(f => f.style.display = 'none');
  btn.classList.add('active');
  document.getElementById(id).style.display = 'block';
}}
</script>
</body>
</html>"""

    index_path = mut_dir / "escher_mutations_index.html"
    index_path.write_text(index_html, encoding="utf-8")
    log.info("  Saved mutation Escher index: %s", index_path)


# ---------------------------------------------------------------------------
# 8. Machine-readable CSV summary
# ---------------------------------------------------------------------------

def export_summary_csv(
    model: "cobra.Model",
    blocked: list[str],
    fba_value: float | None,
    out_dir: Path,
    log: logging.Logger,
    mutation_df: "pd.DataFrame | None" = None,
) -> None:
    """Write a one-row CSV of key model statistics, optionally including mutation analysis KO stats."""
    if not _HAS_PD:
        log.warning("pandas not installed – skipping CSV summary.")
        return

    import pandas as pd
    from collections import Counter
    subsystems = Counter(r.subsystem for r in model.reactions if r.subsystem)

    data: dict = {
        "model_id":          model.id or "model",
        "n_reactions":       len(model.reactions),
        "n_metabolites":     len(model.metabolites),
        "n_genes":           len(model.genes),
        "n_blocked":         len(blocked),
        "pct_blocked":       round(100 * len(blocked) / max(len(model.reactions), 1), 2),
        "n_subsystems":      len(subsystems),
        "fba_objective":     round(fba_value, 6) if fba_value is not None else "N/A",
        "n_compartments":    len(model.compartments),
    }

    # Add mutation analysis statistics if available
    if mutation_df is not None and not mutation_df.empty:
        ko_df = mutation_df[mutation_df["type"] == "knockout"]
        ki_df = mutation_df[mutation_df["type"] == "knockin"]
        data["ko_genes_screened"]   = len(ko_df)
        data["ko_essential"]        = int(ko_df["essential"].sum()) if not ko_df.empty else 0
        data["ko_nonessential"]     = int((~ko_df["essential"]).sum()) if not ko_df.empty else 0
        data["ki_candidates"]       = len(ki_df)
        data["ki_growth_improving"] = int((ki_df["growth_change_pct"] > 1).sum()) if not ki_df.empty else 0

    df = pd.DataFrame([data])
    out_csv = out_dir / "model_summary.csv"
    df.to_csv(out_csv, index=False)
    log.info("  Saved: %s", out_csv)


# ---------------------------------------------------------------------------
# 9. Mutation analysis – knockout & knockin
# ---------------------------------------------------------------------------

def run_mutation_analysis(
    model: "cobra.Model",
    out_dir: Path,
    log: logging.Logger,
    max_genes: int = 200,
) -> "pd.DataFrame | None":
    """
    Perform single-gene knockout and knockin analysis and write an interactive
    HTML report (``mutation_analysis.html``).

    Knockout
    --------
    For every gene, temporarily delete it (using COBRApy's context manager)
    and record the resulting FBA growth.  Genes whose deletion drives growth
    to ≤5 % of wild-type are labelled **essential**; genes with no effect are
    **non-essential**.

    Knockin
    -------
    Genes whose flux is zero in the wild-type FBA solution but become active
    after releasing their bound constraints are candidates for "knock-in"
    rescue experiments.  We identify these by checking which gene-associated
    reactions carry zero flux in the WT solution, then simulating each one
    individually with its lower-bound opened (set to –1000 if reversible, or
    0→1000 if irreversible) and measuring the change in growth.

    Parameters
    ----------
    model    : COBRApy Model
    out_dir  : output directory
    log      : logger
    max_genes: cap on the number of genes screened (for speed; sorted by
               number of associated reactions descending so the most-connected
               genes are always included).

    Returns the results DataFrame (or None on failure).
    """
    if not _HAS_PLOTLY or not _HAS_PD:
        log.warning("plotly/pandas not installed – skipping mutation analysis.")
        return None

    import pandas as pd

    log.info("  Running mutation analysis (knockout + knockin) …")

    # ── wild-type FBA ────────────────────────────────────────────────────────
    try:
        wt_solution = model.optimize()
        if wt_solution.status != "optimal":
            log.warning("  WT FBA not optimal – skipping mutation analysis.")
            return None
        wt_growth = wt_solution.objective_value
        if wt_growth <= 0:
            log.warning("  WT growth ≤ 0 – skipping mutation analysis.")
            return None
    except Exception as exc:
        log.warning("  WT FBA failed: %s", exc)
        return None

    wt_fluxes = wt_solution.fluxes  # Series: rxn_id → flux

    # ── gene list (capped) ───────────────────────────────────────────────────
    genes = sorted(model.genes, key=lambda g: -len(g.reactions))
    if len(genes) > max_genes:
        log.info("    Capping gene screen at %d / %d genes.", max_genes, len(genes))
        genes = genes[:max_genes]

    # ── knockout screen ──────────────────────────────────────────────────────
    ko_rows: list[dict] = []
    for gene in genes:
        reaction_ids = sorted(r.id for r in gene.reactions)
        subsystem_votes = [
            (r.subsystem or "").strip()
            for r in gene.reactions
            if (r.subsystem or "").strip()
        ]
        dominant_subsystem = (
            Counter(subsystem_votes).most_common(1)[0][0]
            if subsystem_votes
            else "Unknown subsystem"
        )
        gene_product = gene.name or gene.id
        try:
            with model:
                gene.knock_out()
                sol = model.optimize()
                ko_growth = sol.objective_value if sol.status == "optimal" else 0.0
        except Exception:
            ko_growth = float("nan")
        ratio = ko_growth / wt_growth if wt_growth else float("nan")
        ko_rows.append({
            "gene_id":    gene.id,
            "gene_name":  gene.name or gene.id,
            "gene_product": gene_product,
            "n_reactions": len(gene.reactions),
            "subsystem": dominant_subsystem,
            "reaction_ids": "|".join(reaction_ids),
            "type":       "knockout",
            "wt_growth":  round(wt_growth, 6),
            "mut_growth": round(ko_growth, 6),
            "growth_ratio": round(ratio, 4),
            "essential":  ratio <= 0.05,
            "growth_change_pct": round((ratio - 1) * 100, 2),
        })

    # ── knockin screen ───────────────────────────────────────────────────────
    # Identify genes whose reactions all carry zero WT flux → candidates
    ki_rows: list[dict] = []
    ki_genes = [
        g for g in genes
        if all(abs(wt_fluxes.get(r.id, 0.0)) < 1e-9 for r in g.reactions)
        and len(g.reactions) > 0
    ]
    log.info("    %d knockout genes screened; %d knockin candidates.", len(genes), len(ki_genes))

    for gene in ki_genes[:max_genes]:
        reaction_ids = sorted(r.id for r in gene.reactions)
        subsystem_votes = [
            (r.subsystem or "").strip()
            for r in gene.reactions
            if (r.subsystem or "").strip()
        ]
        dominant_subsystem = (
            Counter(subsystem_votes).most_common(1)[0][0]
            if subsystem_votes
            else "Unknown subsystem"
        )
        gene_product = gene.name or gene.id
        try:
            with model:
                for rxn in gene.reactions:
                    # Open the reaction: allow flux to pass through
                    if rxn.lower_bound >= 0:
                        rxn.lower_bound = 0
                        rxn.upper_bound = max(rxn.upper_bound, 1000)
                    else:
                        rxn.lower_bound = -1000
                        rxn.upper_bound = max(rxn.upper_bound, 1000)
                sol = model.optimize()
                ki_growth = sol.objective_value if sol.status == "optimal" else 0.0
        except Exception:
            ki_growth = float("nan")
        ratio = ki_growth / wt_growth if wt_growth else float("nan")
        ki_rows.append({
            "gene_id":    gene.id,
            "gene_name":  gene.name or gene.id,
            "gene_product": gene_product,
            "n_reactions": len(gene.reactions),
            "subsystem": dominant_subsystem,
            "reaction_ids": "|".join(reaction_ids),
            "type":       "knockin",
            "wt_growth":  round(wt_growth, 6),
            "mut_growth": round(ki_growth, 6),
            "growth_ratio": round(ratio, 4),
            "essential":  False,
            "growth_change_pct": round((ratio - 1) * 100, 2),
        })

    df = pd.DataFrame(ko_rows + ki_rows)
    if df.empty:
        log.warning("  No mutation results produced.")
        return None

    # ── Plotly figure ────────────────────────────────────────────────────────
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import random

    ko_df = df[df["type"] == "knockout"].copy()
    ki_df = df[df["type"] == "knockin"].copy()

    ko_df["subsystem"] = ko_df["subsystem"].fillna("Unknown subsystem")
    ko_df.loc[ko_df["subsystem"].astype(str).str.strip() == "", "subsystem"] = "Unknown subsystem"

    n_essential = int(ko_df["essential"].sum())
    n_nonessential = int((~ko_df["essential"]).sum())

    subsystem_stats = (
        ko_df.groupby("subsystem", as_index=False)
        .agg(
            n_genes=("gene_id", "count"),
            n_essential=("essential", "sum"),
            median_growth=("growth_ratio", "median"),
        )
    )
    subsystem_stats["essential_pct"] = (
        100.0 * subsystem_stats["n_essential"] / subsystem_stats["n_genes"].clip(lower=1)
    )
    subsystem_stats = subsystem_stats.sort_values(
        ["essential_pct", "median_growth", "n_genes"],
        ascending=[False, True, False],
    )

    subsystem_order = subsystem_stats["subsystem"].tolist()
    sub_to_idx = {s: i for i, s in enumerate(subsystem_order)}

    ko_plot = ko_df.sort_values(["subsystem", "growth_ratio", "gene_id"]).copy()
    rng = random.Random(42)
    ko_plot["x_num"] = [
        float(sub_to_idx.get(sub, 0)) + rng.uniform(-0.26, 0.26)
        for sub in ko_plot["subsystem"].tolist()
    ]
    ko_plot["point_color"] = ko_plot["essential"].map({True: "#dc2626", False: "#059669"})
    ko_plot["point_symbol"] = ko_plot["essential"].map({True: "diamond", False: "circle"})

    has_ki = not ki_df.empty
    n_rows = 3 if has_ki else 2
    titles = [
        "Subsystem Vulnerability Plot (KO growth ratio by metabolic pathway)",
        "Subsystem essentiality burden (% essential knockout genes)",
    ]
    if has_ki:
        titles.append("Knockin: growth change per candidate gene")
    specs = [[{"type": "scatter"}], [{"type": "bar"}]] + ([[{"type": "scatter"}]] if has_ki else [])

    fig = make_subplots(
        rows=n_rows,
        cols=1,
        subplot_titles=titles,
        specs=specs,
        vertical_spacing=0.09,
        row_heights=[0.52, 0.22, 0.26][:n_rows],
    )

    ko_custom = [
        [
            str(row["gene_id"]),
            str(row["gene_name"]),
            str(row.get("gene_product", row["gene_name"])),
            str(row["subsystem"]),
            float(row["growth_ratio"]),
            str(row.get("reaction_ids", "")),
        ]
        for _, row in ko_plot.iterrows()
    ]
    fig.add_trace(
        go.Scatter(
            x=ko_plot["x_num"].tolist(),
            y=ko_plot["growth_ratio"].tolist(),
            mode="markers",
            marker=dict(
                color=ko_plot["point_color"].tolist(),
                symbol=ko_plot["point_symbol"].tolist(),
                size=8,
                opacity=0.85,
                line=dict(width=0.5, color="#ffffff"),
            ),
            customdata=ko_custom,
            hovertemplate=(
                "<b>%{customdata[1]}</b><br>"
                "Gene ID: %{customdata[0]}<br>"
                "Gene product: %{customdata[2]}<br>"
                "Subsystem: %{customdata[3]}<br>"
                "Growth ratio: %{customdata[4]:.4f}<extra></extra>"
            ),
            name="Knockout genes",
            showlegend=False,
        ),
        row=1,
        col=1,
    )

    for _, row in subsystem_stats.iterrows():
        if row["n_genes"] >= 3 and row["essential_pct"] >= 40.0:
            idx = sub_to_idx.get(str(row["subsystem"]), 0)
            fig.add_vrect(
                x0=idx - 0.5,
                x1=idx + 0.5,
                fillcolor="rgba(245, 101, 101, 0.14)",
                line_width=0,
                layer="below",
                row=1,
                col=1,
            )

    fig.add_hline(
        y=0.05,
        line_dash="dash",
        line_color="#dc2626",
        annotation_text="Essentiality threshold (growth < 0.05)",
        row=1,
        col=1,
    )
    fig.add_hline(y=1.0, line_dash="dot", line_color="#64748b", row=1, col=1)

    sub_bar_colors = [
        "#ef4444" if pct >= 40 else ("#f59e0b" if pct >= 20 else "#10b981")
        for pct in subsystem_stats["essential_pct"].tolist()
    ]
    fig.add_trace(
        go.Bar(
            x=subsystem_stats["subsystem"].tolist(),
            y=subsystem_stats["essential_pct"].tolist(),
            marker_color=sub_bar_colors,
            text=[
                f"{pct:.1f}% ({int(ne)}/{int(ng)})"
                for pct, ne, ng in zip(
                    subsystem_stats["essential_pct"],
                    subsystem_stats["n_essential"],
                    subsystem_stats["n_genes"],
                )
            ],
            textposition="outside",
            hovertemplate=(
                "<b>%{x}</b><br>"
                "Essential genes: %{text}<br>"
                "Essentiality burden: %{y:.2f}%<extra></extra>"
            ),
            name="Essentiality burden",
            showlegend=False,
        ),
        row=2,
        col=1,
    )

    if has_ki:
        ki_df_sorted = ki_df.sort_values("growth_change_pct", ascending=False)
        fig.add_trace(
            go.Scatter(
                x=list(range(len(ki_df_sorted))),
                y=ki_df_sorted["growth_change_pct"].tolist(),
                mode="markers",
                marker=dict(
                    color=ki_df_sorted["growth_change_pct"].apply(
                        lambda v: "#2563eb" if v > 1 else "#d97706"
                    ).tolist(),
                    size=6,
                    opacity=0.82,
                    line=dict(width=0.4, color="white"),
                ),
                hovertemplate=(
                    "<b>%{customdata[1]}</b><br>"
                    "Gene ID: %{customdata[0]}<br>"
                    "Subsystem: %{customdata[2]}<br>"
                    "Growth change: %{y:+.2f}%<extra></extra>"
                ),
                customdata=[
                    [
                        str(row["gene_id"]),
                        str(row["gene_name"]),
                        str(row.get("subsystem", "Unknown subsystem")),
                    ]
                    for _, row in ki_df_sorted.iterrows()
                ],
                name="Knockin",
                showlegend=False,
            ),
            row=3,
            col=1,
        )
        fig.add_hline(y=0, line_dash="dot", line_color="#64748b", row=3, col=1)

    fig.update_layout(
        title_text=(
            f"Mutation Analysis — Subsystem Vulnerability ({model.id or 'model'})<br>"
            f"<sup>WT growth = {wt_growth:.4f}  |  "
            f"Genes screened = {len(genes)}  |  "
            f"Essential KO = {n_essential}</sup>"
        ),
        title_font_size=15,
        height=280 * n_rows + 200,
        paper_bgcolor="white",
        plot_bgcolor="white",
        showlegend=False,
        margin=dict(l=80, r=30, t=95, b=40),
    )
    fig.update_xaxes(
        title_text="Metabolic subsystem / pathway",
        row=1,
        col=1,
        tickmode="array",
        tickvals=list(range(len(subsystem_order))),
        ticktext=[s if len(s) <= 28 else s[:25] + "…" for s in subsystem_order],
        tickangle=35,
        range=[-0.7, max(len(subsystem_order) - 0.3, 0.7)],
    )
    fig.update_yaxes(title_text="Growth ratio (KO / WT)", row=1, col=1)
    fig.update_xaxes(title_text="Subsystem", row=2, col=1, tickangle=35)
    fig.update_yaxes(title_text="Essential genes (%)", row=2, col=1, range=[0, 100])
    if has_ki:
        fig.update_xaxes(title_text="Knockin candidate rank", row=3, col=1)
        fig.update_yaxes(title_text="Growth change (%)", row=3, col=1)

    mutation_post_script = """
    (function() {
      const gd = document.getElementsByClassName('plotly-graph-div')[0];
      if (!gd) return;
      gd.on('plotly_click', function(ev) {
        if (!ev || !ev.points || !ev.points.length) return;
        const pt = ev.points[0];
        const cd = pt.customdata;
        if (!cd || !Array.isArray(cd) || cd.length < 6) return;
        const rxnRaw = (cd[5] || '').toString();
        const reactionIds = rxnRaw ? rxnRaw.split('|').filter(Boolean) : [];
        if (!reactionIds.length) return;
        window.parent.postMessage({
          type: 'gsmm-highlight-reactions',
          reactionIds: reactionIds,
          geneId: cd[0],
          geneName: cd[1],
          subsystem: cd[3]
        }, '*');
      });
    })();
    """

    out_html = out_dir / "mutation_analysis.html"
    fig.write_html(
        str(out_html),
        include_plotlyjs="cdn",
        config={"displayModeBar": True, "scrollZoom": True, "displaylogo": False},
        post_script=mutation_post_script,
    )
    log.info("  Saved: %s  (%d KO, %d KI records)", out_html, len(ko_df), len(ki_df))

    return df


# ---------------------------------------------------------------------------
# Master entry point
# ---------------------------------------------------------------------------

def run_all_visualizations(
    model: "cobra.Model",
    blocked: list[str] | None = None,
    fba_value: float | None = None,
    out_dir: Path | None = None,
    log: logging.Logger | None = None,
) -> None:
    """
    Run every visualization function and write outputs to ``out_dir``.

    Parameters
    ----------
    model     : COBRApy Model object
    blocked   : list of blocked reaction IDs (pass result of
                ``find_blocked_reactions``).  Computed here if None.
    fba_value : FBA objective value.  Re-computed here if None.
    out_dir   : destination directory (created if absent)
    log       : a Python logger; a basic one is created if None
    """
    if log is None:
        logging.basicConfig(level=logging.INFO,
                            format="%(asctime)s [%(levelname)s] %(message)s",
                            datefmt="%H:%M:%S")
        log = logging.getLogger("visualize_model")

    if out_dir is None:
        out_dir = Path("output/viz")
    out_dir.mkdir(parents=True, exist_ok=True)

    log.info("Generating visualizations → %s", out_dir)

    # Ensure we have blocked / FBA data
    if blocked is None:
        try:
            from cobra.flux_analysis import find_blocked_reactions
            blocked = find_blocked_reactions(model)
        except Exception:
            blocked = []

    if fba_value is None:
        try:
            sol = model.optimize()
            fba_value = sol.objective_value if sol.status == "optimal" else None
        except Exception:
            fba_value = None

    # Static matplotlib/seaborn plots
    plot_summary_stats(model, blocked, fba_value, out_dir, log)
    plot_compartment_breakdown(model, out_dir, log)
    plot_degree_distribution(model, out_dir, log)
    plot_subsystem_barchart(model, out_dir, log)
    plot_flux_distribution(model, out_dir, log)

    # Interactive Plotly outputs
    plot_reaction_network(model, out_dir, log)
    plot_dashboard(model, blocked, fba_value, out_dir, log)

    # Mutation analysis (knockout + knockin)
    mutation_df = run_mutation_analysis(model, out_dir, log)

    # Escher map (wild-type)
    export_escher_map(model, out_dir, log)

    # Escher maps for gene mutations from gene_config.json (if present)
    # Search order: script dir (project root), then relative to out_dir hierarchy
    _script_dir = Path(__file__).parent
    gene_config_path = _script_dir / "gene_config.json"
    if not gene_config_path.is_file():
        gene_config_path = _script_dir.parent / "config" / "gene_config.json"
    if not gene_config_path.is_file():
        gene_config_path = out_dir.parent.parent / "gene_config.json"
    if not gene_config_path.is_file():
        gene_config_path = out_dir.parent.parent / "config" / "gene_config.json"
    if not gene_config_path.is_file():
        gene_config_path = Path("gene_config.json")
    if not gene_config_path.is_file():
        gene_config_path = Path("config") / "gene_config.json"
    if gene_config_path.is_file():
        try:
            from gene_config import load_gene_config, apply_gene_config
            cfg = load_gene_config(gene_config_path)
            gene_results = apply_gene_config(model, cfg, log=log)
            export_escher_map_mutations(model, out_dir, log, gene_results=gene_results)
        except Exception as exc:
            log.warning("  Gene config Escher maps skipped: %s", exc)

    # CSV summary (pass mutation_df so KO stats can be included)
    export_summary_csv(model, blocked, fba_value, out_dir, log, mutation_df=mutation_df)

    # ---------------------------------------------------------------------------
    # Cross-feeding network (feature [3] – MAG metagenomics inputs)
    # Auto-discover a MAG annotation table in the standard search paths:
    #   <project_root>/mag_annotations.csv  (or .tsv)
    #   <out_dir>/../../mag_annotations.csv
    # If found, run the cross-feeding network analysis.
    # ---------------------------------------------------------------------------
    _script_dir = Path(__file__).parent
    _mag_candidates = [
        _script_dir / "mag_annotations.csv",
        _script_dir / "mag_annotations.tsv",
        out_dir.parent.parent / "mag_annotations.csv",
        out_dir.parent.parent / "mag_annotations.tsv",
        Path("mag_annotations.csv"),
        Path("mag_annotations.tsv"),
    ]
    _mag_table = next((p for p in _mag_candidates if p.is_file()), None)
    if _mag_table is not None:
        try:
            from crossfeed_network import run_crossfeed_analysis
            crossfeed_dir = out_dir / "crossfeed"
            log.info("  Found MAG annotation table: %s", _mag_table)
            run_crossfeed_analysis(
                table_path=_mag_table,
                out_dir=crossfeed_dir,
                log=log,
            )
        except Exception as exc:
            log.warning("  Cross-feeding network skipped: %s", exc)
    else:
        log.info("  No mag_annotations.csv found – cross-feeding network skipped.")

    log.info("All visualizations complete.  Output directory: %s", out_dir)

    # Generate the single consolidated index.html report
    generate_index_html(model, blocked, fba_value, out_dir, log, mutation_df=mutation_df)


# ---------------------------------------------------------------------------
# 10. index.html – single consolidated report with environment editor
# ---------------------------------------------------------------------------

def generate_index_html(
    model: "cobra.Model",
    blocked: list[str],
    fba_value: float | None,
    out_dir: Path,
    log: logging.Logger,
    mutation_df: "pd.DataFrame | None" = None,
) -> None:
    """
    Generate a self-contained ``index.html`` that acts as a unified report hub:

    - **Overview tab** – key model statistics, inline PNG figures (base64),
      and a collapsible mutation analysis summary table.
    - **Interactive reports tab** – iframes embedding ``dashboard.html``,
      ``reaction_network.html``, ``mutation_analysis.html``, and
      ``escher_map.html`` side by side.
    - **Environment editor tab** – an editable form (feature [3]) that lets
      the user inspect and modify the growth medium (gap-fill medium name,
      exchange reaction bounds) and shows a shell command to re-run the
      pipeline with the chosen settings.

    The HTML file is written to ``<out_dir>/index.html``.
    """
    import base64
    import json
    from datetime import datetime

    log.info("  Generating index.html report …")

    # ── helper: embed a PNG as a base64 data URI ──────────────────────────────
    def _b64_img(name: str) -> str:
        p = out_dir / name
        if not p.is_file():
            return ""
        with open(p, "rb") as fh:
            b64 = base64.b64encode(fh.read()).decode()
        return f"data:image/png;base64,{b64}"

    # ── collect available files ───────────────────────────────────────────────
    html_reports = [
        ("dashboard.html",         "Dashboard"),
        ("reaction_network.html",  "Reaction Network"),
        ("mutation_analysis.html", "Mutation Analysis"),
        ("escher_map.html",        "Escher Map (WT)"),
        ("escher_mutations/escher_mutations_index.html", "Escher Mutations"),
        ("gene_config_report.html", "Gene Config Report"),
        ("crossfeed/crossfeed_network.html", "Cross-Feeding Network"),
    ]
    available_html = [(fn, lbl) for fn, lbl in html_reports if (out_dir / fn).is_file()]

    png_files = [
        ("summary_stats.png",        "Summary Statistics"),
        ("compartment_breakdown.png","Compartment Breakdown"),
        ("degree_distribution.png",  "Degree Distribution"),
        ("subsystem_barchart.png",   "Subsystem Barchart"),
        ("flux_distribution.png",    "FBA Flux Distribution"),
    ]
    available_pngs = [(fn, lbl, _b64_img(fn)) for fn, lbl, in png_files if _b64_img(fn)]

    # ── model statistics ──────────────────────────────────────────────────────
    n_rxn   = len(model.reactions)
    n_met   = len(model.metabolites)
    n_gene  = len(model.genes)
    n_blk   = len(blocked) if blocked else 0
    fba_str = f"{fba_value:.4f}" if fba_value is not None else "N/A"
    model_id = model.id or "GSMM model"
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # ── mutation summary table ─────────────────────────────────────────────────
    mut_table_rows = ""
    n_essential = n_nonessential = ki_candidates = 0
    if mutation_df is not None and not mutation_df.empty:
        import pandas as pd
        ko_df = mutation_df[mutation_df["type"] == "knockout"]
        ki_df = mutation_df[mutation_df["type"] == "knockin"]
        n_essential    = int(ko_df["essential"].sum()) if not ko_df.empty else 0
        n_nonessential = int((~ko_df["essential"]).sum()) if not ko_df.empty else 0
        ki_candidates  = len(ki_df)

        top_rows = mutation_df.sort_values("growth_ratio").head(20)
        for _, row in top_rows.iterrows():
            badge = ""
            if row["type"] == "knockout":
                badge = ('<span class="badge badge-red">Essential</span>'
                         if row["essential"]
                         else '<span class="badge badge-green">Non-essential</span>')
            else:
                badge = '<span class="badge badge-blue">Knockin</span>'
            mut_table_rows += (
                f"<tr>"
                f"<td>{row['gene_id']}</td>"
                f"<td>{row['gene_name']}</td>"
                f"<td>{row['type'].capitalize()}</td>"
                f"<td>{row['n_reactions']}</td>"
                f"<td>{row['growth_ratio']:.4f}</td>"
                f"<td>{row['growth_change_pct']:+.2f}%</td>"
                f"<td>{badge}</td>"
                f"</tr>\n"
            )

    # ── medium / exchange reactions ────────────────────────────────────────────
    # Gather exchange reactions and their bounds for the environment editor
    exchange_rxns = [r for r in model.reactions if r.id.startswith("EX_")]
    # Show only those that are open (lb < 0 → can be consumed)
    open_exchanges = [(r.id, r.name or r.id, r.lower_bound, r.upper_bound)
                      for r in exchange_rxns if r.lower_bound < 0]
    open_exchanges.sort(key=lambda x: x[1])
    env_rows_html = ""
    for rxn_id, rxn_name, lb, ub in open_exchanges[:60]:  # cap at 60 for readability
        env_rows_html += (
            f'<tr data-rxn="{rxn_id}">'
            f'<td class="rxn-id">{rxn_id}</td>'
            f'<td>{rxn_name}</td>'
            f'<td><input type="number" class="lb-input" value="{lb}" step="0.1" '
            f'     data-orig="{lb}" title="Lower bound (negative = uptake allowed)" /></td>'
            f'<td><input type="number" class="ub-input" value="{ub}" step="0.1" '
            f'     data-orig="{ub}" title="Upper bound" /></td>'
            f'</tr>\n'
        )

    # JSON of original bounds for reset
    orig_bounds_json = json.dumps(
        {r.id: {"lb": r.lower_bound, "ub": r.upper_bound} for r in exchange_rxns
         if r.lower_bound < 0}
    )

    # ── load environment presets from environments.json ───────────────────────
    env_presets: dict = {}
    # Look for environments.json relative to the model (up two levels from viz dir)
    for candidate in [out_dir.parent.parent / "environments.json",
                      out_dir.parent.parent / "config" / "environments.json",
                      out_dir.parent / "environments.json",
                      out_dir.parent / "config" / "environments.json",
                      Path("config") / "environments.json",
                      Path("environments.json")]:
        if candidate.is_file():
            try:
                with open(candidate) as _fh:
                    _env_data = json.load(_fh)
                env_presets = _env_data.get("presets", {})
            except Exception:
                pass
            break
    env_presets_json = json.dumps(env_presets)

    # ── load cross-feeding network summary (if present) ───────────────────────
    crossfeed_summary: dict = {}
    crossfeed_html_path = out_dir / "crossfeed" / "crossfeed_network.html"
    crossfeed_json_path = out_dir / "crossfeed" / "crossfeed_summary.json"
    if crossfeed_json_path.is_file():
        try:
            with open(crossfeed_json_path, encoding="utf-8") as _cfh:
                crossfeed_summary = json.load(_cfh)
        except Exception:
            pass
    _has_crossfeed = crossfeed_html_path.is_file()

    # Build preset button HTML
    preset_btns_html = ""
    for key, preset in env_presets.items():
        color = preset.get("color", "#4a5568")
        label = preset.get("label", key)
        desc  = preset.get("description", "")
        preset_btns_html += (
            f'<button class="preset-btn" '
            f'style="border-left:4px solid {color};" '
            f'onclick="applyPreset({json.dumps(key)})" '
            f'title="{desc}">{label}</button>\n'
        )

    # ── tab buttons for interactive reports ───────────────────────────────────
    iframe_tabs = ""
    iframe_panels = ""
    for i, (fn, lbl) in enumerate(available_html):
        active = "active" if i == 0 else ""
        iframe_tabs += (
            f'<button id="iframe-tab-{i}" class="tab-btn {active}" '
            f'data-target="iframe-{i}" data-src="{fn}" '
            f'onclick="switchIframe(this,\'iframe-{i}\')">'
            f'{lbl}</button>\n'
        )
        display = "block" if i == 0 else "none"
        iframe_panels += (
            f'<iframe id="iframe-{i}" src="{fn}" style="display:{display};'
            f'width:100%;height:80vh;border:none;"></iframe>\n'
        )

    # ── PNG gallery HTML ──────────────────────────────────────────────────────
    png_gallery = ""
    for fn, lbl, b64 in available_pngs:
        if b64:
            png_gallery += (
                f'<div class="png-card">'
                f'<p class="png-label">{lbl}</p>'
                f'<img src="{b64}" alt="{lbl}" loading="lazy" />'
                f'</div>\n'
            )

    # ── mutation summary stats cards ──────────────────────────────────────────
    mut_stat_cards = ""
    if mutation_df is not None:
        mut_stat_cards = f"""
        <div class="stat-cards">
          <div class="stat-card red">
            <div class="stat-value">{n_essential}</div>
            <div class="stat-label">Essential KO genes</div>
          </div>
          <div class="stat-card green">
            <div class="stat-value">{n_nonessential}</div>
            <div class="stat-label">Non-essential KO genes</div>
          </div>
          <div class="stat-card blue">
            <div class="stat-value">{ki_candidates}</div>
            <div class="stat-label">Knockin candidates</div>
          </div>
        </div>"""

    # ── cross-feeding summary HTML ────────────────────────────────────────────
    crossfeed_tab_btn = ""
    crossfeed_panel_html = ""
    if _has_crossfeed:
        cf_n_mags  = crossfeed_summary.get("n_mags", "—")
        cf_n_edges = crossfeed_summary.get("n_edges", "—")
        cf_n_mets  = crossfeed_summary.get("n_metabolites", "—")
        cf_mag_ids = crossfeed_summary.get("mag_ids", [])
        cf_top_mets = crossfeed_summary.get("top_metabolites", [])

        # Top metabolite table rows
        cf_top_met_rows = ""
        for item in cf_top_mets[:15]:
            cf_top_met_rows += (
                f"<tr>"
                f"<td>{item.get('compound_id','')}</td>"
                f"<td>{item.get('compound_name','')}</td>"
                f"<td>{item.get('n_edges','')}</td>"
                f"</tr>\n"
            )

        # MAG stats rows
        cf_mag_stat_rows = ""
        cf_mag_stats = crossfeed_summary.get("mag_stats", {})
        for mag_id in cf_mag_ids:
            stats = cf_mag_stats.get(mag_id, {})
            cf_mag_stat_rows += (
                f"<tr>"
                f"<td>{mag_id}</td>"
                f"<td>{stats.get('n_products','—')}</td>"
                f"<td>{stats.get('n_substrates','—')}</td>"
                f"</tr>\n"
            )

        crossfeed_tab_btn = (
            '<button class="main-tab" '
            'onclick="switchPanel(this,\'panel-crossfeed\')">Cross-Feeding Network</button>'
        )
        crossfeed_panel_html = f"""
<!-- ══════════════════════════════════════════════════════════
     PANEL 4 – Cross-Feeding Network  (Feature [1] feature_3)
     ══════════════════════════════════════════════════════════ -->
<div id="panel-crossfeed" class="panel">
  <p style="margin-bottom:1rem;font-size:.9rem;color:#4a5568;">
    Metabolic hand-off &amp; cross-feeding network for metagenomics MAG inputs.
    Nodes represent MAGs (left) and exchanged metabolites (right); directed edges
    show which MAG produces a metabolite and which MAG consumes it.
  </p>

  <!-- Summary stat cards -->
  <div class="stat-row" style="margin-bottom:1.2rem;">
    <div class="stat-box"><div class="val">{cf_n_mags}</div><div class="lbl">MAGs</div></div>
    <div class="stat-box"><div class="val">{cf_n_mets}</div><div class="lbl">Exchanged metabolites</div></div>
    <div class="stat-box"><div class="val">{cf_n_edges}</div><div class="lbl">Cross-feeding edges</div></div>
  </div>

  <!-- Interactive bipartite graph iframe -->
  <h2 style="margin-bottom:.5rem;">Interactive Network</h2>
  <iframe src="crossfeed/crossfeed_network.html"
          style="width:100%;height:75vh;border:1px solid #e2e8f0;border-radius:8px;">
  </iframe>

  <!-- MAG statistics table -->
  <h2 style="margin-top:1.5rem;margin-bottom:.5rem;">MAG Metabolic Capacity</h2>
  <div style="overflow-x:auto;margin-bottom:1.5rem;">
  <table class="data-table">
    <thead><tr><th>MAG ID</th><th>Products</th><th>Substrates</th></tr></thead>
    <tbody>{cf_mag_stat_rows}</tbody>
  </table>
  </div>

  <!-- Top exchanged metabolites table -->
  {'<h2 style="margin-bottom:.5rem;">Top Exchanged Metabolites</h2>' if cf_top_met_rows else ''}
  {f'<div style="overflow-x:auto;"><table class="data-table"><thead><tr><th>Compound ID</th><th>Name</th><th>Cross-feeding edges</th></tr></thead><tbody>{cf_top_met_rows}</tbody></table></div>' if cf_top_met_rows else ''}

  <p style="margin-top:1rem;font-size:.82rem;color:#718096;">
    Source: <code>crossfeed/crossfeed_network.html</code> &nbsp;|&nbsp;
    Edge list: <code>crossfeed/crossfeed_edges.csv</code> &nbsp;|&nbsp;
    Summary: <code>crossfeed/crossfeed_summary.json</code>
  </p>
</div>
"""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>GSMM Report – {model_id}</title>
  <style>
    *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
      background: #f4f6f9;
      color: #2d3748;
      line-height: 1.5;
    }}
    header {{
      background: linear-gradient(135deg, #2c5282, #4c72b0);
      color: #fff;
      padding: 1.5rem 2rem;
    }}
    header h1 {{ font-size: 1.6rem; font-weight: 700; }}
    header p  {{ opacity: 0.85; font-size: 0.9rem; margin-top: 0.25rem; }}

    /* ── main tabs ── */
    .main-tabs {{
      display: flex;
      gap: 0;
      background: #fff;
      border-bottom: 2px solid #e2e8f0;
      position: sticky; top: 0; z-index: 100;
      box-shadow: 0 2px 6px rgba(0,0,0,.08);
    }}
    .main-tab {{
      padding: 0.85rem 1.6rem;
      cursor: pointer;
      font-weight: 600;
      font-size: 0.9rem;
      border: none;
      background: none;
      color: #718096;
      border-bottom: 3px solid transparent;
      margin-bottom: -2px;
      transition: color .15s, border-color .15s;
    }}
    .main-tab:hover {{ color: #4c72b0; }}
    .main-tab.active {{ color: #2c5282; border-bottom-color: #4c72b0; }}

    /* ── tab panels ── */
    .panel {{ display: none; padding: 1.5rem 2rem; }}
    .panel.active {{ display: block; }}

    /* ── stat bar ── */
    .stat-row {{
      display: flex; flex-wrap: wrap; gap: 1rem; margin-bottom: 1.5rem;
    }}
    .stat-box {{
      background: #fff; border-radius: 10px; padding: 1rem 1.5rem;
      flex: 1; min-width: 130px;
      box-shadow: 0 1px 4px rgba(0,0,0,.07);
      text-align: center;
    }}
    .stat-box .val {{ font-size: 1.8rem; font-weight: 700; color: #2c5282; }}
    .stat-box .lbl {{ font-size: 0.75rem; color: #718096; text-transform: uppercase;
                      letter-spacing: .05em; margin-top: .25rem; }}

    /* ── section headings ── */
    h2 {{ font-size: 1.1rem; font-weight: 700; color: #2d3748;
          margin: 1.5rem 0 0.75rem; border-left: 4px solid #4c72b0;
          padding-left: 0.6rem; }}

    /* ── PNG gallery ── */
    .png-gallery {{
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(380px, 1fr));
      gap: 1rem;
    }}
    .png-card {{
      background: #fff; border-radius: 10px; padding: 1rem;
      box-shadow: 0 1px 4px rgba(0,0,0,.07);
    }}
    .png-card img {{ width: 100%; height: auto; border-radius: 6px; display: block; }}
    .png-label {{ font-size: 0.8rem; font-weight: 600; color: #4a5568;
                  margin-bottom: .5rem; text-transform: uppercase; letter-spacing:.04em; }}

    /* ── mutation section ── */
    .stat-cards {{
      display: flex; flex-wrap: wrap; gap: 1rem; margin: 0.75rem 0 1rem;
    }}
    .stat-card {{
      flex: 1; min-width: 140px; border-radius: 10px; padding: .9rem 1.2rem;
      text-align: center; color: #fff;
    }}
    .stat-card.red   {{ background: #e53e3e; }}
    .stat-card.green {{ background: #38a169; }}
    .stat-card.blue  {{ background: #4c72b0; }}
    .stat-value {{ font-size: 2rem; font-weight: 700; }}
    .stat-label {{ font-size: 0.75rem; opacity: 0.9; margin-top: .2rem; }}

    /* ── badges ── */
    .badge {{
      display: inline-block; padding: .15em .55em;
      border-radius: 4px; font-size: 0.72rem; font-weight: 600;
    }}
    .badge-red   {{ background:#fed7d7; color:#c53030; }}
    .badge-green {{ background:#c6f6d5; color:#276749; }}
    .badge-blue  {{ background:#bee3f8; color:#2a69ac; }}

    /* ── tables ── */
    .data-table {{ width: 100%; border-collapse: collapse; font-size: 0.85rem;
                   background: #fff; border-radius: 8px; overflow: hidden;
                   box-shadow: 0 1px 4px rgba(0,0,0,.07); }}
    .data-table th {{
      background: #2c5282; color: #fff; padding: .6rem .8rem;
      text-align: left; font-weight: 600; font-size: 0.78rem;
    }}
    .data-table td {{ padding: .5rem .8rem; border-bottom: 1px solid #edf2f7; }}
    .data-table tr:last-child td {{ border-bottom: none; }}
    .data-table tr:hover td {{ background: #f7fafc; }}

    /* ── iframe tab buttons ── */
    .tab-btn {{
      padding: .5rem 1.1rem; border: 1px solid #e2e8f0; background: #fff;
      cursor: pointer; font-size: .85rem; font-weight: 500; border-radius: 6px;
      color: #4a5568; margin-right: .4rem; margin-bottom: .75rem;
      transition: background .12s, color .12s;
    }}
    .tab-btn.active {{ background: #4c72b0; color: #fff; border-color: #4c72b0; }}
    .tab-btn:hover:not(.active) {{ background: #edf2f7; }}

    /* ── environment editor ── */
    .env-toolbar {{
      display: flex; gap: .75rem; flex-wrap: wrap; align-items: center;
      margin-bottom: 1rem;
    }}
    .env-toolbar input[type=text] {{
      padding: .45rem .75rem; border: 1px solid #cbd5e0; border-radius: 6px;
      font-size: .9rem; width: 220px;
    }}
    .btn {{
      padding: .45rem 1rem; border: none; border-radius: 6px;
      font-size: .85rem; font-weight: 600; cursor: pointer;
    }}
    .btn-primary {{ background: #4c72b0; color: #fff; }}
    .btn-primary:hover {{ background: #2c5282; }}
    .btn-secondary {{ background: #e2e8f0; color: #2d3748; }}
    .btn-secondary:hover {{ background: #cbd5e0; }}
    .btn-danger {{ background: #fed7d7; color: #c53030; }}
    .btn-danger:hover {{ background: #feb2b2; }}
    .lb-input, .ub-input {{
      width: 80px; padding: .3rem .4rem; border: 1px solid #cbd5e0;
      border-radius: 4px; font-size: .85rem; text-align: right;
    }}
    .lb-input.changed, .ub-input.changed {{
      border-color: #f6ad55; background: #fffaf0;
    }}
    .cmd-box {{
      background: #1a202c; color: #a0aec0; padding: 1rem 1.2rem;
      border-radius: 8px; font-family: "SFMono-Regular", Menlo, monospace;
      font-size: .82rem; white-space: pre-wrap; word-break: break-all;
      margin-top: 1rem; position: relative;
    }}
    .cmd-copy-btn {{
      position: absolute; top: .5rem; right: .6rem;
      background: #2d3748; color: #a0aec0; border: 1px solid #4a5568;
      border-radius: 4px; padding: .2rem .55rem; font-size: .75rem;
      cursor: pointer;
    }}
    .cmd-copy-btn:hover {{ background: #4a5568; }}
    .filter-row {{
      display: flex; gap: .5rem; align-items: center; margin-bottom: .6rem;
      flex-wrap: wrap;
    }}
    .filter-row input {{ padding: .35rem .6rem; border: 1px solid #cbd5e0;
      border-radius: 5px; font-size: .83rem; width: 220px; }}
    .gap-fill-row {{ margin-bottom: 1rem; }}
    .gap-fill-row label {{ font-size: .85rem; font-weight: 600; margin-right: .5rem; }}
    .gap-fill-row input[type=text] {{
      padding: .35rem .6rem; border: 1px solid #cbd5e0; border-radius: 5px;
      font-size: .85rem; width: 140px;
    }}
    .preset-section {{ margin-bottom: 1.1rem; }}
    .preset-section h3 {{ font-size: .88rem; font-weight: 700; color: #2d3748;
      margin-bottom: .5rem; }}
    .preset-btn {{
      display: inline-block; padding: .35rem .85rem;
      margin: .2rem .3rem .2rem 0; border-radius: 6px;
      border: 1px solid #cbd5e0; background: #f7fafc;
      font-size: .82rem; font-weight: 600; cursor: pointer;
      transition: background .12s, border-color .12s;
    }}
    .preset-btn:hover {{ background: #ebf4ff; border-color: #90cdf4; }}
    .info-chip {{
      display: inline-block; background: #bee3f8; color: #2a69ac;
      font-size: .73rem; font-weight: 600; padding: .12em .5em;
      border-radius: 4px; margin-left: .35rem;
    }}
  </style>
</head>
<body>
<header>
  <h1>GSMM Report — {model_id}</h1>
  <p>Generated: {generated_at}  |  Reactions: {n_rxn:,}  |  Metabolites: {n_met:,}  |  Genes: {n_gene:,}  |  FBA growth: {fba_str}</p>
</header>

<!-- ── main navigation ── -->
<div class="main-tabs">
  <button class="main-tab active" onclick="switchPanel(this,'panel-overview')">Overview</button>
  <button class="main-tab" onclick="switchPanel(this,'panel-reports')">Interactive Reports</button>
  <button class="main-tab" onclick="switchPanel(this,'panel-environment')">Environment Editor</button>
  {crossfeed_tab_btn}
</div>

<!-- ══════════════════════════════════════════════════════════
     PANEL 1 – Overview
     ══════════════════════════════════════════════════════════ -->
<div id="panel-overview" class="panel active">

  <div class="stat-row">
    <div class="stat-box"><div class="val">{n_rxn:,}</div><div class="lbl">Reactions</div></div>
    <div class="stat-box"><div class="val">{n_met:,}</div><div class="lbl">Metabolites</div></div>
    <div class="stat-box"><div class="val">{n_gene:,}</div><div class="lbl">Genes</div></div>
    <div class="stat-box"><div class="val">{n_blk:,}</div><div class="lbl">Blocked</div></div>
    <div class="stat-box"><div class="val">{fba_str}</div><div class="lbl">FBA growth</div></div>
  </div>

  <h2>Static Figures</h2>
  <div class="png-gallery">
{png_gallery}  </div>

  {'<h2>Mutation Analysis Summary</h2>' + mut_stat_cards if mut_stat_cards else ''}
  {'<p style="margin:.5rem 0 .8rem;font-size:.83rem;color:#718096;">Top 20 genes by growth impact (all types). Open the <em>Interactive Reports → Mutation Analysis</em> tab for the full interactive chart.</p>' if mut_table_rows else ''}
  {f'''<div style="overflow-x:auto;margin-top:.5rem;">
  <table class="data-table">
    <thead>
      <tr>
        <th>Gene ID</th><th>Name</th><th>Type</th>
        <th>Reactions</th><th>Growth ratio</th><th>Growth Δ%</th><th>Status</th>
      </tr>
    </thead>
    <tbody>
{mut_table_rows}    </tbody>
  </table>
  </div>''' if mut_table_rows else ''}
</div>

<!-- ══════════════════════════════════════════════════════════
     PANEL 2 – Interactive Reports
     ══════════════════════════════════════════════════════════ -->
<div id="panel-reports" class="panel">
  <div style="margin-bottom:.5rem;">
{iframe_tabs}  </div>
{iframe_panels}
</div>

<!-- ══════════════════════════════════════════════════════════
     PANEL 3 – Environment Editor  (Feature [3])
     ══════════════════════════════════════════════════════════ -->
<div id="panel-environment" class="panel">
  <p style="margin-bottom:1rem;font-size:.9rem;color:#4a5568;">
    Edit exchange reaction bounds to define the growth medium, then copy the
    generated shell command to re-run the pipeline with your changes.
    <span class="info-chip">{len(open_exchanges)} open exchange reactions</span>
  </p>

  <!-- Preset environment buttons -->
  <div class="preset-section">
    <h3>Quick presets</h3>
    {preset_btns_html}
  </div>

  <!-- Gap-fill medium shortcut -->
  <div class="gap-fill-row">
    <label for="gapfill-input">Gap-fill medium</label>
    <input type="text" id="gapfill-input" placeholder="e.g. M9, LB, minimal"
           oninput="updateCmd()" />
    <span style="font-size:.8rem;color:#718096;margin-left:.4rem;">
      (CarveMe medium name — leave blank to skip gap-filling)
    </span>
  </div>

  <!-- Search / filter -->
  <div class="filter-row">
    <input type="text" id="env-search" placeholder="Filter reactions…"
           oninput="filterEnvTable()" />
    <button class="btn btn-secondary" onclick="resetAllBounds()">Reset all</button>
    <button class="btn btn-danger"    onclick="clearAllBounds()">Clear all (set lb=0)</button>
  </div>

  <div style="overflow-x:auto;max-height:55vh;overflow-y:auto;border-radius:8px;">
  <table class="data-table" id="env-table">
    <thead>
      <tr>
        <th>Reaction ID</th>
        <th>Name / metabolite</th>
        <th>Lower bound<br><span style="font-weight:400;opacity:.8">(negative = uptake)</span></th>
        <th>Upper bound</th>
      </tr>
    </thead>
    <tbody id="env-tbody">
{env_rows_html}    </tbody>
  </table>
  </div>

  <!-- Generated command -->
  <h2 style="margin-top:1.25rem;">Generated pipeline command</h2>
  <p style="font-size:.83rem;color:#718096;margin:.4rem 0 .3rem;">
    Copy and run this in your terminal (with the <code>gsmm_pipeline</code> conda env active):
  </p>
  <div class="cmd-box" id="cmd-box">
    <button class="cmd-copy-btn" onclick="copyCmd()">Copy</button>
    <span id="cmd-text"></span>
  </div>
</div>

<script>
let __lastReactionHighlight = null;

// ── panel switching ──────────────────────────────────────────────────────────
function switchPanel(btn, id) {{
  document.querySelectorAll('.main-tab').forEach(b => b.classList.remove('active'));
  document.querySelectorAll('.panel').forEach(p => p.classList.remove('active'));
  btn.classList.add('active');
  document.getElementById(id).classList.add('active');
}}

// ── iframe tab switching ─────────────────────────────────────────────────────
function switchIframe(btn, id) {{
  document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
  document.querySelectorAll('iframe').forEach(f => f.style.display = 'none');
  btn.classList.add('active');
  const frame = document.getElementById(id);
  frame.style.display = 'block';

  if (__lastReactionHighlight && (frame.getAttribute('src') || '').includes('reaction_network.html')) {{
    try {{
      frame.contentWindow.postMessage(__lastReactionHighlight, '*');
      setTimeout(() => frame.contentWindow.postMessage(__lastReactionHighlight, '*'), 250);
    }} catch (e) {{}}
  }}
}}

function _broadcastToReportIframes(message) {{
  document.querySelectorAll('#panel-reports iframe').forEach(frame => {{
    if (!frame || !frame.contentWindow) return;
    try {{
      frame.contentWindow.postMessage(message, '*');
    }} catch (e) {{}}
  }});
}}

window.addEventListener('message', function(event) {{
  const msg = event.data || {{}};
  if (msg.type === 'gsmm-highlight-reactions') {{
    __lastReactionHighlight = msg;
    const reactionBtn = Array.from(document.querySelectorAll('.tab-btn'))
      .find(b => (b.dataset.src || '').includes('reaction_network.html'));
    if (reactionBtn) {{
      switchIframe(reactionBtn, reactionBtn.dataset.target);
    }}
    _broadcastToReportIframes(msg);
  }}
  if (msg.type === 'gsmm-clear-reaction-highlight') {{
    __lastReactionHighlight = null;
    _broadcastToReportIframes(msg);
  }}
}});

// ── environment editor ────────────────────────────────────────────────────────
const ORIG_BOUNDS = {orig_bounds_json};
const ENV_PRESETS = {env_presets_json};

function applyPreset(key) {{
  const preset = ENV_PRESETS[key];
  if (!preset || !preset.bounds) return;
  const bounds = preset.bounds;

  // First reset all to originals, then apply preset overrides
  document.querySelectorAll('#env-tbody tr').forEach(tr => {{
    const rxnId = tr.dataset.rxn;
    const orig  = ORIG_BOUNDS[rxnId];
    if (!orig) return;
    tr.querySelector('.lb-input').value = orig.lb;
    tr.querySelector('.ub-input').value = orig.ub;
    tr.querySelector('.lb-input').classList.remove('changed');
    tr.querySelector('.ub-input').classList.remove('changed');
  }});

  // Apply preset-specific bounds
  Object.entries(bounds).forEach(([rxnId, b]) => {{
    const tr = document.querySelector(`#env-tbody tr[data-rxn="${{rxnId}}"]`);
    if (!tr) return;
    const lbInp = tr.querySelector('.lb-input');
    const ubInp = tr.querySelector('.ub-input');
    const orig  = ORIG_BOUNDS[rxnId];
    if (lbInp && b.lb !== undefined) {{
      lbInp.value = b.lb;
      if (!orig || b.lb !== orig.lb) lbInp.classList.add('changed');
    }}
    if (ubInp && b.ub !== undefined) {{
      ubInp.value = b.ub;
      if (!orig || b.ub !== orig.ub) ubInp.classList.add('changed');
    }}
  }});

  updateCmd();
}}

function filterEnvTable() {{
  const q = document.getElementById('env-search').value.toLowerCase();
  document.querySelectorAll('#env-tbody tr').forEach(tr => {{
    const text = tr.textContent.toLowerCase();
    tr.style.display = text.includes(q) ? '' : 'none';
  }});
}}

function markChanged() {{
  document.querySelectorAll('.lb-input, .ub-input').forEach(inp => {{
    const orig = parseFloat(inp.dataset.orig);
    const cur  = parseFloat(inp.value);
    inp.classList.toggle('changed', !isNaN(cur) && cur !== orig);
  }});
  updateCmd();
}}

document.querySelectorAll('.lb-input, .ub-input').forEach(inp => {{
  inp.addEventListener('input', markChanged);
}});

function resetAllBounds() {{
  document.querySelectorAll('#env-tbody tr').forEach(tr => {{
    const rxnId = tr.dataset.rxn;
    const orig  = ORIG_BOUNDS[rxnId];
    if (!orig) return;
    tr.querySelector('.lb-input').value = orig.lb;
    tr.querySelector('.ub-input').value = orig.ub;
  }});
  document.querySelectorAll('.lb-input, .ub-input').forEach(inp => {{
    inp.classList.remove('changed');
  }});
  updateCmd();
}}

function clearAllBounds() {{
  document.querySelectorAll('.lb-input').forEach(inp => {{
    inp.value = 0;
    inp.classList.add('changed');
  }});
  updateCmd();
}}

function getChangedBounds() {{
  const changes = [];
  document.querySelectorAll('#env-tbody tr').forEach(tr => {{
    const rxnId = tr.dataset.rxn;
    const orig  = ORIG_BOUNDS[rxnId];
    if (!orig) return;
    const newLb = parseFloat(tr.querySelector('.lb-input').value);
    const newUb = parseFloat(tr.querySelector('.ub-input').value);
    if (newLb !== orig.lb || newUb !== orig.ub) {{
      changes.push({{rxnId, newLb, newUb}});
    }}
  }});
  return changes;
}}

function updateCmd() {{
  const gapfill = document.getElementById('gapfill-input').value.trim();
  const changes = getChangedBounds();

  let cmd = 'python generate_model.py \\\\\\n';
  cmd += '    --fasta <path/to/genome.fna> \\\\\\n';
  cmd += '    --gbk   <path/to/genome.gbk> \\\\\\n';
  cmd += '    --output output/model.xml \\\\\\n';
  if (gapfill) {{
    cmd += `    --gap-fill ${{gapfill}} \\\\\\n`;
  }}
  cmd += '    --visualize\\n';

  if (changes.length > 0) {{
    cmd += '\\n# Exchange bound overrides (apply after model is built):\\n';
    cmd += '# python -c "\\n';
    cmd += '#   import cobra\\n';
    cmd += '#   m = cobra.io.read_sbml_model(\\'output/model.xml\\')\\n';
    changes.forEach(c => {{
      cmd += `#   m.reactions.get_by_id(\\'${{c.rxnId}}\\').bounds = (${{c.newLb}}, ${{c.newUb}})\\n`;
    }});
    cmd += '#   cobra.io.write_sbml_model(m, \\'output/model_env.xml\\')\\n';
    cmd += '# "\\n';
  }}

  cmd = cmd.replace('python generate_model.py', 'python src/generate_model.py');
  document.getElementById('cmd-text').textContent = cmd;
}}

function copyCmd() {{
  const text = document.getElementById('cmd-text').textContent;
  navigator.clipboard.writeText(text).then(() => {{
    const btn = document.querySelector('.cmd-copy-btn');
    btn.textContent = 'Copied!';
    setTimeout(() => btn.textContent = 'Copy', 1500);
  }});
}}

// initialise command on load
updateCmd();
</script>
{crossfeed_panel_html}
</body>
</html>
"""

    out_path = out_dir / "index.html"
    out_path.write_text(html, encoding="utf-8")
    log.info("  Saved: %s", out_path)


# ---------------------------------------------------------------------------
# CLI (standalone usage)
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Visualize a COBRApy SBML model comprehensively.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--model",   required=True, help="Path to SBML model (.xml)")
    p.add_argument("--out-dir", default="output/viz", help="Output directory for plots")
    return p.parse_args()


def main() -> None:
    args  = _parse_args()
    model_path = Path(args.model)

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s",
                        datefmt="%H:%M:%S")
    log = logging.getLogger("visualize_model")

    if not model_path.is_file():
        log.error("Model file not found: %s", model_path)
        raise SystemExit(1)

    import cobra
    log.info("Loading model from %s …", model_path)
    # Silence cobra's noisy 'Adding exchange reaction' / 'Ignoring reaction' messages
    logging.getLogger("cobra").setLevel(logging.ERROR)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = cobra.io.read_sbml_model(str(model_path))
    logging.getLogger("cobra").setLevel(logging.WARNING)
    log.info("  Loaded: %d reactions, %d metabolites, %d genes",
             len(model.reactions), len(model.metabolites), len(model.genes))

    run_all_visualizations(model, out_dir=Path(args.out_dir), log=log)


if __name__ == "__main__":
    main()
