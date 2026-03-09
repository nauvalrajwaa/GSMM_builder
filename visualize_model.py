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
  7. escher_map.html            – interactive Escher metabolic map (if escher
                                   is installed and a base-map is available)
  8. model_summary.csv          – machine-readable summary statistics table

Usage (standalone)
------------------
    python visualize_model.py --model output/model.xml --out-dir output/viz

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
    wedges, texts, autotexts = ax.pie(
        sizes,
        labels=labels,
        autopct="%1.1f%%",
        startangle=140,
        colors=palette,
        wedgeprops=dict(linewidth=0.7, edgecolor="white"),
        pctdistance=0.82,
    )
    for at in autotexts:
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

def plot_reaction_network(
    model: "cobra.Model",
    out_dir: Path,
    log: logging.Logger,
    max_nodes: int = 300,
) -> None:
    """
    Build a bipartite reaction-metabolite network and export it as an
    interactive Plotly HTML file.  Large models are sampled to ``max_nodes``
    reactions to keep the browser responsive.
    """
    if not _HAS_PLOTLY:
        log.warning("plotly not installed – skipping interactive network graph.")
        return
    if not _HAS_NX:
        log.warning("networkx not installed – skipping interactive network graph.")
        return

    log.info("  Building reaction network graph …")

    # Sample reactions if model is very large
    reactions = list(model.reactions)
    if len(reactions) > max_nodes:
        import random, math
        random.seed(42)
        reactions = random.sample(reactions, max_nodes)
        log.info("  (sampled %d / %d reactions for legibility)", max_nodes, len(model.reactions))

    G = nx.DiGraph()

    # Add nodes
    for rxn in reactions:
        G.add_node(rxn.id, node_type="reaction",
                   label=rxn.name or rxn.id,
                   subsystem=rxn.subsystem or "unknown")
        for met in rxn.reactants:
            G.add_node(met.id, node_type="metabolite",
                       label=met.name or met.id,
                       compartment=met.compartment or "?")
            G.add_edge(met.id, rxn.id)
        for met in rxn.products:
            G.add_node(met.id, node_type="metabolite",
                       label=met.name or met.id,
                       compartment=met.compartment or "?")
            G.add_edge(rxn.id, met.id)

    # Layout
    try:
        pos = nx.kamada_kawai_layout(G)
    except Exception:
        pos = nx.spring_layout(G, seed=42)

    # Separate node types
    rxn_nodes  = [n for n, d in G.nodes(data=True) if d.get("node_type") == "reaction"]
    met_nodes  = [n for n, d in G.nodes(data=True) if d.get("node_type") == "metabolite"]

    def _make_node_trace(nodes: list[str], color: str, symbol: str,
                         size: int, name: str) -> "go.Scatter":
        xs, ys, texts = [], [], []
        for n in nodes:
            x, y = pos[n]
            xs.append(x); ys.append(y)
            data = G.nodes[n]
            texts.append(
                f"<b>{data.get('label', n)}</b><br>ID: {n}<br>"
                + (f"Subsystem: {data.get('subsystem','')}" if data.get("node_type") == "reaction"
                   else f"Compartment: {data.get('compartment','')}")
            )
        return go.Scatter(
            x=xs, y=ys, mode="markers",
            marker=dict(size=size, color=color, symbol=symbol,
                        line=dict(width=0.5, color="white")),
            text=texts, hoverinfo="text", name=name,
        )

    edge_x, edge_y = [], []
    for u, v in G.edges():
        x0, y0 = pos[u]; x1, y1 = pos[v]
        edge_x += [x0, x1, None]
        edge_y += [y0, y1, None]

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y, mode="lines",
        line=dict(width=0.5, color="#aaa"),
        hoverinfo="none", showlegend=False,
    )

    fig = go.Figure(
        data=[
            edge_trace,
            _make_node_trace(met_nodes,  "#55A868", "circle",  6, "Metabolite"),
            _make_node_trace(rxn_nodes,  "#4C72B0", "square",  9, "Reaction"),
        ],
        layout=go.Layout(
            title=dict(text=f"Reaction–Metabolite Network  ({len(G.nodes())} nodes, "
                            f"{len(G.edges())} edges)",
                       font=dict(size=15)),
            showlegend=True,
            hovermode="closest",
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor="white",
            paper_bgcolor="white",
            margin=dict(l=20, r=20, t=60, b=20),
        ),
    )

    out_html = out_dir / "reaction_network.html"
    fig.write_html(str(out_html))
    log.info("  Saved: %s", out_html)


# ---------------------------------------------------------------------------
# 6b. Interactive Plotly dashboard (overview of all key statistics)
# ---------------------------------------------------------------------------

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
# 8. Machine-readable CSV summary
# ---------------------------------------------------------------------------

def export_summary_csv(
    model: "cobra.Model",
    blocked: list[str],
    fba_value: float | None,
    out_dir: Path,
    log: logging.Logger,
) -> None:
    """Write a one-row CSV of key model statistics."""
    if not _HAS_PD:
        log.warning("pandas not installed – skipping CSV summary.")
        return

    from collections import Counter
    subsystems = Counter(r.subsystem for r in model.reactions if r.subsystem)

    data = {
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

    df = pd.DataFrame([data])
    out_csv = out_dir / "model_summary.csv"
    df.to_csv(out_csv, index=False)
    log.info("  Saved: %s", out_csv)


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

    # Escher map
    export_escher_map(model, out_dir, log)

    # CSV summary
    export_summary_csv(model, blocked, fba_value, out_dir, log)

    log.info("All visualizations complete.  Output directory: %s", out_dir)


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
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = cobra.io.read_sbml_model(str(model_path))
    log.info("  Loaded: %d reactions, %d metabolites, %d genes",
             len(model.reactions), len(model.metabolites), len(model.genes))

    run_all_visualizations(model, out_dir=Path(args.out_dir), log=log)


if __name__ == "__main__":
    main()
