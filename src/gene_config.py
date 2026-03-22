#!/usr/bin/env python3
"""
gene_config.py
==============
Flexible gene-manipulation system for a COBRApy genome-scale metabolic model.

Reads a JSON configuration file that lists genes to **knock out** or **knock in**,
applies the mutations to the model, runs FBA, and writes:
  - ``gene_config_result.json``  – per-gene FBA results
  - ``gene_config_report.html``  – interactive Plotly summary

Config file format (``gene_config.json``)
------------------------------------------
::

    {
      "description": "My experiment",
      "knockouts": [
        {"gene_id": "b0001", "note": "first essential gene"},
        {"gene_id": "b0002"}
      ],
      "knockins": [
        {
          "gene_id": "b3236",
          "reactions": [
            {"reaction_id": "EX_lac__L_e", "lower_bound": -10, "upper_bound": 1000}
          ],
          "note": "force lactate uptake"
        }
      ]
    }

- ``knockouts`` : genes to delete (COBRApy ``gene.knock_out()``).
- ``knockins``  : genes to activate by setting bounds on their associated
  reactions.  Optionally override specific reaction bounds via the
  ``reactions`` list; if omitted, all reactions for that gene have their
  lower bound set to −1000 and upper bound to 1000.

Usage (CLI)
-----------
    python gene_config.py \\
        --model  output/model.xml \\
        --config gene_config.json \\
        --out-dir output/viz

Usage (import)
--------------
    from gene_config import load_gene_config, apply_gene_config
    cfg  = load_gene_config("gene_config.json")
    results = apply_gene_config(model, cfg, log=log)
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import warnings
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import cobra

# ---------------------------------------------------------------------------
# Default config template (written when --init is passed)
# ---------------------------------------------------------------------------

DEFAULT_CONFIG: dict = {
    "description": "Example gene manipulation experiment",
    "knockouts": [
        {"gene_id": "b0001", "note": "replace with a real gene ID from your model"},
    ],
    "knockins": [
        {
            "gene_id": "b0002",
            "note": "replace with a real gene ID from your model",
            "reactions": [
                {
                    "reaction_id": "EX_glc__D_e",
                    "lower_bound": -10,
                    "upper_bound": 1000,
                    "note": "open glucose uptake as an example knockin"
                }
            ]
        }
    ],
}

# ---------------------------------------------------------------------------
# Load / validate config
# ---------------------------------------------------------------------------

def load_gene_config(path: str | Path) -> dict:
    """Load and validate a gene config JSON file.  Returns the parsed dict."""
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"Gene config not found: {p}")
    with open(p) as fh:
        cfg = json.load(fh)
    if not isinstance(cfg, dict):
        raise ValueError("Gene config must be a JSON object.")
    cfg.setdefault("description", "")
    cfg.setdefault("knockouts", [])
    cfg.setdefault("knockins", [])
    return cfg


def write_default_config(path: str | Path) -> None:
    """Write the default config template to *path*."""
    p = Path(path)
    with open(p, "w") as fh:
        json.dump(DEFAULT_CONFIG, fh, indent=2)
    print(f"Default gene config written to: {p}")


# ---------------------------------------------------------------------------
# Apply config to model
# ---------------------------------------------------------------------------

def apply_gene_config(
    model: "cobra.Model",
    cfg: dict,
    log: logging.Logger | None = None,
) -> list[dict]:
    """
    Apply knockout and knockin entries from *cfg* to a **copy** of *model*
    (using COBRApy context managers so the original is never modified) and
    return a list of result dicts, one per gene entry.

    Each result dict contains:
      gene_id, gene_name, type (knockout/knockin), note,
      wt_growth, mut_growth, growth_ratio, growth_change_pct,
      status (optimal/infeasible/error), reactions_affected
    """
    if log is None:
        logging.basicConfig(level=logging.INFO,
                            format="%(asctime)s [%(levelname)s] %(message)s",
                            datefmt="%H:%M:%S")
        log = logging.getLogger("gene_config")

    # Wild-type baseline
    try:
        wt_sol = model.optimize()
        wt_growth = wt_sol.objective_value if wt_sol.status == "optimal" else 0.0
    except Exception as exc:
        log.error("WT FBA failed: %s", exc)
        wt_growth = 0.0

    log.info("  WT growth: %.6f", wt_growth)

    results: list[dict] = []

    # ── knockouts ────────────────────────────────────────────────────────────
    for entry in cfg.get("knockouts", []):
        gid  = entry.get("gene_id", "").strip()
        note = entry.get("note", "")
        if not gid:
            log.warning("  Skipping knockout entry with no gene_id.")
            continue
        gene = model.genes.query(lambda g: g.id == gid or g.name == gid)
        if not gene:
            log.warning("  Knockout: gene '%s' not found in model – skipping.", gid)
            results.append(_err_row(gid, "knockout", note, wt_growth, "gene_not_found"))
            continue
        gene = gene[0]
        rxn_ids = [r.id for r in gene.reactions]
        try:
            with model:
                gene.knock_out()
                sol = model.optimize()
                mut_growth = sol.objective_value if sol.status == "optimal" else 0.0
                status     = sol.status
        except Exception as exc:
            log.warning("  KO %s: FBA error: %s", gid, exc)
            mut_growth, status = 0.0, "error"
        ratio = (mut_growth / wt_growth) if wt_growth > 0 else 0.0
        log.info("  KO  %-12s  growth: %.4f → %.4f  (ratio %.4f)",
                 gid, wt_growth, mut_growth, ratio)
        results.append({
            "gene_id":          gene.id,
            "gene_name":        gene.name or gene.id,
            "type":             "knockout",
            "note":             note,
            "wt_growth":        round(wt_growth, 6),
            "mut_growth":       round(mut_growth, 6),
            "growth_ratio":     round(ratio, 4),
            "growth_change_pct": round((ratio - 1) * 100, 2),
            "essential":        ratio <= 0.05,
            "status":           status,
            "reactions_affected": rxn_ids,
        })

    # ── knockins ──────────────────────────────────────────────────────────────
    for entry in cfg.get("knockins", []):
        gid  = entry.get("gene_id", "").strip()
        note = entry.get("note", "")
        rxn_overrides = entry.get("reactions", [])
        if not gid:
            log.warning("  Skipping knockin entry with no gene_id.")
            continue
        gene = model.genes.query(lambda g: g.id == gid or g.name == gid)
        if not gene:
            log.warning("  Knockin: gene '%s' not found in model – skipping.", gid)
            results.append(_err_row(gid, "knockin", note, wt_growth, "gene_not_found"))
            continue
        gene = gene[0]
        rxn_ids = [r.id for r in gene.reactions]
        try:
            with model:
                if rxn_overrides:
                    # Apply explicit reaction bound overrides
                    for ov in rxn_overrides:
                        rid = ov.get("reaction_id", "")
                        rxn_obj = model.reactions.query(lambda r: r.id == rid)
                        if not rxn_obj:
                            log.warning("    Knockin %s: reaction '%s' not found.", gid, rid)
                            continue
                        rxn_obj = rxn_obj[0]
                        if "lower_bound" in ov:
                            rxn_obj.lower_bound = float(ov["lower_bound"])
                        if "upper_bound" in ov:
                            rxn_obj.upper_bound = float(ov["upper_bound"])
                else:
                    # Default: open all reactions associated with this gene
                    for rxn in gene.reactions:
                        if rxn.lower_bound >= 0:
                            rxn.lower_bound = 0
                            rxn.upper_bound = max(rxn.upper_bound, 1000)
                        else:
                            rxn.lower_bound = -1000
                            rxn.upper_bound = max(rxn.upper_bound, 1000)
                sol = model.optimize()
                mut_growth = sol.objective_value if sol.status == "optimal" else 0.0
                status     = sol.status
        except Exception as exc:
            log.warning("  KI %s: FBA error: %s", gid, exc)
            mut_growth, status = 0.0, "error"
        ratio = (mut_growth / wt_growth) if wt_growth > 0 else 0.0
        log.info("  KI  %-12s  growth: %.4f → %.4f  (ratio %.4f)",
                 gid, wt_growth, mut_growth, ratio)
        results.append({
            "gene_id":          gene.id,
            "gene_name":        gene.name or gene.id,
            "type":             "knockin",
            "note":             note,
            "wt_growth":        round(wt_growth, 6),
            "mut_growth":       round(mut_growth, 6),
            "growth_ratio":     round(ratio, 4),
            "growth_change_pct": round((ratio - 1) * 100, 2),
            "essential":        False,
            "status":           status,
            "reactions_affected": rxn_ids,
        })

    return results


def _err_row(gid: str, typ: str, note: str, wt_growth: float, status: str) -> dict:
    return {
        "gene_id": gid, "gene_name": gid, "type": typ, "note": note,
        "wt_growth": round(wt_growth, 6), "mut_growth": 0.0,
        "growth_ratio": 0.0, "growth_change_pct": -100.0,
        "essential": typ == "knockout", "status": status,
        "reactions_affected": [],
    }


# ---------------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------------

def save_results_json(results: list[dict], out_path: Path) -> None:
    """Write the results list to a JSON file."""
    with open(out_path, "w") as fh:
        json.dump(results, fh, indent=2)


def save_results_html(
    results: list[dict],
    cfg: dict,
    out_path: Path,
    log: logging.Logger,
) -> None:
    """Build an interactive Plotly HTML report from *results*."""
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        import pandas as pd
    except ImportError:
        log.warning("plotly/pandas not available – skipping HTML report.")
        return

    df = pd.DataFrame(results)
    if df.empty:
        log.warning("No results to plot.")
        return

    ko_df = df[df["type"] == "knockout"]
    ki_df = df[df["type"] == "knockin"]
    desc  = cfg.get("description", "")

    has_ko = not ko_df.empty
    has_ki = not ki_df.empty
    n_rows = int(has_ko) + int(has_ki) + 1   # +1 for summary bar

    specs = [[{"type": "bar"}]] * n_rows
    titles = ["Gene mutation summary (growth ratio)"]
    if has_ko:
        titles.append(f"Knockout results  ({len(ko_df)} genes)")
    if has_ki:
        titles.append(f"Knockin results  ({len(ki_df)} genes)")

    fig = make_subplots(
        rows=n_rows, cols=1,
        subplot_titles=titles,
        specs=specs,
        vertical_spacing=0.14,
    )

    # Row 1 – combined bar: all genes, coloured by type + essentiality
    all_genes  = df["gene_name"].tolist()
    all_ratios = df["growth_ratio"].tolist()
    all_colors = []
    for _, row in df.iterrows():
        if row["type"] == "knockout":
            all_colors.append("#C44E52" if row["essential"] else "#55A868")
        else:
            all_colors.append("#4C72B0" if row["growth_change_pct"] > 0 else "#CCB974")

    fig.add_trace(go.Bar(
        x=all_genes,
        y=all_ratios,
        marker_color=all_colors,
        text=[f"{v:.3f}" for v in all_ratios],
        textposition="outside",
        hovertext=[
            f"<b>{r['gene_name']}</b> ({r['gene_id']})<br>"
            f"Type: {r['type']}<br>"
            f"Growth ratio: {r['growth_ratio']:.4f}<br>"
            f"Change: {r['growth_change_pct']:+.2f}%<br>"
            f"Status: {r['status']}<br>"
            f"Note: {r['note'] or '—'}"
            for _, r in df.iterrows()
        ],
        hoverinfo="text",
        showlegend=False,
    ), row=1, col=1)
    fig.add_hline(y=1.0, line_dash="dot", line_color="#888",
                  annotation_text="WT growth", row=1, col=1)
    fig.add_hline(y=0.05, line_dash="dash", line_color="#C44E52",
                  annotation_text="Essential threshold", row=1, col=1)

    cur_row = 2

    # Knockout detail
    if has_ko:
        ko_sorted = ko_df.sort_values("growth_ratio")
        colors_ko = ["#C44E52" if e else "#55A868" for e in ko_sorted["essential"]]
        fig.add_trace(go.Bar(
            x=ko_sorted["gene_name"].tolist(),
            y=ko_sorted["growth_ratio"].tolist(),
            marker_color=colors_ko,
            text=[f"{'ESS' if e else 'non'}" for e in ko_sorted["essential"]],
            textposition="outside",
            hovertext=[
                f"<b>{r['gene_name']}</b><br>Ratio: {r['growth_ratio']:.4f}<br>"
                f"Rxns: {', '.join(r['reactions_affected'][:5])}"
                for _, r in ko_sorted.iterrows()
            ],
            hoverinfo="text",
            showlegend=False,
        ), row=cur_row, col=1)
        fig.add_hline(y=0.05, line_dash="dash", line_color="#C44E52",
                      row=cur_row, col=1)
        cur_row += 1

    # Knockin detail
    if has_ki:
        ki_sorted = ki_df.sort_values("growth_change_pct", ascending=False)
        colors_ki = ["#4C72B0" if v > 0 else "#CCB974"
                     for v in ki_sorted["growth_change_pct"]]
        fig.add_trace(go.Bar(
            x=ki_sorted["gene_name"].tolist(),
            y=ki_sorted["growth_change_pct"].tolist(),
            marker_color=colors_ki,
            text=[f"{v:+.1f}%" for v in ki_sorted["growth_change_pct"]],
            textposition="outside",
            hovertext=[
                f"<b>{r['gene_name']}</b><br>Change: {r['growth_change_pct']:+.2f}%<br>"
                f"Rxns: {', '.join(r['reactions_affected'][:5])}"
                for _, r in ki_sorted.iterrows()
            ],
            hoverinfo="text",
            showlegend=False,
        ), row=cur_row, col=1)
        fig.add_hline(y=0, line_dash="dot", line_color="#888",
                      row=cur_row, col=1)

    title_str = f"Gene Config Report"
    if desc:
        title_str += f" — {desc}"
    n_ko_ess = int(ko_df["essential"].sum()) if has_ko else 0
    title_str += f"<br><sup>{len(df)} genes  |  KO essential: {n_ko_ess}</sup>"

    fig.update_layout(
        title_text=title_str,
        title_font_size=15,
        height=280 * n_rows + 160,
        paper_bgcolor="white",
        plot_bgcolor="white",
    )
    fig.update_yaxes(title_text="Growth ratio (mut/WT)", row=1, col=1)
    if has_ko:
        fig.update_yaxes(title_text="Growth ratio", row=2, col=1)
    if has_ki:
        fig.update_yaxes(title_text="Growth change (%)", row=cur_row - (1 if has_ko else 0), col=1)

    fig.write_html(str(out_path))
    log.info("  Saved: %s", out_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Apply gene knockout/knockin config to a GSMM and report results.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--model",   default="output/model.xml",
                   help="Path to SBML model (.xml)")
    p.add_argument("--config",  default="gene_config.json",
                   help="Path to gene manipulation config JSON")
    p.add_argument("--out-dir", default="output/viz",
                   help="Directory to write output files")
    p.add_argument("--init",    action="store_true",
                   help="Write a default gene_config.json template and exit")
    p.add_argument("--init-path", default="gene_config.json",
                   help="Path for --init template output")
    return p.parse_args()


def main() -> None:
    args = _parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s",
                        datefmt="%H:%M:%S")
    log = logging.getLogger("gene_config")

    if args.init:
        write_default_config(args.init_path)
        sys.exit(0)

    model_path = Path(args.model)
    config_path = Path(args.config)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not model_path.is_file():
        log.error("Model not found: %s", model_path)
        sys.exit(1)

    # Load model
    import cobra
    logging.getLogger("cobra").setLevel(logging.ERROR)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = cobra.io.read_sbml_model(str(model_path))
    logging.getLogger("cobra").setLevel(logging.WARNING)
    log.info("Loaded model: %d reactions, %d genes", len(model.reactions), len(model.genes))

    # Load config
    try:
        cfg = load_gene_config(config_path)
    except FileNotFoundError:
        log.error("Config file not found: %s  (run with --init to create a template)", config_path)
        sys.exit(1)

    log.info("Config: '%s'  |  %d knockouts, %d knockins",
             cfg.get("description", ""), len(cfg["knockouts"]), len(cfg["knockins"]))

    # Apply config
    results = apply_gene_config(model, cfg, log=log)

    # Save outputs
    out_json = out_dir / "gene_config_result.json"
    save_results_json(results, out_json)
    log.info("Saved: %s", out_json)

    out_html = out_dir / "gene_config_report.html"
    save_results_html(results, cfg, out_html, log)

    # Summary
    n_ko = sum(1 for r in results if r["type"] == "knockout")
    n_ki = sum(1 for r in results if r["type"] == "knockin")
    n_ess = sum(1 for r in results if r["type"] == "knockout" and r["essential"])
    log.info("Done. %d knockouts (%d essential), %d knockins.", n_ko, n_ess, n_ki)


if __name__ == "__main__":
    main()
