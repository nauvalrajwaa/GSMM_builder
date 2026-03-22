#!/usr/bin/env python3
"""
crossfeed_network.py
====================
Metabolic Hand-off & Cross-Feeding Network Builder for metagenomics MAG inputs.

Given a set of Metagenome-Assembled Genomes (MAGs) each annotated with KEGG
Orthology (KO) terms or MetaCyc reaction IDs, this module:

  1. Parses each MAG's annotation table (CSV / TSV).
  2. Maps KO terms to KEGG reactions and their substrate / product metabolites
     via the KEGG REST API (with local cache to avoid redundant HTTP calls).
  3. Identifies "metabolic hand-offs": a metabolite that is a *product* in at
     least one MAG and a *substrate* in at least one different MAG.
  4. Builds a NetworkX bipartite graph:
       - Left nodes  → MAGs (type = "mag")
       - Right nodes → metabolites that are exchanged (type = "metabolite")
       - Edges       → directed: MAG_producer → metabolite → MAG_consumer
  5. Exports an interactive Plotly / Cytoscape.js HTML visualisation.

Input format (annotation table)
--------------------------------
A CSV or TSV with at least two columns:

  mag_id    ko_id
  MAG_001   K00001
  MAG_001   K00002
  MAG_002   K00003
  …

Optional extra columns (``gene_id``, ``gene_name``, ``description``) are
preserved but not required.

The ``mag_id`` column names the bin; ``ko_id`` holds the KO accession
(``K#####`` format).  MetaCyc reaction IDs (``RXN-*``) are also accepted and
passed through a MetaCyc-to-KEGG mapping layer.

CLI usage
---------
    # Build network from an annotation table
    python src/crossfeed_network.py --input mag_annotations.csv --out-dir output/crossfeed

    # Use a pre-populated KEGG cache (avoids API calls in offline/test mode)
    python src/crossfeed_network.py --input mag_annotations.csv --cache kegg_cache.json

    # Generate sample fixture and exit
    python src/crossfeed_network.py --generate-sample sample_mag_annotations.csv

Outputs (written to ``out_dir``)
---------------------------------
  crossfeed_network.html  – interactive bipartite graph (Plotly)
  crossfeed_edges.csv     – tabular edge list: producer_mag, metabolite, consumer_mag
  crossfeed_summary.json  – machine-readable summary stats

Dependencies
------------
Required:  (all in gsmm_pipeline conda env)
    requests, networkx, plotly, pandas

Optional (soft-import, gracefully skipped if absent):
    no additional hard requirements

Author note: KEGG API calls are subject to KEGG's terms of use.  Results are
cached in ``<out_dir>/kegg_cache.json`` to minimise requests.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Soft imports
# ---------------------------------------------------------------------------
try:
    import requests
    import requests as _requests
    _HAS_REQUESTS = True
except ImportError:
    _requests = None
    requests = None
    _HAS_REQUESTS = False

try:
    import networkx as nx
    import networkx as _networkx
    _HAS_NX = True
except ImportError:
    _networkx = None
    nx = None
    _HAS_NX = False

try:
    import plotly.graph_objects as go
    import plotly.graph_objects as _plotly_go
    from plotly.subplots import make_subplots
    _HAS_PLOTLY = True
except ImportError:
    _plotly_go = None
    go = None
    make_subplots = None
    _HAS_PLOTLY = False

try:
    import pandas as pd
    _HAS_PD = True
except ImportError:
    _HAS_PD = False

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

KEGG_API_BASE = "https://rest.kegg.jp"
_API_DELAY    = 0.35   # seconds between API calls (KEGG policy: ≤ 3 req/s)
_KO_ID_PATTERN = re.compile(r"^K\d{5}$", re.IGNORECASE)
_EC_ID_PATTERN = re.compile(r"^\d+\.\d+\.\d+\.(?:\d+|-)$")


def _normalize_annotation_id(raw: str) -> str | None:
    token = str(raw).strip().upper().rstrip(";,")
    if not token:
        return None
    if token.startswith("KO:"):
        token = token[3:]
    elif token.startswith("KEGG:"):
        token = token[5:]

    if _KO_ID_PATTERN.match(token):
        return token

    if token.startswith("EC:"):
        ec = token[3:]
        if _EC_ID_PATTERN.match(ec):
            return f"EC:{ec}"
        return None

    if _EC_ID_PATTERN.match(token):
        return f"EC:{token}"

    return None

# ---------------------------------------------------------------------------
# KEGG cache helpers
# ---------------------------------------------------------------------------

def _load_cache(cache_path: Path) -> dict[str, Any]:
    if cache_path.is_file():
        try:
            with open(cache_path, encoding="utf-8") as fh:
                return json.load(fh)
        except Exception:
            pass
    return {}


def _save_cache(cache: dict[str, Any], cache_path: Path) -> None:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w", encoding="utf-8") as fh:
        json.dump(cache, fh, indent=2)


# ---------------------------------------------------------------------------
# KEGG REST helpers
# ---------------------------------------------------------------------------

def _kegg_get(endpoint: str, log: logging.Logger) -> str | None:
    if not _HAS_REQUESTS or requests is None:
        log.warning("requests not installed – KEGG API calls disabled.")
        return None
    url = f"{KEGG_API_BASE}/{endpoint}"
    try:
        resp = requests.get(url, timeout=20)
        if resp.status_code == 200:
            return resp.text
        log.debug("  KEGG GET %s → HTTP %d", url, resp.status_code)
        return None
    except Exception as exc:
        log.debug("  KEGG GET %s failed: %s", url, exc)
        return None


def _fetch_ko_reactions(ko_id: str, cache: dict, log: logging.Logger) -> list[str]:
    normalized = _normalize_annotation_id(ko_id)
    if normalized is None:
        return []

    if normalized.startswith("EC:"):
        ec_id = normalized[3:]
        cache_key = f"ec_rxn:{ec_id}"
        endpoint = f"link/reaction/ec:{ec_id}"
    else:
        cache_key = f"ko_rxn:{normalized}"
        endpoint = f"link/reaction/{normalized}"

    if cache_key in cache:
        return cache[cache_key]

    text = _kegg_get(endpoint, log)
    time.sleep(_API_DELAY)
    rxn_ids: list[str] = []
    if text:
        for line in text.strip().splitlines():
            parts = line.split("\t")
            if len(parts) == 2:
                # format:  ko:K00001\trn:R00001
                rxn_id = parts[1].replace("rn:", "").strip()
                if rxn_id:
                    rxn_ids.append(rxn_id)
    cache[cache_key] = rxn_ids
    return rxn_ids


def _fetch_reaction_compounds(rxn_id: str, cache: dict, log: logging.Logger
                              ) -> tuple[list[str], list[str]]:
    """Return (substrates, products) compound IDs for a KEGG reaction."""
    cache_key = f"rxn:{rxn_id}"
    if cache_key in cache:
        entry = cache[cache_key]
        return entry.get("substrates", []), entry.get("products", [])

    text = _kegg_get(f"get/{rxn_id}", log)
    time.sleep(_API_DELAY)
    substrates: list[str] = []
    products:   list[str] = []

    if text:
        in_equation = False
        equation_text = ""
        for line in text.splitlines():
            if line.startswith("EQUATION"):
                in_equation = True
                equation_text = line[len("EQUATION"):].strip()
            elif in_equation and line.startswith(" "):
                equation_text += " " + line.strip()
            elif in_equation:
                break

        # Parse compound IDs from equation: "C00001 + C00002 <=> C00003 + C00004"
        if "<=>" in equation_text:
            lhs, rhs = equation_text.split("<=>", 1)
            substrates = _parse_compound_ids(lhs)
            products   = _parse_compound_ids(rhs)

    cache[cache_key] = {"substrates": substrates, "products": products}
    return substrates, products


def _parse_compound_ids(side: str) -> list[str]:
    """Extract KEGG compound IDs (C#####) from one side of a reaction equation."""
    import re
    return re.findall(r"C\d{5}", side)


def _fetch_compound_name(cpd_id: str, cache: dict, log: logging.Logger) -> str:
    """Return a human-readable name for a KEGG compound."""
    cache_key = f"cpd_name:{cpd_id}"
    if cache_key in cache:
        return cache[cache_key]

    text = _kegg_get(f"get/{cpd_id}", log)
    time.sleep(_API_DELAY)
    name = cpd_id  # fallback: use raw ID
    if text:
        for line in text.splitlines():
            if line.startswith("NAME"):
                raw = line[len("NAME"):].strip().rstrip(";")
                name = raw.split(";")[0].strip() or cpd_id
                break
    cache[cache_key] = name
    return name


# ---------------------------------------------------------------------------
# Annotation table parser
# ---------------------------------------------------------------------------

def load_mag_annotations(
    table_path: Path,
    log: logging.Logger,
    mag_col:  str = "mag_id",
    ko_col:   str = "ko_id",
) -> dict[str, list[str]]:
    import re
    path = Path(table_path)
    if not path.is_file():
        raise FileNotFoundError(f"Annotation table not found: {path}")

    delimiter = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","

    mag_to_kos: dict[str, list[str]] = defaultdict(list)
    n_rows = 0

    with open(path, encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError("Annotation table appears empty or has no header row.")

        # Case-insensitive column matching
        fields_lower = {f.lower(): f for f in (reader.fieldnames or [])}

        actual_mag_col = fields_lower.get(mag_col.lower())
        actual_ko_col  = fields_lower.get(ko_col.lower())

        if actual_mag_col is None:
            # Fallback: first column
            actual_mag_col = (reader.fieldnames or ["mag_id"])[0]
            log.warning("  Column '%s' not found; using '%s' as MAG ID column.",
                        mag_col, actual_mag_col)
        if actual_ko_col is None:
            # Fallback: second column
            actual_ko_col = (reader.fieldnames or ["ko_id", "ko_id"])[1]
            log.warning("  Column '%s' not found; using '%s' as KO column.",
                        ko_col, actual_ko_col)

        for row in reader:
            n_rows += 1
            mag_id = str(row.get(actual_mag_col, "")).strip()
            ko_raw = str(row.get(actual_ko_col,  "")).strip()
            if not mag_id or not ko_raw:
                continue
            parts = re.split(r"[;,|]", ko_raw)
            for part in parts:
                ann_id = _normalize_annotation_id(part)
                if ann_id:
                    mag_to_kos[mag_id].append(ann_id)

    log.info("  Loaded %d rows; %d MAGs, %d total KO annotations.",
             n_rows, len(mag_to_kos), sum(len(v) for v in mag_to_kos.values()))
    return dict(mag_to_kos)


# ---------------------------------------------------------------------------
# Core: build substrate/product sets per MAG
# ---------------------------------------------------------------------------

def build_mag_metabolite_sets(
    mag_to_kos:  dict[str, list[str]],
    cache:       dict[str, Any],
    log:         logging.Logger,
    max_api_kos: int = 500,
) -> dict[str, dict[str, set[str]]]:
    """
    For each MAG, derive the sets of compounds it can *produce* and *consume*.

    Parameters
    ----------
    mag_to_kos  : output of ``load_mag_annotations``
    cache       : KEGG cache dict (modified in place)
    log         : logger
    max_api_kos : safety cap on the number of distinct KOs to look up via
                  the KEGG API (avoids runaway requests in large datasets)

    Returns
    -------
    dict  mag_id → {"substrates": set[cpd_id], "products": set[cpd_id]}
    """
    all_kos: set[str] = set()
    for kos in mag_to_kos.values():
        all_kos.update(kos)

    log.info("  Resolving %d unique KO/EC terms via KEGG API (cache + live) …",
             len(all_kos))

    # KO → reaction IDs (cached)
    ko_to_rxns: dict[str, list[str]] = {}
    for i, ko in enumerate(sorted(all_kos)[:max_api_kos]):
        ko_to_rxns[ko] = _fetch_ko_reactions(ko, cache, log)
        if (i + 1) % 50 == 0:
            log.debug("    … %d / %d KOs resolved", i + 1,
                      min(len(all_kos), max_api_kos))

    # Reaction → substrates / products (cached)
    all_rxns: set[str] = set()
    for rxns in ko_to_rxns.values():
        all_rxns.update(rxns)

    log.info("  Resolving %d unique KEGG reactions …", len(all_rxns))
    rxn_to_subs: dict[str, list[str]] = {}
    rxn_to_prds: dict[str, list[str]] = {}
    for rxn in sorted(all_rxns):
        subs, prds = _fetch_reaction_compounds(rxn, cache, log)
        rxn_to_subs[rxn] = subs
        rxn_to_prds[rxn] = prds

    # Per-MAG compound sets
    mag_sets: dict[str, dict[str, set[str]]] = {}
    for mag_id, kos in mag_to_kos.items():
        substrates: set[str] = set()
        products:   set[str] = set()
        for ko in kos:
            for rxn in ko_to_rxns.get(ko, []):
                substrates.update(rxn_to_subs.get(rxn, []))
                products.update(rxn_to_prds.get(rxn, []))
        mag_sets[mag_id] = {"substrates": substrates, "products": products}

    return mag_sets


# ---------------------------------------------------------------------------
# Cross-feeding detection
# ---------------------------------------------------------------------------

def find_cross_feeding_edges(
    mag_sets: dict[str, dict[str, set[str]]],
    cache:    dict[str, Any],
    log:      logging.Logger,
) -> list[dict[str, str]]:
    """
    Identify metabolic hand-offs between MAGs.

    A "cross-feeding edge" exists when MAG A can *produce* compound C and
    MAG B (A ≠ B) can *consume* compound C.

    Returns
    -------
    List of dicts with keys:
        producer_mag, compound_id, compound_name, consumer_mag
    """
    log.info("  Identifying cross-feeding metabolite hand-offs …")

    # Map compound → set of producer MAGs and consumer MAGs
    cpd_producers: dict[str, set[str]] = defaultdict(set)
    cpd_consumers: dict[str, set[str]] = defaultdict(set)

    for mag_id, sets in mag_sets.items():
        for cpd in sets["products"]:
            cpd_producers[cpd].add(mag_id)
        for cpd in sets["substrates"]:
            cpd_consumers[cpd].add(mag_id)

    edges: list[dict[str, str]] = []
    exchanged_cpds: set[str] = set()

    for cpd, producers in cpd_producers.items():
        consumers = cpd_consumers.get(cpd, set())
        # Filter: must have at least one consumer different from any producer
        cross_consumers = consumers - producers
        if not cross_consumers:
            continue
        exchanged_cpds.add(cpd)
        for producer in producers:
            for consumer in cross_consumers:
                edges.append({
                    "producer_mag":   producer,
                    "compound_id":    cpd,
                    "compound_name":  "",   # filled below
                    "consumer_mag":   consumer,
                })

    # Fetch compound names (batch, with cache)
    log.info("  Fetching names for %d exchanged compounds …", len(exchanged_cpds))
    cpd_name_map: dict[str, str] = {}
    for cpd in sorted(exchanged_cpds):
        cpd_name_map[cpd] = _fetch_compound_name(cpd, cache, log)

    for edge in edges:
        edge["compound_name"] = cpd_name_map.get(edge["compound_id"],
                                                   edge["compound_id"])

    log.info("  Found %d cross-feeding edges across %d exchanged metabolites.",
             len(edges), len(exchanged_cpds))
    return edges


# ---------------------------------------------------------------------------
# NetworkX bipartite graph builder
# ---------------------------------------------------------------------------

def build_bipartite_graph(
    edges: list[dict[str, str]],
    mag_sets: dict[str, dict[str, set[str]]],
    log: logging.Logger,
) -> Any | None:
    """
    Build a directed bipartite graph: MAG → metabolite → MAG.

    Node attributes:
      - type      : "mag" | "metabolite"
      - n_kos     : (MAG nodes) number of KO annotations
      - n_products: (MAG nodes) number of produced compounds
      - n_substrates: (MAG nodes) number of consumed compounds
      - compound_id: (metabolite nodes)

    Edge attributes:
      - role : "produces" | "consumed_by"
    """
    if not _HAS_NX or nx is None:
        log.warning("networkx not installed – skipping graph construction.")
        return None

    G = nx.DiGraph()

    # Add MAG nodes
    for mag_id, sets in mag_sets.items():
        G.add_node(
            mag_id,
            type="mag",
            n_products=len(sets["products"]),
            n_substrates=len(sets["substrates"]),
        )

    # Add metabolite nodes and edges
    seen_pairs: set[tuple[str, str, str]] = set()
    for edge in edges:
        producer = edge["producer_mag"]
        consumer = edge["consumer_mag"]
        cpd_id   = edge["compound_id"]
        cpd_name = edge["compound_name"] or cpd_id

        # Metabolite node (keyed by compound id)
        if cpd_id not in G:
            G.add_node(cpd_id, type="metabolite",
                       compound_id=cpd_id, compound_name=cpd_name)

        # producer → metabolite
        key_prod = (producer, cpd_id, "produces")
        if key_prod not in seen_pairs:
            G.add_edge(producer, cpd_id, role="produces")
            seen_pairs.add(key_prod)

        # metabolite → consumer
        key_cons = (cpd_id, consumer, "consumed_by")
        if key_cons not in seen_pairs:
            G.add_edge(cpd_id, consumer, role="consumed_by")
            seen_pairs.add(key_cons)

    log.info("  Graph: %d nodes (%d MAG, %d metabolite), %d edges",
             G.number_of_nodes(),
             sum(1 for n, d in G.nodes(data=True) if d.get("type") == "mag"),
             sum(1 for n, d in G.nodes(data=True) if d.get("type") == "metabolite"),
             G.number_of_edges())
    return G


# ---------------------------------------------------------------------------
# Plotly bipartite visualisation
# ---------------------------------------------------------------------------

def _mag_color_palette(n: int) -> list[str]:
    """Return n distinct colors for MAG nodes."""
    palette = [
        "#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2",
        "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD",
        "#3A87C8", "#E6824B", "#4CAF50", "#E53935", "#7B61B5",
    ]
    if n <= len(palette):
        return palette[:n]
    # Repeat with slight brightness variation if many MAGs
    return [palette[i % len(palette)] for i in range(n)]


def plot_crossfeed_network(
    G:       Any,
    edges:   list[dict[str, str]],
    mag_sets: dict[str, dict[str, set[str]]],
    out_dir: Path,
    log:     logging.Logger,
) -> Path | None:
    if not _HAS_PLOTLY or go is None:
        log.warning("plotly not installed – skipping cross-feeding network plot.")
        return None
    if not _HAS_NX or nx is None:
        log.warning("networkx not installed – skipping cross-feeding network plot.")
        return None
    if G is None or G.number_of_nodes() == 0:
        log.warning("  Empty graph – no cross-feeding network to visualise.")
        return None

    mag_nodes = sorted([n for n, d in G.nodes(data=True) if d.get("type") == "mag"])
    met_nodes = sorted([n for n, d in G.nodes(data=True) if d.get("type") == "metabolite"])
    n_mag = len(mag_nodes)
    n_met = len(met_nodes)
    if n_mag == 0 or n_met == 0:
        log.warning("  Cross-feeding graph lacks MAG or metabolite nodes – skipping plot.")
        return None

    import math

    mag_colors = _mag_color_palette(n_mag)
    mag_color_map = dict(zip(mag_nodes, mag_colors))

    cpd_label: dict[str, str] = {}
    for n, d in G.nodes(data=True):
        if d.get("type") != "metabolite":
            continue
        compound_id = str(d.get("compound_id", n))
        compound_name_raw = d.get("compound_name")
        compound_name = str(compound_name_raw).strip() if compound_name_raw else compound_id
        cpd_label[compound_id] = compound_name

    def _to_float(val: Any) -> float | None:
        try:
            if val is None:
                return None
            return float(val)
        except Exception:
            return None

    compound_support: dict[str, int] = defaultdict(int)
    for e in edges:
        compound_support[str(e.get("compound_id", ""))] += 1

    raw_values: list[float] = []
    edge_metric: dict[tuple[str, str, str], float] = {}
    edge_metric_is_flux: dict[tuple[str, str, str], bool] = {}
    for e in edges:
        producer = str(e.get("producer_mag", ""))
        consumer = str(e.get("consumer_mag", ""))
        cpd_id = str(e.get("compound_id", ""))
        metric = None
        found_flux = False
        for k in ("flux", "flux_value", "value", "weight", "score", "mmol_gdcw_h"):
            metric = _to_float(e.get(k))
            if metric is not None:
                found_flux = True
                break
        if metric is None:
            metric = float(max(1, compound_support.get(cpd_id, 1)))
        metric = abs(metric)
        raw_values.append(metric)
        edge_metric[(producer, cpd_id, consumer)] = metric
        edge_metric_is_flux[(producer, cpd_id, consumer)] = found_flux

    min_metric = min(raw_values) if raw_values else 1.0
    max_metric = max(raw_values) if raw_values else 1.0

    def _edge_width(metric: float) -> float:
        if max_metric <= min_metric + 1e-12:
            return 3.8
        norm = (metric - min_metric) / (max_metric - min_metric)
        return 1.8 + 7.2 * norm

    mag_radius = max(5.6, 4.8 + 0.95 * math.sqrt(n_mag))
    met_radius = max(1.4, mag_radius * 0.26)
    mag_pos: dict[str, tuple[float, float]] = {}
    for i, m in enumerate(mag_nodes):
        a = (2.0 * math.pi * i) / max(n_mag, 1)
        mag_pos[m] = (mag_radius * math.cos(a), mag_radius * math.sin(a))

    met_pos: dict[str, tuple[float, float]] = {}
    for j, c in enumerate(met_nodes):
        producers = [p for p in G.predecessors(c) if p in mag_pos]
        consumers = [q for q in G.successors(c) if q in mag_pos]
        linked = producers + consumers
        if linked:
            mx = sum(mag_pos[m][0] for m in linked) / len(linked)
            my = sum(mag_pos[m][1] for m in linked) / len(linked)
            scale = 0.31 + 0.14 / max(len(linked), 1)
            jitter_ang = (2.0 * math.pi * (j + 1)) / max(n_met, 1)
            jitter = 0.045 * met_radius
            met_pos[c] = (
                max(-met_radius * 1.6, min(met_radius * 1.6, mx * scale + jitter * math.cos(jitter_ang))),
                max(-met_radius * 1.6, min(met_radius * 1.6, my * scale + jitter * math.sin(jitter_ang))),
            )
        else:
            a = (2.0 * math.pi * j) / max(n_met, 1)
            met_pos[c] = (met_radius * math.cos(a), met_radius * math.sin(a))

    export_color = "#f97316"
    import_color = "#0ea5e9"

    edge_traces: list[Any] = []
    edge_hover_x: list[float] = []
    edge_hover_y: list[float] = []
    edge_hover_text: list[str] = []
    annotations: list[dict[str, Any]] = []

    for e in edges:
        producer = str(e.get("producer_mag", ""))
        consumer = str(e.get("consumer_mag", ""))
        cpd_id = str(e.get("compound_id", ""))
        cpd_name = str(e.get("compound_name", "") or cpd_label.get(cpd_id, cpd_id))
        if producer not in mag_pos or consumer not in mag_pos or cpd_id not in met_pos:
            continue

        edge_key = (producer, cpd_id, consumer)
        metric = edge_metric.get(edge_key, 1.0)
        metric_label = "Flux" if edge_metric_is_flux.get(edge_key, False) else "Support proxy"
        w = _edge_width(metric)

        x0, y0 = mag_pos[producer]
        x1, y1 = met_pos[cpd_id]
        edge_traces.append(go.Scatter(
            x=[x0, x1, None], y=[y0, y1, None],
            mode="lines",
            line=dict(width=w, color=export_color),
            hoverinfo="none",
            showlegend=False,
        ))
        annotations.append(dict(
            x=x1, y=y1, ax=x0, ay=y0,
            xref="x", yref="y", axref="x", ayref="y",
            showarrow=True,
            arrowhead=2, arrowsize=1.08,
            arrowwidth=max(1.0, 0.62 * w),
            arrowcolor=export_color,
            opacity=0.92,
        ))
        edge_hover_x.append((x0 + x1) / 2)
        edge_hover_y.append((y0 + y1) / 2)
        edge_hover_text.append(
            f"<b>Export</b><br>MAG: {producer}<br>Metabolite: {cpd_name} ({cpd_id})"
            f"<br>{metric_label}: {metric:.4g}"
        )

        x2, y2 = met_pos[cpd_id]
        x3, y3 = mag_pos[consumer]
        edge_traces.append(go.Scatter(
            x=[x2, x3, None], y=[y2, y3, None],
            mode="lines",
            line=dict(width=max(1.0, 0.76 * w), color=import_color),
            hoverinfo="none",
            showlegend=False,
        ))
        annotations.append(dict(
            x=x3, y=y3, ax=x2, ay=y2,
            xref="x", yref="y", axref="x", ayref="y",
            showarrow=True,
            arrowhead=2, arrowsize=1.08,
            arrowwidth=max(1.0, 0.52 * w),
            arrowcolor=import_color,
            opacity=0.92,
        ))
        edge_hover_x.append((x2 + x3) / 2)
        edge_hover_y.append((y2 + y3) / 2)
        edge_hover_text.append(
            f"<b>Import</b><br>Metabolite: {cpd_name} ({cpd_id})<br>MAG: {consumer}"
            f"<br>{metric_label}: {metric:.4g}"
        )

    edge_hover_trace = go.Scatter(
        x=edge_hover_x,
        y=edge_hover_y,
        mode="markers",
        marker=dict(size=13, color="rgba(0,0,0,0)", line=dict(width=0)),
        hovertext=edge_hover_text,
        hoverinfo="text",
        showlegend=False,
        name="Edge details",
    )

    mag_x = [mag_pos[m][0] for m in mag_nodes]
    mag_y = [mag_pos[m][1] for m in mag_nodes]
    mag_sizes: list[float] = []
    mag_hover: list[str] = []
    for m in mag_nodes:
        products = len(mag_sets.get(m, {}).get("products", set()))
        substrates = len(mag_sets.get(m, {}).get("substrates", set()))
        out_e = sum(1 for _ in G.successors(m))
        in_e = sum(1 for p in G.predecessors(m) if G.nodes[p].get("type") == "metabolite")
        complexity = products + substrates
        mag_sizes.append(min(100.0, max(70.0, 70.0 + 2.2 * math.sqrt(max(complexity, 1)))))
        mag_hover.append(
            f"<b>{m}</b><br>"
            f"Total exports: {products}<br>"
            f"Total imports: {substrates}<br>"
            f"Outgoing edges: {out_e}<br>"
            f"Incoming edges: {in_e}"
        )

    mag_trace = go.Scatter(
        x=mag_x, y=mag_y,
        mode="markers+text",
        text=mag_nodes,
        textposition="middle center",
        textfont=dict(size=12, color="#0f172a"),
        marker=dict(
            size=mag_sizes,
            color=[mag_color_map[m] for m in mag_nodes],
            line=dict(width=3, color="#f8fafc"),
            symbol="circle",
            opacity=0.95,
        ),
        hovertext=mag_hover,
        hoverinfo="text",
        name="MAG (genome)",
    )

    met_x = [met_pos[c][0] for c in met_nodes]
    met_y = [met_pos[c][1] for c in met_nodes]
    met_sizes = [max(13.0, min(16.0, 14.0 + 0.25 * (G.in_degree(c) + G.out_degree(c)))) for c in met_nodes]
    met_hover = [
        f"<b>{cpd_label.get(c, c)}</b><br>ID: {c}<br>"
        f"Produced by MAGs: {', '.join(list(G.predecessors(c))[:10]) or '—'}<br>"
        f"Consumed by MAGs: {', '.join(list(G.successors(c))[:10]) or '—'}"
        for c in met_nodes
    ]
    def _met_label(c: str) -> str:
        label = cpd_label.get(c)
        if not label:
            label = c
        return label if len(label) < 20 else label[:17] + "…"

    met_trace = go.Scatter(
        x=met_x, y=met_y,
        mode="markers+text",
        text=[_met_label(c) for c in met_nodes],
        textposition="middle right",
        textfont=dict(size=10, color="#334155"),
        marker=dict(
            size=met_sizes,
            color="#fde68a",
            line=dict(width=1.8, color="#b45309"),
            symbol="diamond",
        ),
        hovertext=met_hover,
        hoverinfo="text",
        name="Exchanged metabolite",
    )

    flow_legend_traces = [
        go.Scatter(x=[None], y=[None], mode="lines", line=dict(color=export_color, width=3),
                   name="Export (MAG → metabolite)", hoverinfo="none"),
        go.Scatter(x=[None], y=[None], mode="lines", line=dict(color=import_color, width=3),
                   name="Import (metabolite → MAG)", hoverinfo="none"),
    ]

    n_edges_str = f"{len(edges)} hand-off pairs | {n_mag} MAGs | {n_met} metabolites"
    fig = go.Figure(data=edge_traces + [edge_hover_trace, met_trace, mag_trace] + flow_legend_traces)
    lim = mag_radius + 1.2
    fig.update_layout(
        title=dict(
            text=(
                "Community Cross-Feeding Network — Giant Factory Layout<br>"
                f"<sup>{n_edges_str}</sup>"
            ),
            font=dict(size=19),
            x=0.5,
        ),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-lim, lim]),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-lim, lim]),
        annotations=annotations,
        plot_bgcolor="#f8fafc",
        paper_bgcolor="#f8fafc",
        hovermode="closest",
        legend=dict(
            title="Flow semantics",
            borderwidth=1,
            bordercolor="#e2e8f0",
            bgcolor="rgba(255,255,255,0.9)",
            x=1.02, y=1.0,
        ),
        margin=dict(l=30, r=210, t=95, b=30),
        height=max(780, 120 * max(4, n_mag)),
    )

    fig.add_annotation(
        x=0.5,
        y=1.08,
        xref="paper",
        yref="paper",
        text=(
            "<b>Design:</b> giant MAG circles are core players; small diamonds are exchanged metabolites. "
            "Edge color indicates export/import direction; edge width scales with measured flux when available "
            "or structural support otherwise."
        ),
        showarrow=False,
        font=dict(size=11, color="#334155"),
    )

    out_path = out_dir / "crossfeed_network.html"
    fig.write_html(
        str(out_path),
        include_plotlyjs="cdn",
        config={"displayModeBar": True, "scrollZoom": True, "displaylogo": False},
    )
    log.info("  Saved: %s", out_path)
    return out_path


# ---------------------------------------------------------------------------
# Edge list CSV export
# ---------------------------------------------------------------------------

def export_edge_csv(
    edges:   list[dict[str, str]],
    out_dir: Path,
    log:     logging.Logger,
) -> Path:
    """Write crossfeed_edges.csv with columns: producer_mag, compound_id, compound_name, consumer_mag."""
    out_path = out_dir / "crossfeed_edges.csv"
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["producer_mag", "compound_id", "compound_name", "consumer_mag"],
        )
        writer.writeheader()
        writer.writerows(edges)
    log.info("  Saved: %s  (%d edges)", out_path, len(edges))
    return out_path


# ---------------------------------------------------------------------------
# Summary JSON export
# ---------------------------------------------------------------------------

def export_summary_json(
    edges:    list[dict[str, str]],
    mag_sets: dict[str, dict[str, set[str]]],
    G:        Any | None,
    out_dir:  Path,
    log:      logging.Logger,
) -> Path:
    """Write crossfeed_summary.json with high-level statistics."""
    from collections import Counter

    cpd_counts = Counter(e["compound_id"] for e in edges)
    top_metabolites = [
        {"compound_id": cid, "compound_name": edges[[e["compound_id"] for e in edges].index(cid)]["compound_name"],
         "n_edges": cnt}
        for cid, cnt in cpd_counts.most_common(20)
    ]

    summary = {
        "n_mags":        len(mag_sets),
        "n_edges":       len(edges),
        "n_metabolites": len(cpd_counts),
        "mag_ids":       sorted(mag_sets.keys()),
        "top_metabolites": top_metabolites,
        "mag_stats": {
            mag_id: {
                "n_products":   len(sets.get("products",   set())),
                "n_substrates": len(sets.get("substrates", set())),
            }
            for mag_id, sets in mag_sets.items()
        },
    }
    if G is not None and _HAS_NX:
        summary["graph_stats"] = {
            "n_nodes": G.number_of_nodes(),
            "n_graph_edges": G.number_of_edges(),
        }

    out_path = out_dir / "crossfeed_summary.json"
    out_path.write_text(json.dumps(summary, indent=2, default=list),
                        encoding="utf-8")
    log.info("  Saved: %s", out_path)
    return out_path


# ---------------------------------------------------------------------------
# Taxonomy colour helper
# ---------------------------------------------------------------------------

_TAXONOMY_COLORS: dict[str, str] = {
    "bifidobacterium": "#4CAF50",
    "lactobacillus":   "#8BC34A",
    "firmicutes":      "#2196F3",
    "akkermansia":     "#FF9800",
    "bacteroides":     "#9C27B0",
    "ruminococcus":    "#F44336",
    "faecalibacterium":"#00BCD4",
    "methanobrevibacter":"#FF5722",
    "desulfovibrio":   "#795548",
    "prevotella":      "#607D8B",
    "clostridium":     "#3F51B5",
    "eubacterium":     "#009688",
    "roseburia":       "#CDDC39",
    "blautia":         "#FFC107",
}


def _taxonomy_color(mag_id: str, fallback: str = "#607D8B") -> str:
    """Return a colour for a MAG based on the genus token in its name."""
    token = mag_id.split("_")[-1].lower()
    return _TAXONOMY_COLORS.get(token, fallback)


def _taxonomy_colors_for_mags(mag_ids: list) -> dict:
    """Return {mag_id: hex_color} using taxonomy names embedded in IDs."""
    palette = _mag_color_palette(len(mag_ids))
    result = {}
    for i, m in enumerate(mag_ids):
        result[m] = _taxonomy_color(m, fallback=palette[i % len(palette)])
    return result


# ---------------------------------------------------------------------------
# Sankey diagram
# ---------------------------------------------------------------------------



def _hex_to_rgba(hex_color: str, alpha: float = 0.65) -> str:
    """Convert #rrggbb to rgba(r,g,b,alpha) string for Plotly compatibility."""
    h = hex_color.lstrip("#")
    if len(h) == 6:
        r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
        return f"rgba({r},{g},{b},{alpha})"
    return hex_color  # fallback: return as-is

def plot_sankey_diagram(
    edges: list,
    mag_sets: dict,
    out_dir: Path,
    log: logging.Logger,
    abundance: dict | None = None,
) -> "Path | None":
    """Write crossfeed_sankey.html — metabolic flow Sankey diagram."""
    if not _HAS_PLOTLY or go is None:
        log.warning("plotly not available – skipping Sankey diagram")
        return None
    if len(edges) < 2:
        log.warning("Too few edges (%d) for Sankey diagram", len(edges))
        return None

    # Build node list (deduplicated, stable order)
    seen: dict[str, int] = {}
    labels: list[str] = []

    def _idx(name: str) -> int:
        if name not in seen:
            seen[name] = len(labels)
            labels.append(name)
        return seen[name]

    mag_ids = sorted(mag_sets.keys())
    color_map = _taxonomy_colors_for_mags(mag_ids)

    link_source: list[int] = []
    link_target: list[int] = []
    link_value:  list[float] = []
    link_color:  list[str] = []
    link_label:  list[str] = []

    for e in edges:
        prod = e["producer_mag"]
        cpd  = e.get("compound_name") or e["compound_id"]
        cons = e["consumer_mag"]
        weight = float(abundance.get(prod, 1)) if abundance else 1.0

        # producer → metabolite
        link_source.append(_idx(prod))
        link_target.append(_idx(cpd))
        link_value.append(weight)
        link_color.append(_hex_to_rgba(color_map.get(prod, "#888888")))
        link_label.append(f"{prod} → {cpd}")

        # metabolite → consumer
        link_source.append(_idx(cpd))
        link_target.append(_idx(cons))
        link_value.append(weight)
        link_color.append("rgba(187,187,187,0.5)")
        link_label.append(f"{cpd} → {cons}")

    # Node colours: MAGs use taxonomy color, metabolites use soft gold
    node_colors = []
    for lbl in labels:
        if lbl in color_map:
            node_colors.append(color_map[lbl])
        else:
            node_colors.append("#ecc94b")

    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(
            pad=18,
            thickness=22,
            line=dict(color="#333", width=0.8),
            label=labels,
            color=node_colors,
            hovertemplate="<b>%{label}</b><br>Flow: %{value:.1f}<extra></extra>",
        ),
        link=dict(
            source=link_source,
            target=link_target,
            value=link_value,
            color=link_color,
            label=link_label,
            hovertemplate="%{label}<br>Flow: %{value:.1f}<extra></extra>",
        ),
    ))
    fig.update_layout(
        title_text="Metabolic Flow — Community Cross-Feeding",
        title_font_size=16,
        font_size=12,
        margin=dict(l=20, r=20, t=60, b=20),
        paper_bgcolor="#fafafa",
    )

    out_path = out_dir / "crossfeed_sankey.html"
    fig.write_html(str(out_path), include_plotlyjs="cdn")
    log.info("  Saved Sankey diagram: %s", out_path)
    return out_path


# ---------------------------------------------------------------------------
# Pathway completeness heatmap
# ---------------------------------------------------------------------------


def plot_pathway_heatmap(
    mag_to_kos: dict,
    out_dir: Path,
    log: logging.Logger,
    ko_descriptions: dict | None = None,
) -> "Path | None":
    """Write crossfeed_heatmap.html — KO presence mosaic across MAGs."""
    if not _HAS_PLOTLY or go is None:
        log.warning("plotly not available – skipping heatmap")
        return None
    if not mag_to_kos:
        log.warning("No MAG annotations – skipping heatmap")
        return None

    ko_descriptions = ko_descriptions or {}
    mag_ids = sorted(mag_to_kos.keys())

    # Collect all KOs; sort by prevalence (most common first)
    ko_counter: dict[str, int] = {}
    for kos in mag_to_kos.values():
        for k in kos:
            ko_counter[k] = ko_counter.get(k, 0) + 1
    all_kos = sorted(ko_counter, key=lambda k: -ko_counter[k])

    # Build presence matrix: rows=KOs, cols=MAGs
    z: list[list[int]] = []
    hover: list[list[str]] = []
    for ko in all_kos:
        row_z = []
        row_h = []
        for mag in mag_ids:
            present = 1 if ko in mag_to_kos.get(mag, set()) else 0
            row_z.append(present)
            desc = ko_descriptions.get(ko, "")
            status = "Present ✓" if present else "Absent ✗"
            row_h.append(f"MAG: {mag}<br>KO: {ko}<br>Gene: {desc or '—'}<br>{status}")
        z.append(row_z)
        hover.append(row_h)

    # Completeness % per MAG (bottom annotation)
    completeness = [
        round(sum(z[r][c] for r in range(len(all_kos))) / max(1, len(all_kos)) * 100)
        for c in range(len(mag_ids))
    ]

    fig = go.Figure(go.Heatmap(
        z=z,
        x=mag_ids,
        y=all_kos,
        text=hover,
        hoverinfo="text",
        colorscale=[[0, "#f7fafc"], [1, "#2dd4bf"]],
        zmin=0, zmax=1,
        showscale=False,
        xgap=2, ygap=1,
    ))

    # Completeness annotations at the bottom
    for i, (mag, pct) in enumerate(zip(mag_ids, completeness)):
        fig.add_annotation(
            x=mag, y=-0.5,
            text=f"<b>{pct}%</b>",
            showarrow=False,
            yref="paper",
            font=dict(size=10, color="#555"),
        )

    n_kos = len(all_kos)
    n_mags = len(mag_ids)
    fig.update_layout(
        title_text="Pathway Completeness Heatmap — KO Presence by MAG",
        title_font_size=15,
        xaxis=dict(side="top", tickangle=-35, title=""),
        yaxis=dict(autorange="reversed", title="KO (Enzymatic Step)", tickfont=dict(size=9)),
        height=max(500, 20 * n_kos + 120),
        width=max(600, 80 * n_mags + 120),
        margin=dict(l=120, r=40, t=130, b=60),
        paper_bgcolor="#fafafa",
        plot_bgcolor="#fafafa",
    )

    out_path = out_dir / "crossfeed_heatmap.html"
    fig.write_html(str(out_path), include_plotlyjs="cdn")
    log.info("  Saved pathway heatmap: %s", out_path)
    return out_path


# ---------------------------------------------------------------------------
# Keystone species & dependency report
# ---------------------------------------------------------------------------


def compute_keystone_report(
    edges: list,
    mag_sets: dict,
    out_dir: Path,
    log: logging.Logger,
) -> "dict | None":
    """Compute Keystone Species and Dependency scores; write CSV + HTML report."""
    from collections import defaultdict
    import csv as _csv
    from datetime import datetime

    if not edges:
        log.warning("No edges – skipping keystone report")
        return None

    mag_ids = sorted(mag_sets.keys())

    # cpd → set of producers / consumers
    cpd_producers: dict[str, set] = defaultdict(set)
    cpd_consumers: dict[str, set] = defaultdict(set)
    cpd_names: dict[str, str] = {}
    for e in edges:
        cpd_producers[e["compound_id"]].add(e["producer_mag"])
        cpd_consumers[e["compound_id"]].add(e["consumer_mag"])
        cpd_names[e["compound_id"]] = e.get("compound_name") or e["compound_id"]

    # Dependency score: fraction of substrates supplied by OTHER MAGs
    dep_scores: dict[str, float] = {}
    for mag in mag_ids:
        subs = mag_sets[mag].get("substrates", set())
        received = sum(
            1 for cid, prods in cpd_producers.items()
            if mag in cpd_consumers.get(cid, set())
            and (prods - {mag})  # at least one OTHER producer
        )
        dep_scores[mag] = round(received / max(1, len(subs)), 3)

    # Provision score: (metabolite, consumer) pairs provided
    provision: dict[str, int] = defaultdict(int)
    for e in edges:
        provision[e["producer_mag"]] += 1

    # Keystone detection: sole producer of ≥1 metabolite consumed by ≥1 other MAG
    keystone_flags: dict[str, bool] = {m: False for m in mag_ids}
    keystone_metabolites: dict[str, list] = {m: [] for m in mag_ids}
    keystone_alerts: list[str] = []
    for cid, producers in cpd_producers.items():
        consumers = cpd_consumers.get(cid, set())
        cross_consumers = consumers - producers
        if len(producers) == 1 and cross_consumers:
            sole = next(iter(producers))
            keystone_flags[sole] = True
            keystone_metabolites[sole].append(cpd_names.get(cid, cid))
            alert = (
                f"⚠️  WARNING: {sole} is a Keystone species. "
                f"Removal would cut supply of '{cpd_names.get(cid, cid)}' "
                f"to {len(cross_consumers)} MAG(s): {', '.join(sorted(cross_consumers))}"
            )
            keystone_alerts.append(alert)
            log.warning(alert)

    # Aggregate stats
    rows = []
    for mag in mag_ids:
        n_provided = provision.get(mag, 0)
        n_received = sum(
            1 for cid in cpd_consumers
            if mag in cpd_consumers[cid] and (cpd_producers[cid] - {mag})
        )
        rows.append({
            "mag_id": mag,
            "dependency_score": dep_scores[mag],
            "provision_score": round(n_provided / max(1, len(mag_ids) - 1), 3),
            "n_metabolites_provided": n_provided,
            "n_metabolites_received": n_received,
            "is_keystone": keystone_flags[mag],
            "keystone_metabolites": "|".join(keystone_metabolites[mag]),
        })

    # Write CSV
    csv_path = out_dir / "keystone_report.csv"
    fieldnames = [
        "mag_id", "dependency_score", "provision_score",
        "n_metabolites_provided", "n_metabolites_received",
        "is_keystone", "keystone_metabolites",
    ]
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = _csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    log.info("  Saved keystone CSV: %s", csv_path)

    # Write HTML report
    n_keystone = sum(1 for r in rows if r["is_keystone"])
    ts = datetime.now().strftime("%Y-%m-%d %H:%M")

    alert_html = ""
    for alert in keystone_alerts:
        alert_html += f'<div class="alert">{alert}</div>\n'

    table_rows = ""
    for i, r in enumerate(rows):
        ks_badge = ('<span class="badge">Keystone</span>' if r["is_keystone"] else "")
        row_class = "even" if i % 2 == 0 else "odd"
        table_rows += (
            f'<tr class="{row_class}">\n'
            f'  <td>{r["mag_id"]}</td>\n'
            f'  <td>{r["dependency_score"]:.3f}</td>\n'
            f'  <td>{r["provision_score"]:.3f}</td>\n'
            f'  <td>{r["n_metabolites_provided"]}</td>\n'
            f'  <td>{r["n_metabolites_received"]}</td>\n'
            f'  <td>{ks_badge}</td>\n'
            f'  <td class="ks-mets">{r["keystone_metabolites"].replace("|", ", ") or "—"}</td>\n'
            f'</tr>\n'
        )

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Keystone Species Report</title>
<style>
  body {{font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; margin: 0; background:#f0f4f8; color:#2d3748;}}
  header {{background:#1a202c; color:#fff; padding:24px 32px;}}
  header h1 {{margin:0; font-size:1.4rem;}}
  header p {{margin:4px 0 0; font-size:.85rem; color:#a0aec0;}}
  .content {{padding:24px 32px;}}
  .summary {{background:#fff; border-radius:8px; padding:16px 24px; margin-bottom:20px;
             box-shadow:0 1px 4px rgba(0,0,0,.08);}}
  .alert {{background:#fff5f5; border-left:4px solid #fc8181; border-radius:6px;
           padding:12px 18px; margin-bottom:10px; font-size:.88rem; line-height:1.5;}}
  table {{width:100%; border-collapse:collapse; background:#fff;
          border-radius:8px; overflow:hidden; box-shadow:0 1px 4px rgba(0,0,0,.08);}}
  th {{background:#2d3748; color:#e2e8f0; padding:10px 14px; text-align:left; font-size:.82rem;}}
  tr.even td {{background:#f7fafc;}}
  tr.odd  td {{background:#fff;}}
  td {{padding:9px 14px; font-size:.84rem; border-bottom:1px solid #e2e8f0;}}
  .badge {{background:#e53e3e; color:#fff; border-radius:4px; padding:2px 7px; font-size:.73rem; font-weight:700;}}
  .ks-mets {{color:#6b46c1; font-style:italic;}}
  footer {{text-align:center; color:#a0aec0; font-size:.78rem; padding:20px;}}
</style>
</head>
<body>
<header>
  <h1>🔑 Keystone Species &amp; Dependency Report</h1>
  <p>Generated: {ts} &nbsp;|&nbsp; {len(mag_ids)} MAG(s) &nbsp;|&nbsp; {len(edges)} cross-feeding edge(s)</p>
</header>
<div class="content">
  <div class="summary">
    <strong>Community summary:</strong> {len(mag_ids)} MAGs detected.
    <strong>{n_keystone} Keystone species</strong> identified.
    Keystone MAGs are the sole suppliers of one or more critical metabolites.
  </div>
{alert_html}
  <table>
    <thead><tr>
      <th>MAG ID</th>
      <th title="Fraction of substrates sourced from other MAGs">Dependency Score ↑</th>
      <th title="Cross-feeding pairs provided per potential partner">Provision Score ↑</th>
      <th>Metabolites Provided</th>
      <th>Metabolites Received</th>
      <th>Keystone?</th>
      <th>Keystone Metabolites</th>
    </tr></thead>
    <tbody>
{table_rows}
    </tbody>
  </table>
</div>
<footer>Cross-Feeding Network Analyzer — {ts}</footer>
</body>
</html>
"""

    html_path = out_dir / "keystone_report.html"
    html_path.write_text(html_content, encoding="utf-8")
    log.info("  Saved keystone HTML report: %s", html_path)

    return {
        "rows": rows,
        "alerts": keystone_alerts,
        "n_keystone": n_keystone,
        "csv_path": csv_path,
        "html_path": html_path,
    }


# ---------------------------------------------------------------------------
# Master orchestrator
# ---------------------------------------------------------------------------

def run_crossfeed_analysis(
    table_path:  Path,
    out_dir:     Path,
    log:         logging.Logger | None = None,
    cache_path:  Path | None = None,
    mag_col:     str = "mag_id",
    ko_col:      str = "ko_id",
    max_api_kos: int = 500,
) -> dict[str, Any]:
    """
    Full cross-feeding network pipeline:
      load annotations → KEGG lookup → graph build → visualise → export.

    Parameters
    ----------
    table_path  : Path to MAG annotation CSV/TSV.
    out_dir     : Directory to write outputs.
    log         : Python logger (created if None).
    cache_path  : Path for KEGG JSON cache.
                  Defaults to ``<out_dir>/kegg_cache.json``.
    mag_col     : Name of the MAG ID column.
    ko_col      : Name of the KO ID column.
    max_api_kos : Maximum number of KOs to resolve via the KEGG API.

    Returns
    -------
    dict with keys:
        "edges"    – list of cross-feeding edge dicts
        "mag_sets" – per-MAG substrate/product sets
        "graph"    – NetworkX DiGraph (or None)
        "outputs"  – list of produced file Paths
    """
    if log is None:
        logging.basicConfig(level=logging.INFO,
                            format="%(asctime)s [%(levelname)s] %(message)s",
                            datefmt="%H:%M:%S")
        log = logging.getLogger("crossfeed_network")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if cache_path is None:
        cache_path = out_dir / "kegg_cache.json"

    log.info("Cross-Feeding Network Builder → %s", out_dir)
    log.info("  Annotation table : %s", table_path)
    log.info("  KEGG cache       : %s", cache_path)

    # Load cache
    cache = _load_cache(cache_path)

    # 1. Parse annotations
    mag_to_kos = load_mag_annotations(table_path, log, mag_col=mag_col, ko_col=ko_col)
    if not mag_to_kos:
        log.error("  No MAG annotations loaded – aborting.")
        return {"edges": [], "mag_sets": {}, "graph": None, "outputs": []}

    # Load KO descriptions from the annotation table (ko_id -> description)
    ko_descriptions: dict[str, str] = {}
    keystone_data: "dict | None" = None
    json_path: "Path | None" = None
    try:
        import csv as _csv_mod
        with open(table_path, newline="", encoding="utf-8") as _fh:
            _sniff = _fh.read(2048); _fh.seek(0)
            _sep = "\t" if "\t" in _sniff else ","
            for _row in _csv_mod.DictReader(_fh, delimiter=_sep):
                _norm = {k.strip().lower(): v.strip() for k, v in _row.items() if k}
                _ko = _norm.get(ko_col.lower(), "")
                _desc = _norm.get("description", "") or _norm.get("gene_name", "")
                if _ko and _desc and _ko not in ko_descriptions:
                    ko_descriptions[_ko] = _desc
    except Exception as _e:
        log.debug("Could not load ko_descriptions: %s", _e)

    # 2. Resolve KEGG compounds per MAG
    mag_sets = build_mag_metabolite_sets(mag_to_kos, cache, log,
                                         max_api_kos=max_api_kos)

    # 3. Find cross-feeding edges
    edges = find_cross_feeding_edges(mag_sets, cache, log)

    outputs: list[Path] = []
    G = None  # type: ignore[assignment]

    if not edges:
        log.warning("  No cross-feeding edges found.  "
                    "Check that KO terms map to valid KEGG reactions.")
    else:
        # 4. Build graph
        G = build_bipartite_graph(edges, mag_sets, log)

        # 5. Visualise
        html_path = plot_crossfeed_network(G, edges, mag_sets, out_dir, log)
        if html_path:
            outputs.append(html_path)

        # 6. Export CSV
        csv_path = export_edge_csv(edges, out_dir, log)
        outputs.append(csv_path)

        # 7. Export JSON summary
        json_path = export_summary_json(edges, mag_sets, G, out_dir, log)
        outputs.append(json_path)

        # 8. Sankey flow diagram
        sankey_path = plot_sankey_diagram(edges, mag_sets, out_dir, log)
        if sankey_path:
            outputs.append(sankey_path)

        # 9. Pathway completeness heatmap
        heatmap_path = plot_pathway_heatmap(
            mag_to_kos, out_dir, log, ko_descriptions=ko_descriptions
        )
        if heatmap_path:
            outputs.append(heatmap_path)

        # 10. Keystone & dependency report
        keystone_data = compute_keystone_report(edges, mag_sets, out_dir, log)
        if keystone_data:
            outputs.append(keystone_data["csv_path"])
            outputs.append(keystone_data["html_path"])

    # Save cache (updated with any new API results)
    _save_cache(cache, cache_path)
    log.info("  KEGG cache saved: %s  (%d entries)", cache_path, len(cache))

    log.info("Cross-feeding analysis complete.  %d output(s) written.", len(outputs))

    return {
        "edges":        edges,
        "mag_sets":     mag_sets,
        "graph":        G,
        "outputs":      outputs,
        "keystone_data": keystone_data,
    }


# ---------------------------------------------------------------------------
# Sample fixture generator
# ---------------------------------------------------------------------------

SAMPLE_ANNOTATIONS = """\
mag_id,ko_id,gene_id,description
MAG_001_Bifidobacterium,K01192,gene_001,Beta-galactosidase
MAG_001_Bifidobacterium,K01193,gene_002,Beta-fructosidase
MAG_001_Bifidobacterium,K00016,gene_003,L-lactate dehydrogenase
MAG_001_Bifidobacterium,K01834,gene_004,Phosphoglycerate mutase
MAG_001_Bifidobacterium,K00873,gene_005,Pyruvate kinase
MAG_002_Firmicutes,K00016,gene_006,L-lactate dehydrogenase (consumer)
MAG_002_Firmicutes,K00108,gene_007,Choline dehydrogenase
MAG_002_Firmicutes,K01640,gene_008,Homocitrate synthase
MAG_002_Firmicutes,K00382,gene_009,Dihydrolipoamide dehydrogenase
MAG_002_Firmicutes,K00242,gene_010,Succinate dehydrogenase
MAG_003_Akkermansia,K00382,gene_011,Dihydrolipoamide dehydrogenase
MAG_003_Akkermansia,K01190,gene_012,Beta-galactosidase
MAG_003_Akkermansia,K00600,gene_013,Serine hydroxymethyltransferase
MAG_003_Akkermansia,K01915,gene_014,Glutamine synthetase
MAG_003_Akkermansia,K01207,gene_015,Pullulanase
MAG_004_Bacteroides,K01186,gene_016,Neuraminidase
MAG_004_Bacteroides,K01190,gene_017,Beta-galactosidase
MAG_004_Bacteroides,K00626,gene_018,Acetyl-CoA acetyltransferase
MAG_004_Bacteroides,K01595,gene_019,Phosphoenolpyruvate carboxykinase
MAG_004_Bacteroides,K00239,gene_020,Succinate dehydrogenase subunit A
"""


def generate_sample_fixture(out_path: Path, log: logging.Logger) -> None:
    """Write a sample MAG annotation CSV to ``out_path``."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(SAMPLE_ANNOTATIONS, encoding="utf-8")
    log.info("Sample annotation fixture written to: %s", out_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Metabolic Hand-off & Cross-Feeding Network Builder for MAG inputs.\n\n"
            "Accepts a CSV/TSV annotation table with 'mag_id' and 'ko_id' columns, "
            "resolves KEGG reactions, identifies metabolic hand-offs between MAGs, "
            "and produces an interactive bipartite network visualisation."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--input", "-i",
        metavar="TABLE",
        help="Path to MAG annotation CSV/TSV (required unless --generate-sample).",
    )
    p.add_argument(
        "--out-dir", "-o",
        default="output/crossfeed",
        metavar="DIR",
        help="Output directory for generated files.",
    )
    p.add_argument(
        "--cache",
        metavar="FILE",
        help="Path to KEGG JSON cache file. Defaults to <out-dir>/kegg_cache.json.",
    )
    p.add_argument(
        "--mag-col",
        default="mag_id",
        metavar="COL",
        help="Name of the MAG ID column in the annotation table.",
    )
    p.add_argument(
        "--ko-col",
        default="ko_id",
        metavar="COL",
        help="Name of the KO ID column in the annotation table.",
    )
    p.add_argument(
        "--max-api-kos",
        type=int,
        default=500,
        metavar="N",
        help="Maximum number of unique KO terms to resolve via KEGG API.",
    )
    p.add_argument(
        "--generate-sample",
        metavar="FILE",
        help="Write a sample MAG annotation CSV to FILE and exit.",
    )
    p.add_argument(
        "--verbose", "-v",
        action="store_true",
        default=False,
        help="Enable DEBUG logging.",
    )
    return p.parse_args()


def main() -> None:
    args = _parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    log = logging.getLogger("crossfeed_network")

    # ── generate sample fixture ───────────────────────────────────────────────
    if args.generate_sample:
        generate_sample_fixture(Path(args.generate_sample), log)
        return

    if not args.input:
        log.error("--input is required (or use --generate-sample to create a fixture).")
        raise SystemExit(1)

    cache_path = Path(args.cache) if args.cache else None

    run_crossfeed_analysis(
        table_path=Path(args.input),
        out_dir=Path(args.out_dir),
        log=log,
        cache_path=cache_path,
        mag_col=args.mag_col,
        ko_col=args.ko_col,
        max_api_kos=args.max_api_kos,
    )


if __name__ == "__main__":
    main()
