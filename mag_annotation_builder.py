#!/usr/bin/env python3
"""
mag_annotation_builder.py
=========================
Build MAG annotation tables from downloaded GenBank files.

This module extracts KO/EC identifiers from CDS features and writes a
crossfeed-compatible annotation table with columns:

    mag_id,ko_id,gene_id,description,accession,organism,source_gbk

The output can be consumed directly by crossfeed_network.py.
"""

from __future__ import annotations

import csv
import logging
import re
from pathlib import Path
from typing import Any

from Bio import SeqIO

_KO_TOKEN = re.compile(r"K\d{5}", re.IGNORECASE)
_EC_TOKEN = re.compile(r"^\d+\.\d+\.\d+\.(?:\d+|-)$")


def _slug(text: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")
    return cleaned or "Unknown"


def _normalise_ec(raw: str) -> str | None:
    token = raw.strip().upper().replace(" ", "")
    if token.startswith("EC:"):
        token = token[3:]
    token = token.rstrip(";,")
    if _EC_TOKEN.match(token):
        return f"EC:{token}"
    return None


def _extract_identifiers(cds_qualifiers: dict[str, list[str]]) -> list[str]:
    identifiers: list[str] = []
    seen: set[str] = set()

    for xref in cds_qualifiers.get("db_xref", []):
        for match in _KO_TOKEN.findall(str(xref)):
            ko = match.upper()
            if ko not in seen:
                seen.add(ko)
                identifiers.append(ko)

    for ec_raw in cds_qualifiers.get("EC_number", []):
        parts = re.split(r"[;,\s]+", str(ec_raw))
        for part in parts:
            ec = _normalise_ec(part)
            if ec and ec not in seen:
                seen.add(ec)
                identifiers.append(ec)

    return identifiers


def _derive_mag_id(index: int, accession: str, organism: str) -> str:
    genus = organism.split()[0] if organism else "Unknown"
    return f"MAG_{index:03d}_{_slug(genus)}_{_slug(accession)}"


def build_mag_annotations_from_gbks(
    gbk_paths: list[Path],
    out_csv: Path,
    log: logging.Logger | None = None,
    max_unique_ids_per_mag: int = 220,
) -> dict[str, Any]:
    if log is None:
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s [%(levelname)s] %(message)s",
            datefmt="%H:%M:%S",
        )
        log = logging.getLogger("mag_annotation_builder")

    unique_paths: list[Path] = []
    seen_paths: set[Path] = set()
    for p in gbk_paths:
        rp = Path(p).resolve()
        if rp in seen_paths:
            continue
        seen_paths.add(rp)
        unique_paths.append(rp)

    rows: list[dict[str, str]] = []
    cds_total = 0
    ko_rows = 0
    ec_rows = 0
    mags_with_rows = 0

    for i, gbk_path in enumerate(unique_paths, start=1):
        if not gbk_path.is_file():
            log.warning("GBK not found, skipping: %s", gbk_path)
            continue

        records = list(SeqIO.parse(str(gbk_path), "genbank"))
        if not records:
            log.warning("No GenBank records in: %s", gbk_path)
            continue

        first_record = records[0]
        accession = str(first_record.id or gbk_path.stem)
        organism = str(first_record.annotations.get("organism", "unknown"))
        mag_id = _derive_mag_id(i, accession, organism)
        selected_ids: set[str] = set()

        before_count = len(rows)
        cds_index = 0
        for record in records:
            record_accession = str(record.id or accession)
            for feat in record.features:
                if feat.type != "CDS":
                    continue
                cds_total += 1
                cds_index += 1
                ids = _extract_identifiers(feat.qualifiers)
                if not ids:
                    continue

                gene_id = (
                    feat.qualifiers.get("locus_tag", [None])[0]
                    or feat.qualifiers.get("protein_id", [None])[0]
                    or feat.qualifiers.get("gene", [None])[0]
                    or f"{record_accession}_CDS_{cds_index:05d}"
                )
                description = feat.qualifiers.get("product", ["hypothetical protein"])[0]

                for identifier in ids:
                    if identifier not in selected_ids and len(selected_ids) >= max_unique_ids_per_mag:
                        continue
                    selected_ids.add(identifier)
                    if identifier.startswith("K"):
                        ko_rows += 1
                    elif identifier.startswith("EC:"):
                        ec_rows += 1
                    rows.append(
                        {
                            "mag_id": mag_id,
                            "ko_id": identifier,
                            "gene_id": str(gene_id),
                            "description": str(description),
                            "accession": accession,
                            "organism": organism,
                            "source_gbk": str(gbk_path),
                        }
                    )

        if len(rows) > before_count:
            mags_with_rows += 1

    deduped_rows: list[dict[str, str]] = []
    seen_keys: set[tuple[str, str, str]] = set()
    for row in rows:
        key = (row["mag_id"], row["ko_id"], row["gene_id"])
        if key in seen_keys:
            continue
        seen_keys.add(key)
        deduped_rows.append(row)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "mag_id",
                "ko_id",
                "gene_id",
                "description",
                "accession",
                "organism",
                "source_gbk",
            ],
        )
        writer.writeheader()
        writer.writerows(deduped_rows)

    stats = {
        "out_csv": str(out_csv),
        "n_input_gbks": len(unique_paths),
        "n_rows": len(deduped_rows),
        "n_cds": cds_total,
        "n_mags_with_annotations": mags_with_rows,
        "n_ko_rows": ko_rows,
        "n_ec_rows": ec_rows,
    }
    log.info(
        "MAG annotation table written: %s (rows=%d, MAGs=%d, KO rows=%d, EC rows=%d)",
        out_csv,
        stats["n_rows"],
        stats["n_mags_with_annotations"],
        stats["n_ko_rows"],
        stats["n_ec_rows"],
    )
    return stats
