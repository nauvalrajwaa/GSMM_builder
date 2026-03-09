#!/usr/bin/env python3
"""
fetch_ncbi.py
=============
Download complete RefSeq genome records (FASTA nucleotide + GenBank annotation)
from NCBI for use as test input to the GSMM pipeline.

Features
--------
* Fetches nucleotide FASTA  (.fna) and GenBank annotation (.gbk) for any
  RefSeq accession (NC_*, NZ_*, etc.) using Biopython's Entrez interface.
* Retries with exponential back-off on transient NCBI failures.
* Verifies the downloaded files are non-empty and well-formed.
* Supports batch fetching of multiple accessions in one command.
* Writes a JSON manifest (ncbi_manifest.json) recording what was downloaded.

Built-in test accessions (E. coli K-12 only, as requested)
------------------------------------------------------------
    NC_000913.3   E. coli str. K-12 substr. MG1655   (complete chromosome)

Usage
-----
    # Fetch the default test set (E. coli K-12)
    python fetch_ncbi.py

    # Fetch one specific accession into a custom directory
    python fetch_ncbi.py --accessions NC_000913.3 --out-dir data/ecoli

    # Fetch multiple accessions with your NCBI email
    python fetch_ncbi.py --accessions NC_000913.3 NC_002695.2 \\
                         --email you@example.com --out-dir data/

    # Skip download if files already exist (idempotent)
    python fetch_ncbi.py --skip-existing

Requirements
------------
    biopython >= 1.81   (already in environment.yml)
    requests            (already in environment.yml via pip)
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any

# Biopython Entrez – always available in this environment
from Bio import Entrez, SeqIO

# ---------------------------------------------------------------------------
# Built-in test dataset
# ---------------------------------------------------------------------------

# Each entry: (accession, organism_name, strain, notes)
TEST_ACCESSIONS: list[dict[str, str]] = [
    {
        "accession": "NC_000913.3",
        "organism":  "Escherichia coli",
        "strain":    "K-12 substr. MG1655",
        "notes":     "Complete chromosome; NCBI RefSeq benchmark genome.",
    },
]

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def _setup_logging(verbose: bool = False) -> logging.Logger:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    return logging.getLogger("fetch_ncbi")


# ---------------------------------------------------------------------------
# Core download helpers
# ---------------------------------------------------------------------------

_RETRY_DELAYS = [5, 15, 30, 60]   # seconds between attempts


def _entrez_fetch(
    db: str,
    accession: str,
    rettype: str,
    retmode: str,
    log: logging.Logger,
) -> str:
    """
    Fetch a single record from NCBI Entrez with automatic retries.
    Returns the raw text content.
    """
    for attempt, delay in enumerate([0] + _RETRY_DELAYS, start=1):
        if delay:
            log.debug("  Retry %d/%d – waiting %ds …", attempt, len(_RETRY_DELAYS) + 1, delay)
            time.sleep(delay)
        try:
            handle = Entrez.efetch(
                db=db,
                id=accession,
                rettype=rettype,
                retmode=retmode,
            )
            data = handle.read()
            handle.close()
            if data:
                return data
            log.warning("  Empty response for %s (rettype=%s) – retrying …", accession, rettype)
        except Exception as exc:
            log.warning("  Entrez error (attempt %d): %s", attempt, exc)
    raise RuntimeError(
        f"Failed to fetch {accession} (rettype={rettype}) after {len(_RETRY_DELAYS)+1} attempts."
    )


def _write_file(content: str, path: Path, log: logging.Logger) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    log.debug("  Wrote %d bytes → %s", len(content.encode()), path)


def _verify_fasta(path: Path, log: logging.Logger) -> bool:
    """Return True if the file contains at least one valid FASTA record."""
    try:
        records = list(SeqIO.parse(str(path), "fasta"))
        if not records:
            log.error("  FASTA verification failed: no records in %s", path)
            return False
        log.info("  FASTA OK: %d records, first id = %s", len(records), records[0].id)
        return True
    except Exception as exc:
        log.error("  FASTA verification error: %s", exc)
        return False


def _verify_gbk(path: Path, log: logging.Logger) -> bool:
    """Return True if the file contains at least one valid GenBank record with CDS features."""
    try:
        records = list(SeqIO.parse(str(path), "genbank"))
        if not records:
            log.error("  GBK verification failed: no records in %s", path)
            return False
        cds_total = sum(
            1 for r in records for f in r.features if f.type == "CDS"
        )
        log.info(
            "  GBK OK: %d record(s), %d CDS feature(s), organism = %s",
            len(records),
            cds_total,
            records[0].annotations.get("organism", "unknown"),
        )
        return True
    except Exception as exc:
        log.error("  GBK verification error: %s", exc)
        return False


# ---------------------------------------------------------------------------
# Public API: fetch one accession
# ---------------------------------------------------------------------------

def fetch_genome(
    accession: str,
    out_dir: Path,
    email: str,
    log: logging.Logger,
    skip_existing: bool = False,
) -> dict[str, Any]:
    """
    Download the FASTA and GenBank files for a single RefSeq accession.

    Returns a result dict:
        {
            "accession": str,
            "fasta":     Path,
            "gbk":       Path,
            "skipped":   bool,
            "ok":        bool,
            "error":     str | None,
        }
    """
    Entrez.email = email
    out_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = out_dir / f"{accession}.fna"
    gbk_path   = out_dir / f"{accession}.gbk"

    result: dict[str, Any] = {
        "accession": accession,
        "fasta":     fasta_path,
        "gbk":       gbk_path,
        "skipped":   False,
        "ok":        False,
        "error":     None,
    }

    # --- skip if both files exist and look valid ---
    if skip_existing and fasta_path.is_file() and gbk_path.is_file():
        log.info("[%s] Files exist – skipping download.", accession)
        result["skipped"] = True
        result["ok"] = True
        return result

    log.info("[%s] Fetching FASTA (nucleotide) …", accession)
    try:
        fasta_text = _entrez_fetch("nucleotide", accession, "fasta", "text", log)
        _write_file(fasta_text, fasta_path, log)
        if not _verify_fasta(fasta_path, log):
            raise ValueError("FASTA verification failed")
    except Exception as exc:
        result["error"] = f"FASTA fetch failed: {exc}"
        log.error("[%s] %s", accession, result["error"])
        return result

    log.info("[%s] Fetching GenBank annotation …", accession)
    try:
        gbk_text = _entrez_fetch("nucleotide", accession, "gbwithparts", "text", log)
        _write_file(gbk_text, gbk_path, log)
        if not _verify_gbk(gbk_path, log):
            raise ValueError("GenBank verification failed")
    except Exception as exc:
        result["error"] = f"GenBank fetch failed: {exc}"
        log.error("[%s] %s", accession, result["error"])
        return result

    result["ok"] = True
    log.info("[%s] Download complete.  FASTA=%s  GBK=%s", accession, fasta_path, gbk_path)
    return result


# ---------------------------------------------------------------------------
# Batch fetch + manifest
# ---------------------------------------------------------------------------

def fetch_all(
    accessions: list[str],
    out_dir: Path,
    email: str,
    log: logging.Logger,
    skip_existing: bool = False,
) -> list[dict[str, Any]]:
    """Fetch multiple accessions and write a JSON manifest."""
    results = []
    for acc in accessions:
        log.info("=" * 55)
        res = fetch_genome(acc, out_dir, email, log, skip_existing=skip_existing)
        results.append(res)

    # Serialise Paths to strings for JSON
    manifest: list[dict[str, Any]] = []
    for r in results:
        manifest.append({
            "accession": r["accession"],
            "fasta":     str(r["fasta"]),
            "gbk":       str(r["gbk"]),
            "skipped":   r["skipped"],
            "ok":        r["ok"],
            "error":     r["error"],
        })

    manifest_path = out_dir / "ncbi_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    log.info("Manifest written → %s", manifest_path)

    # Summary
    n_ok  = sum(1 for r in results if r["ok"])
    n_err = len(results) - n_ok
    log.info("=" * 55)
    log.info("Fetch summary: %d OK, %d failed (total %d)", n_ok, n_err, len(results))
    if n_err:
        for r in results:
            if not r["ok"]:
                log.error("  FAILED: %s – %s", r["accession"], r["error"])

    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Download RefSeq FASTA + GenBank files from NCBI for GSMM pipeline testing."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--accessions",
        nargs="+",
        default=None,
        metavar="ACC",
        help=(
            "RefSeq accession(s) to download (e.g. NC_000913.3). "
            "Omit to use the built-in E. coli K-12 test set."
        ),
    )
    p.add_argument(
        "--out-dir",
        default="data/ncbi",
        metavar="DIR",
        help="Directory to store downloaded files.",
    )
    p.add_argument(
        "--email",
        default="user@example.com",
        help="E-mail address passed to NCBI Entrez (required by NCBI policy).",
    )
    p.add_argument(
        "--skip-existing",
        action="store_true",
        default=False,
        help="Skip download if output files already exist.",
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
    log  = _setup_logging(args.verbose)

    # Determine accession list
    if args.accessions:
        accessions = args.accessions
    else:
        accessions = [t["accession"] for t in TEST_ACCESSIONS]
        log.info("No accessions specified – using built-in test set:")
        for t in TEST_ACCESSIONS:
            log.info("  %s  %s %s", t["accession"], t["organism"], t["strain"])

    results = fetch_all(
        accessions  = accessions,
        out_dir     = Path(args.out_dir),
        email       = args.email,
        log         = log,
        skip_existing=args.skip_existing,
    )

    any_failed = any(not r["ok"] for r in results)
    sys.exit(1 if any_failed else 0)


if __name__ == "__main__":
    main()
