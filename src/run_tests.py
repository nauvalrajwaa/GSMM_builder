#!/usr/bin/env python3
"""
run_tests.py
============
End-to-end test runner for the GSMM pipeline.

For each configured test case the runner:
  1. Downloads the genome from NCBI (FASTA + GenBank) via fetch_ncbi.py
     – skips the download if the files already exist.
  2. Runs generate_model.py  (CarveMe reconstruction + COBRApy QC).
  3. Runs visualize_model.py (comprehensive plots and interactive maps).
  4. Records pass / fail for every step and prints a final summary table.

Test cases (E. coli K-12 focused, as requested)
------------------------------------------------
  NC_000913.3   E. coli K-12 substr. MG1655   complete chromosome

Usage
-----
    # Run all tests (downloads data automatically)
    python src/run_tests.py

    # Skip download if data already exists
    python src/run_tests.py --skip-download

    # Run only a specific accession
    python src/run_tests.py --accessions NC_000913.3

    # Keep intermediate files (proteins.faa, logs) after run
    python src/run_tests.py --keep-intermediates

    # Run with gap-filling on M9 medium
    python src/run_tests.py --gap-fill M9

    # Set NCBI email (required by NCBI policy)
    python src/run_tests.py --email your@email.com

Requirements
------------
    conda activate gsmm_pipeline
    # run from repository root:
    # python src/run_tests.py
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from mag_annotation_builder import build_mag_annotations_from_gbks


SCRIPT_DIR = Path(__file__).resolve().parent
FETCH_SCRIPT = SCRIPT_DIR / "fetch_ncbi.py"
GENERATE_SCRIPT = SCRIPT_DIR / "generate_model.py"
VISUALIZE_SCRIPT = SCRIPT_DIR / "visualize_model.py"
CROSSFEED_SCRIPT = SCRIPT_DIR / "crossfeed_network.py"


# ---------------------------------------------------------------------------
# Resolve the Python interpreter to use for sub-scripts.
# Priority:
#   1. GSMM_PYTHON env var (user override)
#   2. gsmm_pipeline conda env  (auto-discovered from CONDA_PREFIX or common paths)
#   3. sys.executable fallback  (works if already inside gsmm_pipeline env)
# ---------------------------------------------------------------------------

def _find_pipeline_python() -> str:
    """Return the absolute path to the Python interpreter to use for sub-scripts."""
    # 1. Explicit user override
    override = os.environ.get("GSMM_PYTHON")
    if override and Path(override).is_file():
        return override

    # 2. Already inside the right conda env?
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    if "gsmm_pipeline" in conda_prefix:
        candidate = Path(conda_prefix) / "bin" / "python"
        if candidate.is_file():
            return str(candidate)

    # 3. Scan common conda base locations for gsmm_pipeline env
    conda_base_candidates = [
        Path.home() / "miniforge3",
        Path.home() / "miniconda3",
        Path.home() / "anaconda3",
        Path.home() / "opt" / "miniforge3",
        Path.home() / "opt" / "miniconda3",
        Path("/opt/conda"),
        Path("/opt/miniforge3"),
        Path("/opt/miniconda3"),
    ]
    # Also try the base of the current conda prefix
    if conda_prefix:
        cp = Path(conda_prefix)
        # e.g. /Users/user/miniforge3/envs/some_env → /Users/user/miniforge3
        if "envs" in cp.parts:
            idx = cp.parts.index("envs")
            conda_base_candidates.insert(0, Path(*cp.parts[:idx]))

    for base in conda_base_candidates:
        candidate = base / "envs" / "gsmm_pipeline" / "bin" / "python"
        if candidate.is_file():
            return str(candidate)

    # 4. Fall back to the current interpreter (may work if env already active)
    return sys.executable

# ---------------------------------------------------------------------------
# Test registry
# ---------------------------------------------------------------------------

@dataclass
class TestCase:
    accession:  str
    organism:   str
    strain:     str
    gap_fill:   str | None = None   # optional CarveMe gap-fill medium

TEST_CASES: list[TestCase] = [
    TestCase(
        accession = "NC_000913.3",
        organism  = "Escherichia coli",
        strain    = "K-12 substr. MG1655",
    ),
]

# ---------------------------------------------------------------------------
# Step result tracking
# ---------------------------------------------------------------------------

@dataclass
class StepResult:
    name:    str
    ok:      bool
    elapsed: float = 0.0
    notes:   str   = ""


@dataclass
class TestResult:
    case:    TestCase
    steps:   list[StepResult] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return all(s.ok for s in self.steps)

    def add(self, name: str, ok: bool, elapsed: float = 0.0, notes: str = "") -> None:
        self.steps.append(StepResult(name, ok, elapsed, notes))


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
    return logging.getLogger("run_tests")


# ---------------------------------------------------------------------------
# Helper: run a subprocess and return (ok, elapsed, stdout+stderr)
# ---------------------------------------------------------------------------

def _run(
    cmd: list[str],
    log: logging.Logger,
    label: str,
    timeout: int = 3600,
) -> tuple[bool, float, str]:
    log.info("  $ %s", " ".join(cmd))
    t0 = time.monotonic()
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        elapsed = time.monotonic() - t0
        combined = proc.stdout + proc.stderr
        if proc.returncode != 0:
            log.error("  [%s] exited with code %d", label, proc.returncode)
            for line in combined.splitlines()[-30:]:   # last 30 lines on error
                log.error("    %s", line)
            return False, elapsed, combined
        for line in combined.splitlines():
            log.debug("    %s", line)
        return True, elapsed, combined
    except subprocess.TimeoutExpired:
        elapsed = time.monotonic() - t0
        log.error("  [%s] timed out after %ds", label, timeout)
        return False, elapsed, "TIMEOUT"
    except Exception as exc:
        elapsed = time.monotonic() - t0
        log.error("  [%s] exception: %s", label, exc)
        return False, elapsed, str(exc)


# ---------------------------------------------------------------------------
# Per-accession test logic
# ---------------------------------------------------------------------------

def run_test_case(
    case: TestCase,
    base_dir: Path,
    email: str,
    log: logging.Logger,
    skip_download: bool = False,
    gap_fill: str | None = None,
    keep_intermediates: bool = False,
    community_fetch_dir: Path | None = None,
    community_mag_table: Path | None = None,
) -> TestResult:
    result = TestResult(case=case)

    acc     = case.accession
    safe    = acc.replace(".", "_")
    work    = base_dir / safe
    work.mkdir(parents=True, exist_ok=True)

    fasta_path = work / f"{acc}.fna"
    gbk_path   = work / f"{acc}.gbk"
    model_path = work / "output" / "model.xml"
    viz_dir    = work / "output" / "viz"
    log_path   = work / "pipeline.log"

    log.info("=" * 60)
    log.info("Test case: %s  %s %s", acc, case.organism, case.strain)
    log.info("Work dir : %s", work)
    log.info("=" * 60)

    pipeline_python = _find_pipeline_python()
    log.debug("Using Python interpreter: %s", pipeline_python)
    runtime_mag_table = base_dir / "community_mag_annotations_runtime.csv"

    # ------------------------------------------------------------------
    # Step 1 – Download from NCBI
    # ------------------------------------------------------------------
    step_name = "fetch_ncbi"
    if skip_download and fasta_path.is_file() and gbk_path.is_file():
        log.info("[%s] Skipping download (--skip-download, files exist).", acc)
        result.add(step_name, ok=True, notes="skipped (files exist)")
    elif community_fetch_dir is not None:
        src_fasta = community_fetch_dir / f"{acc}.fna"
        src_gbk = community_fetch_dir / f"{acc}.gbk"
        if src_fasta.is_file() and src_gbk.is_file():
            t0 = time.monotonic()
            shutil.copy2(src_fasta, fasta_path)
            shutil.copy2(src_gbk, gbk_path)
            elapsed = time.monotonic() - t0
            result.add(step_name, ok=True, elapsed=elapsed, notes="copied from shared accession fetch")
        else:
            log.info("[%s] Shared fetch files missing; falling back to direct fetch.", acc)
            fetch_cmd = [
                pipeline_python, str(FETCH_SCRIPT),
                "--accessions", acc,
                "--out-dir", str(work),
                "--email", email,
            ]
            if skip_download:
                fetch_cmd.append("--skip-existing")
            ok, elapsed, _ = _run(fetch_cmd, log, label=step_name)
            result.add(step_name, ok=ok, elapsed=elapsed)
            if not ok:
                log.error("[%s] Download failed – aborting test case.", acc)
                return result
    else:
        log.info("[%s] Step 1/3: Downloading genome from NCBI …", acc)
        fetch_cmd = [
            pipeline_python, str(FETCH_SCRIPT),
            "--accessions", acc,
            "--out-dir", str(work),
            "--email", email,
        ]
        if skip_download:
            fetch_cmd.append("--skip-existing")
        ok, elapsed, _ = _run(
            fetch_cmd,
            log,
            label=step_name,
        )
        result.add(step_name, ok=ok, elapsed=elapsed)
        if not ok:
            log.error("[%s] Download failed – aborting test case.", acc)
            return result

    # Verify downloaded files exist
    if not (fasta_path.is_file() and gbk_path.is_file()):
        log.error("[%s] Expected files missing after download step.", acc)
        result.add("file_check", ok=False, notes="FASTA or GBK missing")
        return result

    # ------------------------------------------------------------------
    # Step 2 – Run GSMM pipeline (generate_model.py)
    # ------------------------------------------------------------------
    log.info("[%s] Step 2/3: Running GSMM pipeline …", acc)
    pipeline_cmd = [
        pipeline_python, str(GENERATE_SCRIPT),
        "--fasta",  str(fasta_path),
        "--gbk",    str(gbk_path),
        "--output", str(model_path),
        "--log-file", str(log_path),
    ]
    gf = gap_fill or case.gap_fill
    if gf:
        pipeline_cmd += ["--gap-fill", gf]

    ok, elapsed, output_text = _run(
        pipeline_cmd,
        log,
        label="generate_model",
        timeout=7200,   # CarveMe can be slow on large genomes
    )

    # Parse key metrics from pipeline log
    notes = ""
    for keyword in ("Reactions", "Metabolites", "Genes", "FBA growth", "Blocked rxns"):
        for line in output_text.splitlines():
            if keyword in line:
                notes += line.strip() + "  |  "
                break

    result.add("generate_model", ok=ok, elapsed=elapsed, notes=notes.rstrip(" |"))
    if not ok:
        log.error("[%s] Pipeline failed – skipping visualization.", acc)
        return result

    if not model_path.is_file():
        result.add("model_file_check", ok=False, notes="model.xml not found")
        return result

    # ------------------------------------------------------------------
    # Step 3 – Visualize
    # ------------------------------------------------------------------
    effective_mag_table = community_mag_table if (community_mag_table is not None and community_mag_table.is_file()) else None
    if effective_mag_table is None:
        dynamic_gbks = sorted(base_dir.glob("*/*.gbk"))
        if dynamic_gbks:
            stats_dynamic = build_mag_annotations_from_gbks(dynamic_gbks, runtime_mag_table, log=log)
            if stats_dynamic["n_rows"] > 0:
                effective_mag_table = runtime_mag_table

    if effective_mag_table is not None and effective_mag_table.is_file():
        shutil.copy2(effective_mag_table, work / "mag_annotations.csv")

    log.info("[%s] Step 3/3: Generating visualizations …", acc)
    ok, elapsed, _ = _run(
        [
            pipeline_python, str(VISUALIZE_SCRIPT),
            "--model",   str(model_path),
            "--out-dir", str(viz_dir),
        ],
        log,
        label="visualize_model",
        timeout=600,
    )
    # Check expected output files.
    # "required" files must always be present; "optional" files depend on
    # whether the model has the relevant annotations (e.g. subsystem info).
    required = [
        "summary_stats.png",
        "compartment_breakdown.png",
        "degree_distribution.png",
        "flux_distribution.png",
        "reaction_network.html",
        "dashboard.html",
        "mutation_analysis.html",
        "index.html",
        "model_summary.csv",
    ]
    optional = [
        "subsystem_barchart.png",   # only produced when model has subsystem annotations
        "escher_map.html",          # only produced when escher map download succeeds
    ]
    present_req = [f for f in required if (viz_dir / f).is_file()]
    missing_req = [f for f in required if f not in present_req]
    present_opt = [f for f in optional if (viz_dir / f).is_file()]
    total_present = len(present_req) + len(present_opt)
    total_checked = len(required) + len(optional)
    notes_v = f"{total_present}/{total_checked} output files present"
    if missing_req:
        notes_v += f"  |  missing (required): {', '.join(missing_req)}"
    if present_opt:
        notes_v += f"  |  optional present: {', '.join(present_opt)}"

    result.add("visualize_model", ok=ok and not missing_req, elapsed=elapsed, notes=notes_v)

    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    mag_table = work / "mag_annotations.csv"
    if mag_table.is_file():
        log.info("[%s] Step 4/4: Running cross-feeding network analysis …", acc)
        crossfeed_dir = viz_dir / "crossfeed"
        ok_cf, elapsed_cf, _ = _run(
            [
                pipeline_python, str(CROSSFEED_SCRIPT),
                "--input",   str(mag_table),
                "--out-dir", str(crossfeed_dir),
                "--max-api-kos", "220",
            ],
            log,
            label="crossfeed_network",
            timeout=1800,
        )
        cf_html    = crossfeed_dir / "crossfeed_network.html"
        cf_sankey  = crossfeed_dir / "crossfeed_sankey.html"
        cf_heatmap = crossfeed_dir / "crossfeed_heatmap.html"
        cf_report  = crossfeed_dir / "keystone_report.html"
        cf_csv     = crossfeed_dir / "crossfeed_edges.csv"
        cf_json    = crossfeed_dir / "crossfeed_summary.json"
        cf_ks_csv  = crossfeed_dir / "keystone_report.csv"
        _cf_all    = [cf_html, cf_sankey, cf_heatmap, cf_report,
                      cf_csv, cf_json, cf_ks_csv]
        cf_outputs = [f.name for f in _cf_all if f.is_file()]
        notes_cf = f"{len(cf_outputs)}/7 output files present"
        if cf_outputs:
            notes_cf += f"  |  {', '.join(cf_outputs)}"
        result.add("crossfeed_network", ok=ok_cf, elapsed=elapsed_cf, notes=notes_cf)
    else:
        log.info("[%s] mag_annotations.csv not found – skipping crossfeed step.", acc)
        result.add("crossfeed_network", ok=True, notes="skipped (no MAG annotation table)")

    # ------------------------------------------------------------------
    # Cleanup (optional)
    # ------------------------------------------------------------------
    if not keep_intermediates:
        proteins_faa = work / "output" / "proteins.faa"
        if proteins_faa.is_file():
            proteins_faa.unlink()
        log.debug("[%s] Intermediate proteins.faa removed.", acc)

    return result


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

def _print_summary(results: list[TestResult], log: logging.Logger) -> None:
    log.info("")
    log.info("=" * 70)
    log.info("  TEST SUMMARY")
    log.info("=" * 70)

    col_acc  = max(len(r.case.accession) for r in results) + 2
    col_step = 20
    col_stat = 8
    col_time = 8

    header = (
        f"{'Accession':<{col_acc}}"
        f"{'Step':<{col_step}}"
        f"{'Status':<{col_stat}}"
        f"{'Elapsed':>{col_time}}  Notes"
    )
    log.info(header)
    log.info("-" * 70)

    for tr in results:
        acc = tr.case.accession
        for s in tr.steps:
            status = "PASS" if s.ok else "FAIL"
            log.info(
                "%-*s %-*s %-*s %*s  %s",
                col_acc,  acc,
                col_step, s.name,
                col_stat, status,
                col_time, f"{s.elapsed:.1f}s",
                s.notes[:70],
            )
            acc = ""   # only print accession on first row per test case
        if not tr.steps:
            log.info("%-*s  (no steps recorded)", col_acc, tr.case.accession)

    log.info("-" * 70)
    n_pass = sum(1 for r in results if r.ok)
    n_fail = len(results) - n_pass
    log.info("Overall: %d / %d test cases PASSED", n_pass, len(results))
    if n_fail:
        log.error("%d test case(s) FAILED.", n_fail)
    log.info("=" * 70)


# ---------------------------------------------------------------------------
# Write JSON report
# ---------------------------------------------------------------------------

def _write_report(results: list[TestResult], out_path: Path) -> None:
    report: list[dict[str, Any]] = []
    for tr in results:
        report.append({
            "accession": tr.case.accession,
            "organism":  tr.case.organism,
            "strain":    tr.case.strain,
            "passed":    tr.ok,
            "steps": [
                {
                    "name":    s.name,
                    "ok":      s.ok,
                    "elapsed": round(s.elapsed, 2),
                    "notes":   s.notes,
                }
                for s in tr.steps
            ],
        })
    out_path.write_text(json.dumps(report, indent=2), encoding="utf-8")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="End-to-end GSMM pipeline tester with NCBI data fetching.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--accessions",
        nargs="+",
        default=None,
        metavar="ACC",
        help="RefSeq accession(s) to test. Defaults to built-in E. coli K-12 set.",
    )
    p.add_argument(
        "--base-dir",
        default="tests",
        metavar="DIR",
        help="Root directory for per-accession work directories.",
    )
    p.add_argument(
        "--email",
        default="user@example.com",
        help="E-mail for NCBI Entrez (required by NCBI policy).",
    )
    p.add_argument(
        "--skip-download",
        action="store_true",
        default=False,
        help="Skip NCBI download if FASTA + GBK files already exist.",
    )
    p.add_argument(
        "--gap-fill",
        default=None,
        metavar="MEDIUM",
        help="Gap-fill medium passed to CarveMe (e.g. M9).",
    )
    p.add_argument(
        "--keep-intermediates",
        action="store_true",
        default=False,
        help="Keep intermediate files (proteins.faa) after each test.",
    )
    p.add_argument(
        "--report",
        default="tests/test_report.json",
        metavar="FILE",
        help="Path for JSON test report.",
    )
    p.add_argument(
        "--verbose", "-v",
        action="store_true",
        default=False,
        help="Enable DEBUG logging.",
    )
    return p.parse_args()


def main() -> None:
    args    = _parse_args()
    log     = _setup_logging(args.verbose)
    base    = Path(args.base_dir)
    report  = Path(args.report)

    # Resolve test cases
    if args.accessions:
        # Build minimal TestCase objects for user-supplied accessions
        test_cases = [
            TestCase(accession=acc, organism="user-supplied", strain="")
            for acc in args.accessions
        ]
    else:
        test_cases = TEST_CASES

    log.info("GSMM pipeline test runner")
    log.info("Test cases  : %d", len(test_cases))
    log.info("Base dir    : %s", base)
    log.info("Gap-fill    : %s", args.gap_fill or "none")
    log.info("Skip dl     : %s", args.skip_download)

    community_fetch_dir = base / "_community_fetch"
    community_mag_table = base / "community_mag_annotations.csv"
    accessions = [tc.accession for tc in test_cases]

    if not args.skip_download:
        log.info("Preparing shared fetch for all accessions (%d) …", len(accessions))
        pipeline_python = _find_pipeline_python()
        ok_shared, elapsed_shared, _ = _run(
            [
                pipeline_python,
                str(FETCH_SCRIPT),
                "--accessions",
                *accessions,
                "--out-dir",
                str(community_fetch_dir),
                "--email",
                args.email,
                "--mag-table",
                str(community_mag_table),
            ],
            log,
            label="fetch_ncbi_shared",
            timeout=3600,
        )
        if not ok_shared:
            log.warning("Shared fetch failed after %.1fs; per-case fetch fallback will be used.", elapsed_shared)

    gbk_candidates = [community_fetch_dir / f"{acc}.gbk" for acc in accessions if (community_fetch_dir / f"{acc}.gbk").is_file()]
    if not gbk_candidates:
        for acc in accessions:
            safe = acc.replace(".", "_")
            candidate = base / safe / f"{acc}.gbk"
            if candidate.is_file():
                gbk_candidates.append(candidate)

    if gbk_candidates:
        stats = build_mag_annotations_from_gbks(gbk_candidates, community_mag_table, log=log)
        log.info(
            "Community MAG table prepared: %s (rows=%d, mags=%d)",
            community_mag_table,
            stats["n_rows"],
            stats["n_mags_with_annotations"],
        )
    else:
        log.warning("No GenBank files available to build community MAG annotation table.")

    results: list[TestResult] = []
    for case in test_cases:
        tr = run_test_case(
            case               = case,
            base_dir           = base,
            email              = args.email,
            log                = log,
            skip_download      = args.skip_download,
            gap_fill           = args.gap_fill,
            keep_intermediates = args.keep_intermediates,
            community_fetch_dir = community_fetch_dir if community_fetch_dir.is_dir() else None,
            community_mag_table = community_mag_table if community_mag_table.is_file() else None,
        )
        results.append(tr)

    _print_summary(results, log)

    report.parent.mkdir(parents=True, exist_ok=True)
    _write_report(results, report)
    log.info("JSON report → %s", report)

    any_failed = any(not r.ok for r in results)
    sys.exit(1 if any_failed else 0)


if __name__ == "__main__":
    main()
