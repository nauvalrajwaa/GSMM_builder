#!/usr/bin/env python3
"""
main.py
=======
Single-entry orchestrator for the GSMM workflow:

  fetch_ncbi.py -> generate_model.py -> visualize_model.py -> crossfeed_network.py

Run from repository root:

    python main.py --email you@example.com

This script is fail-fast by default. Use --continue-on-error to keep processing
remaining accessions if one fails.
"""

from __future__ import annotations

import argparse
import logging
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
SRC_DIR = SCRIPT_DIR / "src"

FETCH_SCRIPT = SRC_DIR / "fetch_ncbi.py"
GENERATE_SCRIPT = SRC_DIR / "generate_model.py"
VISUALIZE_SCRIPT = SRC_DIR / "visualize_model.py"
CROSSFEED_SCRIPT = SRC_DIR / "crossfeed_network.py"

DEFAULT_ACCESSIONS = ["NC_000913.3"]


@dataclass
class AccessionResult:
    accession: str
    ok: bool = True
    notes: list[str] = field(default_factory=list)


def _setup_logging(verbose: bool = False) -> logging.Logger:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    return logging.getLogger("main")


def _require_scripts() -> None:
    missing = [
        p for p in (FETCH_SCRIPT, GENERATE_SCRIPT, VISUALIZE_SCRIPT, CROSSFEED_SCRIPT)
        if not p.is_file()
    ]
    if missing:
        raise FileNotFoundError(
            "Required script(s) not found:\n" + "\n".join(f"- {m}" for m in missing)
        )


def _run(
    cmd: list[str],
    log: logging.Logger,
    label: str,
    timeout: int,
) -> tuple[bool, str]:
    log.info("[%s] $ %s", label, shlex.join(cmd))
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        log.error("[%s] timed out after %ds", label, timeout)
        return False, "TIMEOUT"
    except Exception as exc:
        log.error("[%s] exception: %s", label, exc)
        return False, str(exc)

    combined = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0:
        log.error("[%s] exited with code %d", label, proc.returncode)
        tail = combined.splitlines()[-30:]
        for line in tail:
            log.error("    %s", line)
        return False, combined

    if combined.strip():
        for line in combined.splitlines():
            log.debug("    %s", line)
    return True, combined


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Run full GSMM pipeline in one command (fetch -> model -> viz -> crossfeed).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--accessions",
        nargs="+",
        default=DEFAULT_ACCESSIONS,
        metavar="ACC",
        help="RefSeq accession(s) to process.",
    )
    p.add_argument(
        "--email",
        default="user@example.com",
        help="E-mail for NCBI Entrez (required by NCBI policy).",
    )
    p.add_argument(
        "--run-dir",
        default="output/all_in_one",
        metavar="DIR",
        help="Root output directory for this all-in-one run.",
    )
    p.add_argument(
        "--python",
        default=sys.executable,
        metavar="PYTHON",
        help="Python interpreter used to invoke src scripts.",
    )
    p.add_argument(
        "--gap-fill",
        default=None,
        metavar="MEDIUM",
        help="Gap-fill medium passed to generate_model.py (e.g. M9, LB).",
    )
    p.add_argument(
        "--max-api-kos",
        type=int,
        default=220,
        metavar="N",
        help="Maximum KO terms to resolve in crossfeed step.",
    )
    p.add_argument("--skip-download", action="store_true", help="Skip fetch_ncbi step.")
    p.add_argument("--skip-generate", action="store_true", help="Skip generate_model step.")
    p.add_argument("--skip-visualize", action="store_true", help="Skip visualize_model step.")
    p.add_argument("--skip-crossfeed", action="store_true", help="Skip crossfeed_network step.")
    p.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Continue remaining accessions if one fails.",
    )
    p.add_argument("--verbose", "-v", action="store_true", help="Enable DEBUG logging.")
    return p.parse_args()


def main() -> None:
    args = _parse_args()
    log = _setup_logging(args.verbose)

    try:
        _require_scripts()
    except Exception as exc:
        log.error("%s", exc)
        raise SystemExit(1)

    run_dir = Path(args.run_dir)
    fetch_dir = run_dir / "fetch"
    mag_table = run_dir / "mag_annotations.csv"
    kegg_cache = run_dir / "kegg_cache.json"

    run_dir.mkdir(parents=True, exist_ok=True)
    fetch_dir.mkdir(parents=True, exist_ok=True)

    accessions = list(dict.fromkeys(args.accessions))
    log.info("All-in-one GSMM run")
    log.info("Accessions : %s", ", ".join(accessions))
    log.info("Run dir    : %s", run_dir)
    log.info("Python     : %s", args.python)

    if not args.skip_download:
        fetch_cmd = [
            args.python,
            str(FETCH_SCRIPT),
            "--accessions",
            *accessions,
            "--out-dir",
            str(fetch_dir),
            "--email",
            args.email,
            "--mag-table",
            str(mag_table),
        ]
        ok, _ = _run(fetch_cmd, log, label="fetch_ncbi", timeout=3600)
        if not ok:
            raise SystemExit(1)
    else:
        log.info("[fetch_ncbi] skipped by --skip-download")

    results: list[AccessionResult] = []

    for acc in accessions:
        safe = acc.replace(".", "_")
        acc_dir = run_dir / safe
        acc_dir.mkdir(parents=True, exist_ok=True)

        fasta_path = fetch_dir / f"{acc}.fna"
        gbk_path = fetch_dir / f"{acc}.gbk"
        model_path = acc_dir / "output" / "model.xml"
        viz_dir = acc_dir / "output" / "viz"
        log_file = acc_dir / "pipeline.log"

        r = AccessionResult(accession=acc)

        if not args.skip_generate:
            if not fasta_path.is_file() or not gbk_path.is_file():
                msg = f"Missing fetched inputs for {acc}: {fasta_path.name}/{gbk_path.name}"
                log.error(msg)
                r.ok = False
                r.notes.append(msg)
                results.append(r)
                if args.continue_on_error:
                    continue
                raise SystemExit(1)

            gen_cmd = [
                args.python,
                str(GENERATE_SCRIPT),
                "--fasta",
                str(fasta_path),
                "--gbk",
                str(gbk_path),
                "--output",
                str(model_path),
                "--log-file",
                str(log_file),
            ]
            if args.gap_fill:
                gen_cmd += ["--gap-fill", args.gap_fill]

            ok, _ = _run(gen_cmd, log, label=f"generate_model:{acc}", timeout=7200)
            if not ok:
                r.ok = False
                r.notes.append("generate_model failed")
                results.append(r)
                if args.continue_on_error:
                    continue
                raise SystemExit(1)
        else:
            log.info("[generate_model:%s] skipped by --skip-generate", acc)

        if not args.skip_visualize:
            if not model_path.is_file():
                msg = f"Missing model.xml for {acc}: {model_path}"
                log.error(msg)
                r.ok = False
                r.notes.append(msg)
                results.append(r)
                if args.continue_on_error:
                    continue
                raise SystemExit(1)

            if not args.skip_crossfeed and mag_table.is_file():
                shutil.copy2(mag_table, acc_dir / "mag_annotations.csv")

            viz_cmd = [
                args.python,
                str(VISUALIZE_SCRIPT),
                "--model",
                str(model_path),
                "--out-dir",
                str(viz_dir),
            ]
            ok, _ = _run(viz_cmd, log, label=f"visualize_model:{acc}", timeout=2400)
            if not ok:
                r.ok = False
                r.notes.append("visualize_model failed")
                results.append(r)
                if args.continue_on_error:
                    continue
                raise SystemExit(1)
        else:
            log.info("[visualize_model:%s] skipped by --skip-visualize", acc)

        if not args.skip_crossfeed:
            if mag_table.is_file():
                crossfeed_dir = viz_dir / "crossfeed"
                crossfeed_html = crossfeed_dir / "crossfeed_network.html"

                should_run_crossfeed = args.skip_visualize or not crossfeed_html.is_file()

                if should_run_crossfeed:
                    cf_cmd = [
                        args.python,
                        str(CROSSFEED_SCRIPT),
                        "--input",
                        str(mag_table),
                        "--out-dir",
                        str(crossfeed_dir),
                        "--cache",
                        str(kegg_cache),
                        "--max-api-kos",
                        str(args.max_api_kos),
                    ]
                    ok, _ = _run(cf_cmd, log, label=f"crossfeed_network:{acc}", timeout=3600)
                    if not ok:
                        r.ok = False
                        r.notes.append("crossfeed_network failed")
                        results.append(r)
                        if args.continue_on_error:
                            continue
                        raise SystemExit(1)
                else:
                    log.info("[crossfeed_network:%s] generated during visualize_model auto-discovery", acc)
            else:
                msg = f"MAG table not found, skipping crossfeed: {mag_table}"
                log.warning("[%s] %s", acc, msg)
                r.notes.append(msg)
        else:
            log.info("[crossfeed_network:%s] skipped by --skip-crossfeed", acc)

        results.append(r)

    failed = [x for x in results if not x.ok]
    log.info("-" * 70)
    log.info("Run summary: %d/%d accession(s) succeeded", len(results) - len(failed), len(results))
    for x in results:
        status = "PASS" if x.ok else "FAIL"
        msg = f"{x.accession}: {status}"
        if x.notes:
            msg += f" | {'; '.join(x.notes)}"
        log.info(msg)
    log.info("Outputs root: %s", run_dir)
    log.info("-" * 70)

    raise SystemExit(1 if failed else 0)


if __name__ == "__main__":
    main()
