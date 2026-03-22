#!/usr/bin/env python3
"""
generate_model.py
=================
Automated pipeline to build a Genome-Scale Metabolic Model (GSMM)
from a FASTA protein/nucleotide file and a GenBank annotation file.

Workflow
--------
1. Validate input files (genomics.fasta, genomics.gbk)
2. Extract / translate CDS sequences from the GenBank file → proteins.faa
3. Run CarveMe to reconstruct a draft SBML model
4. Load model with COBRApy, run basic quality checks
5. Export final model as SBML (.xml)
6. (Optional) Generate comprehensive visualizations via visualize_model.py

Usage
-----
    python src/generate_model.py [--fasta genomics.fasta] [--gbk genomics.gbk]
                             [--output output/model.xml] [--gap-fill medium]
                             [--visualize] [--viz-dir output/viz]

Requirements
------------
    conda env create -f environment.yml
    conda activate gsmm_pipeline
"""

import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path

import colorlog
from Bio import SeqIO

# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------

def setup_logging(log_file: Path | None = None) -> logging.Logger:
    fmt = "%(log_color)s%(asctime)s [%(levelname)s] %(message)s%(reset)s"
    handler = colorlog.StreamHandler()
    handler.setFormatter(colorlog.ColoredFormatter(fmt, datefmt="%H:%M:%S"))
    handlers: list[logging.Handler] = [handler]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    logging.basicConfig(level=logging.INFO, handlers=handlers)
    return logging.getLogger("gsmm_pipeline")


# ---------------------------------------------------------------------------
# Step 1 – Validate inputs
# ---------------------------------------------------------------------------

def validate_inputs(fasta: Path, gbk: Path, log: logging.Logger) -> None:
    log.info("Validating input files …")
    missing = [p for p in (fasta, gbk) if not p.is_file()]
    if missing:
        for m in missing:
            log.error("Missing input file: %s", m)
        sys.exit(1)
    log.info("  FASTA : %s  (%d bytes)", fasta, fasta.stat().st_size)
    log.info("  GBK   : %s  (%d bytes)", gbk, gbk.stat().st_size)


# ---------------------------------------------------------------------------
# Step 2 – Extract protein sequences from GenBank
# ---------------------------------------------------------------------------

def extract_proteins(gbk: Path, out_faa: Path, log: logging.Logger) -> int:
    """
    Parse the GenBank file and write translated CDS sequences to a FASTA file.
    Returns the number of proteins extracted.
    """
    log.info("Extracting CDS protein sequences from GenBank …")
    records_written = 0
    with out_faa.open("w") as fh:
        for record in SeqIO.parse(str(gbk), "genbank"):
            for feat in record.features:
                if feat.type != "CDS":
                    continue
                # Prefer pre-translated sequence when available
                translation = feat.qualifiers.get("translation", [None])[0]
                locus_tag = feat.qualifiers.get("locus_tag", [f"CDS_{records_written}"])[0]
                gene = feat.qualifiers.get("gene", [""])[0]
                product = feat.qualifiers.get("product", ["hypothetical protein"])[0]

                if translation is None:
                    # Translate from nucleotide sequence
                    try:
                        nt_seq = feat.extract(record.seq)
                        translation = str(nt_seq.translate(to_stop=True))
                    except Exception as exc:
                        log.warning("  Could not translate %s: %s", locus_tag, exc)
                        continue

                if len(translation) < 10:
                    continue  # skip fragments

                header = f">{locus_tag} {gene} {product}"
                fh.write(f"{header.strip()}\n{translation}\n")
                records_written += 1

    log.info("  Proteins written: %d → %s", records_written, out_faa)
    if records_written == 0:
        log.error("No CDS features found in %s. Check the GenBank file.", gbk)
        sys.exit(1)
    return records_written


# ---------------------------------------------------------------------------
# Step 3 – Run CarveMe
# ---------------------------------------------------------------------------

def run_carveme(
    faa: Path,
    output_xml: Path,
    gap_fill: str | None,
    log: logging.Logger,
) -> None:
    log.info("Running CarveMe reconstruction …")
    # Resolve carve binary relative to the current Python interpreter so it
    # works regardless of whether the conda env bin dir is on PATH.
    carve_bin = str(Path(sys.executable).parent / "carve")
    cmd = [carve_bin, str(faa), "--output", str(output_xml)]
    if gap_fill:
        cmd += ["--gapfill", gap_fill]

    log.info("  Command: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True, env={**os.environ, "PATH": str(Path(sys.executable).parent) + ":" + os.environ.get("PATH", "")})

    if result.stdout:
        for line in result.stdout.splitlines():
            log.info("  [carve] %s", line)
    if result.stderr:
        for line in result.stderr.splitlines():
            log.warning("  [carve stderr] %s", line)

    if result.returncode != 0:
        log.error("CarveMe exited with code %d", result.returncode)
        sys.exit(result.returncode)

    if not output_xml.is_file():
        log.error("CarveMe did not produce output file: %s", output_xml)
        sys.exit(1)

    log.info("  CarveMe finished → %s", output_xml)


# ---------------------------------------------------------------------------
# Step 4 – QC with COBRApy
# ---------------------------------------------------------------------------

def qc_model(xml_path: Path, log: logging.Logger) -> tuple:
    """
    Load the SBML model and run basic quality checks.

    Returns
    -------
    tuple : (cobra.Model, list[str] blocked_rxn_ids, float|None fba_value)
    """
    log.info("Running COBRApy quality checks …")
    try:
        import cobra
        # Silence cobra's verbose 'Adding exchange reaction' / 'Ignoring reaction' messages
        logging.getLogger("cobra").setLevel(logging.ERROR)
        model = cobra.io.read_sbml_model(str(xml_path))
        logging.getLogger("cobra").setLevel(logging.WARNING)
    except Exception as exc:
        log.error("Failed to load SBML model: %s", exc)
        sys.exit(1)

    log.info("  Reactions   : %d", len(model.reactions))
    log.info("  Metabolites : %d", len(model.metabolites))
    log.info("  Genes       : %d", len(model.genes))

    # FBA
    fba_value = None
    try:
        solution = model.optimize()
        fba_value = solution.objective_value if solution.status == "optimal" else None
        log.info("  FBA growth  : %.6f  (status: %s)", solution.objective_value, solution.status)
        if solution.status != "optimal":
            log.warning("  FBA did not return an optimal solution – model may need gap-filling.")
    except Exception as exc:
        log.warning("  FBA failed: %s", exc)

    # Blocked reactions
    blocked: list[str] = []
    try:
        from cobra.flux_analysis import find_blocked_reactions
        blocked = find_blocked_reactions(model)
        pct = 100 * len(blocked) / max(len(model.reactions), 1)
        log.info("  Blocked rxns: %d / %d  (%.1f%%)", len(blocked), len(model.reactions), pct)
    except Exception as exc:
        log.warning("  Could not compute blocked reactions: %s", exc)

    return model, blocked, fba_value


# ---------------------------------------------------------------------------
# Step 5 – Export
# ---------------------------------------------------------------------------

def export_model(xml_path: Path, log: logging.Logger) -> None:
    """
    The primary output is the SBML file already written by CarveMe.
    Confirm it exists and log its size.
    """
    log.info("Exporting final model …")
    if xml_path.is_file():
        log.info("  SBML model  : %s  (%d bytes)", xml_path, xml_path.stat().st_size)
    else:
        log.error("Expected SBML output not found: %s", xml_path)
        sys.exit(1)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Automated GSMM pipeline using CarveMe + COBRApy",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--fasta", default="genomics.fasta", help="Input FASTA file")
    p.add_argument("--gbk",   default="genomics.gbk",   help="Input GenBank file")
    p.add_argument("--output", default="output/model.xml", help="Output SBML model path")
    p.add_argument(
        "--gap-fill",
        default=None,
        metavar="MEDIUM",
        help="Gap-fill on a medium, e.g. 'M9' or 'LB'. Pass the medium name accepted by CarveMe.",
    )
    p.add_argument("--log-file", default=None, help="Optional log file path")
    p.add_argument(
        "--visualize",
        action="store_true",
        default=False,
        help="Generate comprehensive visualizations after model reconstruction.",
    )
    p.add_argument(
        "--viz-dir",
        default=None,
        metavar="DIR",
        help="Directory for visualization outputs (default: <output-dir>/viz).",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    log = setup_logging(Path(args.log_file) if args.log_file else None)

    fasta  = Path(args.fasta)
    gbk    = Path(args.gbk)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    log.info("=" * 60)
    log.info("  GSMM Pipeline  –  CarveMe + COBRApy")
    log.info("=" * 60)

    # Intermediate file: extracted protein sequences
    proteins_faa = output.parent / "proteins.faa"

    # --- Pipeline steps ---
    validate_inputs(fasta, gbk, log)
    extract_proteins(gbk, proteins_faa, log)
    run_carveme(proteins_faa, output, args.gap_fill, log)
    model, blocked, fba_value = qc_model(output, log)
    export_model(output, log)

    # --- Optional visualizations ---
    viz_dir: Path | None = None
    if args.visualize:
        viz_dir = Path(args.viz_dir) if args.viz_dir else output.parent / "viz"
        log.info("=" * 60)
        log.info("  Generating visualizations → %s", viz_dir)
        log.info("=" * 60)
        try:
            from visualize_model import run_all_visualizations
            run_all_visualizations(
                model,
                blocked=blocked,
                fba_value=fba_value,
                out_dir=viz_dir,
                log=log,
            )
        except ImportError:
            log.error(
                "visualize_model.py not found or missing dependencies. "
                "Ensure it is in the same directory and the conda env is active."
            )

    log.info("=" * 60)
    log.info("  Pipeline complete!")
    log.info("  Model saved to: %s", output)
    if args.visualize:
        log.info("  Visualizations: %s", viz_dir)
    log.info("=" * 60)


if __name__ == "__main__":
    main()
