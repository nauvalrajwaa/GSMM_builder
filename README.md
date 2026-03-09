# GSMM Pipeline

Automated pipeline for generating a **Genome-Scale Metabolic Model (GSMM)**
from a FASTA + GenBank genomic input using **CarveMe** and **COBRApy**,
with optional **comprehensive visualizations**.

---

## Requirements

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda
- GLPK solver (bundled via conda) **or** CPLEX (optional, faster)

---

## Setup

```bash
# 1. Create and activate the conda environment
conda env create -f environment.yml
conda activate gsmm_pipeline

# 2. Download the CarveMe universal reference model (one-time)
carve --download
```

---

## Input files

| File | Description |
|------|-------------|
| `genomics.fasta` | Genomic nucleotide or protein FASTA |
| `genomics.gbk` | Annotated GenBank file (CDS features required) |

Place both files in the project root (same directory as `generate_model.py`).

---

## Usage

### Minimal (default paths, no gap-filling)
```bash
python generate_model.py
```

### Custom paths + gap-filling on M9 minimal medium
```bash
python generate_model.py \
    --fasta  genomics.fasta \
    --gbk    genomics.gbk \
    --output output/my_model.xml \
    --gap-fill M9
```

### With a log file
```bash
python generate_model.py --log-file pipeline.log
```

### With comprehensive visualizations
```bash
python generate_model.py --visualize
```

### Custom visualization output directory
```bash
python generate_model.py --visualize --viz-dir output/my_viz
```

### All options
```
usage: generate_model.py [-h] [--fasta FASTA] [--gbk GBK]
                         [--output OUTPUT] [--gap-fill MEDIUM]
                         [--log-file LOG_FILE]
                         [--visualize] [--viz-dir DIR]

  --fasta      Input FASTA file          (default: genomics.fasta)
  --gbk        Input GenBank file        (default: genomics.gbk)
  --output     Output SBML model path    (default: output/model.xml)
  --gap-fill   Gap-fill medium name      (e.g. M9, LB)
  --log-file   Save log to file
  --visualize  Generate visualizations after reconstruction
  --viz-dir    Directory for viz output  (default: output/viz)
```

---

## Pipeline steps

```
1. validate_inputs   – check that both input files exist and are non-empty
2. extract_proteins  – parse GenBank CDS features → output/proteins.faa
3. run_carveme       – reconstruct draft SBML model via CarveMe
4. qc_model          – load model with COBRApy; report reactions /
                       metabolites / genes, run FBA, count blocked reactions
5. export_model      – confirm final SBML file at --output path
6. visualize (opt.)  – generate all plots and interactive maps (--visualize)
```

---

## Output

| File | Description |
|------|-------------|
| `output/model.xml` | Final GSMM in SBML Level 3 format |
| `output/proteins.faa` | Intermediate protein FASTA used by CarveMe |
| `pipeline.log` | Full run log (if `--log-file` is used) |

### Visualization outputs (when `--visualize` is used)

All written to `output/viz/` by default (override with `--viz-dir`).

| File | Description |
|------|-------------|
| `summary_stats.png` | Bar chart of reactions, metabolites, genes, blocked reactions |
| `compartment_breakdown.png` | Pie chart of metabolite counts per compartment |
| `degree_distribution.png` | Histogram of metabolites per reaction (connectivity) |
| `subsystem_barchart.png` | Horizontal bar chart of top-25 subsystems by reaction count |
| `flux_distribution.png` | Bar chart of top-30 reactions by absolute FBA flux |
| `reaction_network.html` | Interactive Plotly bipartite reaction–metabolite network |
| `dashboard.html` | Single-page interactive dashboard (stats + subsystems + fluxes) |
| `escher_map.html` | Interactive Escher metabolic map with FBA flux overlay |
| `model_summary.csv` | Machine-readable table of key model statistics |

---

## Standalone visualization (existing model)

You can also visualize any existing SBML model without re-running the full pipeline:

```bash
python visualize_model.py --model output/model.xml --out-dir output/viz
```

Or from Python:

```python
import cobra
from pathlib import Path
from visualize_model import run_all_visualizations

model = cobra.io.read_sbml_model("output/model.xml")
run_all_visualizations(model, out_dir=Path("output/viz"))
```

---

## Downstream analysis (COBRApy)

```python
import cobra

model = cobra.io.read_sbml_model("output/model.xml")
print(model.optimize())              # FBA
print(model.summary())               # flux summary
```

---

## Troubleshooting

| Symptom | Fix |
|---------|-----|
| `carve: command not found` | Activate the conda env: `conda activate gsmm_pipeline` |
| `No CDS features found` | Ensure the `.gbk` file contains annotated CDS features |
| FBA status is `infeasible` | Re-run with `--gap-fill M9` (or another medium) |
| Diamond DB errors | Re-run `carve --download` to refresh the reference DB |
| Visualizations not generated | Ensure `--visualize` flag is passed and `visualize_model.py` is in the same directory |
| `plotly` / `escher` not found | Re-create the conda env: `conda env create -f environment.yml` |
| Escher map shows wrong organism | Supply your own Escher JSON map (see [escher docs](https://escher.readthedocs.io)) |
