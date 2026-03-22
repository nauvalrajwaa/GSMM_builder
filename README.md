# GSMM Pipeline

Automated pipeline for generating a **Genome-Scale Metabolic Model (GSMM)**
from a FASTA + GenBank genomic input using **CarveMe** and **COBRApy**,
with optional **comprehensive visualizations**, built-in **NCBI data fetching**,
and a **Metabolic Cross-Feeding Network Builder** for metagenomics MAG inputs.

---

## Requirements

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda
- GLPK solver (bundled via conda) **or** CPLEX (optional, faster)
- Internet access for NCBI downloads and KEGG API calls

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

## Scripts overview

| Script | Purpose |
|--------|---------|
| `main.py` | **All-in-one runner**: fetch → model → visualize → crossfeed in one command |
| `fetch_ncbi.py` | Download FASTA + GenBank from NCBI RefSeq by accession |
| `generate_model.py` | Full GSMM pipeline: extract proteins → CarveMe → COBRApy QC |
| `visualize_model.py` | Comprehensive model visualization (standalone or via pipeline) |
| `gene_config.py` | Apply knockout / knockin gene manipulations and generate Plotly report |
| `crossfeed_network.py` | MAG cross-feeding network builder: KO annotations → bipartite metabolite flow graph |
| `run_tests.py` | End-to-end test runner: fetch → pipeline → visualize → crossfeed, with pass/fail report |

---

## Quick start: test with real data (E. coli K-12)

The fastest all-in-one way to run fetch + model + visualization + cross-feeding
using the built-in *E. coli* K-12 test dataset:

```bash
# One command – downloads genome, builds model, generates all visualizations
python main.py --email your@email.com
```

This will:
1. Download `NC_000913.3` (*E. coli* K-12 MG1655) from NCBI into `output/all_in_one/fetch/`
2. Run GSMM reconstruction for each accession
3. Generate full visualization outputs per accession
4. Run cross-feeding network analysis with giant-factory interactive view
5. Print a per-accession run summary

---

## One-command pipeline runner (`main.py`)

`main.py` orchestrates the complete workflow from repository root:

```bash
python main.py --email your@email.com
```

### Common `main.py` commands

```bash
# Custom accession list
python main.py --accessions NC_000913.3 NC_002695.2 --email your@email.com

# Set output root directory
python main.py --email your@email.com --run-dir output/my_run

# Run gap-filling during reconstruction
python main.py --email your@email.com --gap-fill M9

# Skip selected stages
python main.py --email your@email.com --skip-download
python main.py --email your@email.com --skip-crossfeed

# Continue processing remaining accessions after one fails
python main.py --accessions NC_000913.3 NC_002695.2 --email your@email.com --continue-on-error
```

### `main.py` options

```
  --accessions ACC [ACC ...]   RefSeq accession(s) to process (default: NC_000913.3)
  --email EMAIL                NCBI email for fetch step
  --run-dir DIR                Output root for the full run (default: output/all_in_one)
  --python PYTHON              Python interpreter used to invoke src scripts
  --gap-fill MEDIUM            Gap-fill medium for generate_model.py (e.g. M9)
  --max-api-kos N              Max KO terms resolved in crossfeed step (default: 220)
  --skip-download              Skip fetch step
  --skip-generate              Skip model generation step
  --skip-visualize             Skip visualization step
  --skip-crossfeed             Skip crossfeed step
  --continue-on-error          Continue remaining accessions after one fails
  --verbose / -v               Enable DEBUG logging
```

---

## Fetching NCBI data manually (`fetch_ncbi.py`)

```bash
# Fetch the built-in E. coli K-12 test genome
python src/fetch_ncbi.py --email your@email.com

# Fetch a specific accession
python src/fetch_ncbi.py --accessions NC_000913.3 --out-dir data/ecoli \
                     --email your@email.com

# Fetch multiple accessions
python src/fetch_ncbi.py --accessions NC_000913.3 NC_002695.2 \
                     --out-dir data/ --email your@email.com

# Skip download if files already exist (idempotent)
python src/fetch_ncbi.py --skip-existing --email your@email.com
```

### fetch_ncbi.py options

```
  --accessions ACC [ACC ...]   RefSeq accession(s) (default: NC_000913.3)
  --out-dir DIR                Output directory    (default: data/ncbi)
  --email EMAIL                NCBI email address  (required by NCBI policy)
  --skip-existing              Skip if FASTA+GBK already present
  --mag-table FILE             Output MAG annotation CSV path
  --skip-mag-table             Disable automatic MAG annotation generation
  --verbose / -v               Enable DEBUG logging
```

### Fetch outputs

| File | Description |
|------|-------------|
| `<outdir>/<accession>.fna` | Nucleotide FASTA |
| `<outdir>/<accession>.gbk` | Full GenBank annotation (with translations) |
| `<outdir>/mag_annotations.csv` | Auto-generated MAG annotation table for crossfeed (`mag_id,ko_id,gene_id,...`) |
| `<outdir>/ncbi_manifest.json` | Download manifest (accession, paths, status) |

---

## Running the full test suite (`run_tests.py`)

```bash
# Run all built-in tests
python src/run_tests.py --email your@email.com

# Re-run without re-downloading (data already in tests/)
python src/run_tests.py --skip-download

# Run with gap-filling on M9 medium
python src/run_tests.py --email your@email.com --gap-fill M9

# Test a custom accession
python src/run_tests.py --accessions NC_002695.2 --email your@email.com

# Keep intermediate files (proteins.faa) after run
python src/run_tests.py --keep-intermediates --email your@email.com
```

### run_tests.py options

```
  --accessions ACC [ACC ...]   Accession(s) to test (default: E. coli K-12 set)
  --base-dir DIR               Root dir for test work dirs  (default: tests)
  --email EMAIL                NCBI email for fetch_ncbi.py
  --skip-download              Skip NCBI fetch if files exist
  --gap-fill MEDIUM            CarveMe gap-fill medium (e.g. M9)
  --keep-intermediates         Keep proteins.faa after each test
  --report FILE                JSON report path  (default: tests/test_report.json)
  --verbose / -v               Enable DEBUG logging
```

### Test outputs (per accession)

```
tests/
└── NC_000913_3/
    ├── NC_000913.3.fna          ← downloaded FASTA
    ├── NC_000913.3.gbk          ← downloaded GenBank
    ├── ncbi_manifest.json       ← download manifest
    ├── pipeline.log             ← full pipeline log
    └── output/
        ├── model.xml            ← GSMM in SBML format
        ├── proteins.faa         ← intermediate (deleted unless --keep-intermediates)
        └── viz/
            ├── summary_stats.png
            ├── compartment_breakdown.png
            ├── degree_distribution.png
            ├── subsystem_barchart.png
            ├── flux_distribution.png
            ├── reaction_network.html
            ├── dashboard.html
            ├── escher_map.html
            ├── model_summary.csv
            └── crossfeed/              ← cross-feeding network (if mag_annotations.csv present)
                ├── index.html              ← MAG hub: 4-tab navigator
                ├── crossfeed_network.html  ← bipartite network (taxonomy colors + click-highlight)
                ├── crossfeed_sankey.html   ← metabolic flow Sankey diagram
                ├── crossfeed_heatmap.html  ← KO pathway completeness heatmap
                ├── keystone_report.html    ← Keystone species & dependency report
                ├── keystone_report.csv     ← Keystone CSV (dependency + provision scores)
                ├── crossfeed_edges.csv     ← edge list: producer → metabolite → consumer
                ├── crossfeed_summary.json  ← summary statistics
                └── kegg_cache.json         ← KEGG API response cache
tests/test_report.json           ← JSON pass/fail report
```

---

## Input files (manual use)

| File | Description |
|------|-------------|
| `genomics.fasta` | Genomic nucleotide or protein FASTA |
| `genomics.gbk` | Annotated GenBank file (CDS features required) |

Place both files in the project root,
or use `fetch_ncbi.py` to download them automatically.

---

## Usage: pipeline only

### Minimal (default paths, no gap-filling)
```bash
python src/generate_model.py
```

### Custom paths + gap-filling on M9 minimal medium
```bash
python src/generate_model.py \
    --fasta  genomics.fasta \
    --gbk    genomics.gbk \
    --output output/my_model.xml \
    --gap-fill M9
```

### With a log file
```bash
python src/generate_model.py --log-file pipeline.log
```

### With comprehensive visualizations
```bash
python src/generate_model.py --visualize
```

### Custom visualization output directory
```bash
python src/generate_model.py --visualize --viz-dir output/my_viz
```

### All options
```
usage: src/generate_model.py [-h] [--fasta FASTA] [--gbk GBK]
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
| `mutation_analysis.html` | Subsystem vulnerability KO plot + linked click-to-network highlighting |
| `escher_map.html` | Interactive Escher metabolic map with FBA flux overlay |
| `gene_config_report.html` | Gene config KO/KI Plotly report (if `gene_config.json` present) |
| `model_summary.csv` | Machine-readable table of key model statistics |
| `index.html` | Unified single-page hub: Overview, Interactive Reports, Environment Editor, Cross-Feeding Network |
| `crossfeed/crossfeed_network.html` | Bipartite cross-feeding network — taxonomy colors, click-to-highlight |
| `crossfeed/crossfeed_sankey.html` | Sankey metabolic flow diagram (optional abundance weighting) |
| `crossfeed/crossfeed_heatmap.html` | Pathway completeness heatmap: KO presence per MAG |
| `crossfeed/keystone_report.html` | Keystone Species & Dependency HTML report |
| `crossfeed/keystone_report.csv` | Keystone CSV: dependency score, provision score, keystone flag per MAG |
| `crossfeed/crossfeed_edges.csv` | Tabular edge list: producer MAG → metabolite → consumer MAG |
| `crossfeed/crossfeed_summary.json` | Machine-readable cross-feeding summary statistics |

---

## Standalone visualization (existing model)

```bash
python src/visualize_model.py --model output/model.xml --out-dir output/viz
```

**Auto-detected extras** — place these files in the project root alongside
the project root and they are picked up automatically:

| File | Effect |
|------|--------|
| `gene_config.json` | Runs KO/KI analysis and generates Escher mutation maps |
| `mag_annotations.csv` | Runs cross-feeding network and adds it to `index.html` (auto-generated by `fetch_ncbi.py` by default) |
| `sample_mag_annotations.csv` | Bundled 8-MAG gut-microbiome fixture for demos/manual experiments |
| `sample_kegg_cache.json` | Pre-built KEGG cache for offline/fast sample runs |

Or from Python:

```python
import cobra
from pathlib import Path
from visualize_model import run_all_visualizations

model = cobra.io.read_sbml_model("output/model.xml")
run_all_visualizations(model, out_dir=Path("output/viz"))
```

---

## MAG cross-feeding network (`crossfeed_network.py`)

Builds a **Metabolic Hand-off & Cross-Feeding Network** from metagenomics
Metagenome-Assembled Genomes (MAGs) annotated with KEGG Orthology (KO) terms.

For each MAG the tool resolves KO → KEGG reactions → substrate/product
compound sets, then finds *cross-feeding edges*: metabolites produced by one
MAG that are consumed by a different MAG. Results are exported as an
interactive Plotly bipartite graph, an edge-list CSV, and a JSON summary.

### Input format

A CSV or TSV file with at least two columns:

```
mag_id,ko_id,gene_id,description
MAG_001_Bifidobacterium,K00016,gene_003,L-lactate dehydrogenase
MAG_002_Firmicutes,K00016,gene_006,L-lactate dehydrogenase
MAG_002_Firmicutes,K00242,gene_010,Succinate dehydrogenase
...
```

Column names `mag_id` and `ko_id` are the defaults; override with
`--mag-col` / `--ko-col`.  Extra columns are ignored.

### Quickstart

```bash
# Generate the built-in sample fixture (8 gut-microbiome MAGs, 57 KO terms)
python src/crossfeed_network.py --generate-sample my_mags.csv

# Build cross-feeding network from the sample
python src/crossfeed_network.py --input my_mags.csv --out-dir output/crossfeed

# Use your own annotation table
python src/crossfeed_network.py --input my_real_mags.csv --out-dir output/crossfeed
```

### Common commands

```bash
# Minimal: just an input table (output goes to output/crossfeed)
python src/crossfeed_network.py --input mag_annotations.csv

# Custom output directory
python src/crossfeed_network.py \
    --input  mag_annotations.csv \
    --out-dir results/crossfeed

# TSV input with non-default column names
python src/crossfeed_network.py \
    --input   annotation_table.tsv \
    --mag-col bin_id \
    --ko-col  kegg_ko

# Reuse a previously downloaded KEGG cache (avoids repeated API calls)
python src/crossfeed_network.py \
    --input  mag_annotations.csv \
    --cache  output/crossfeed/kegg_cache.json

# Limit the number of unique KO terms sent to the KEGG API
# (useful for large datasets or rate-limited environments)
python src/crossfeed_network.py \
    --input       mag_annotations.csv \
    --max-api-kos 200

# Fully offline: set max-api-kos to 0 and supply a pre-built cache
python src/crossfeed_network.py \
    --input       mag_annotations.csv \
    --cache       my_kegg_cache.json \
    --max-api-kos 0

# Verbose logging (shows every KEGG API call and graph detail)
python src/crossfeed_network.py \
    --input mag_annotations.csv \
    --verbose

# Generate the sample fixture and immediately analyse it
python src/crossfeed_network.py --generate-sample /tmp/sample.csv
python src/crossfeed_network.py --input /tmp/sample.csv --out-dir /tmp/crossfeed_out
```

### Integrate with an existing GSMM visualization run

Place `mag_annotations.csv` in the project root.
It is auto-discovered and the cross-feeding network is appended as a new tab in
`index.html`:

```bash
# GSMM pipeline + visualizations + cross-feeding network (all in one)
python src/generate_model.py \
    --fasta  genome.fna \
    --gbk    genome.gbk \
    --output output/model.xml \
    --visualize

# The file mag_annotations.csv must exist in the project root.
# Results appear in output/viz/crossfeed/ and as a tab in output/viz/index.html
```

Alternatively, trigger only the cross-feeding step from within Python:

```python
from pathlib import Path
from crossfeed_network import run_crossfeed_analysis

results = run_crossfeed_analysis(
    table_path = Path("mag_annotations.csv"),
    out_dir    = Path("output/crossfeed"),
)
print(f"Found {len(results['edges'])} cross-feeding edges")
for edge in results["edges"][:5]:
    print(f"  {edge['producer_mag']} → {edge['compound_name']} → {edge['consumer_mag']}")
```

### crossfeed_network.py options

```
  --input FILE, -i FILE     MAG annotation CSV/TSV (required)
  --out-dir DIR, -o DIR     Output directory           (default: output/crossfeed)
  --cache FILE              Path to KEGG JSON cache    (default: <out-dir>/kegg_cache.json)
  --mag-col COL             MAG ID column name         (default: mag_id)
  --ko-col COL              KO ID column name          (default: ko_id)
  --max-api-kos N           Max unique KOs to resolve  (default: 500)
  --generate-sample FILE    Write sample CSV to FILE and exit
  --verbose / -v            Enable DEBUG logging
```

### Cross-feeding outputs

All written to `--out-dir` (default `output/crossfeed/`):

| File | Description |
|------|-------------|
| `index.html` | MAG hub: 4-tab navigator (Network, Sankey, Heatmap, Report) with stat cards |
| `crossfeed_network.html` | Bipartite graph: MAGs colored by taxonomy, click a metabolite to highlight its producers/consumers |
| `crossfeed_sankey.html` | Sankey metabolic flow: metabolite cascade from producers through metabolites to consumers |
| `crossfeed_heatmap.html` | Pathway completeness heatmap: KO presence by MAG with completeness % per column |
| `keystone_report.html` | Keystone Species & Dependency report: dependency score, provision score, keystone alerts |
| `keystone_report.csv` | CSV: `mag_id, dependency_score, provision_score, n_metabolites_provided, n_metabolites_received, is_keystone, keystone_metabolites` |
| `crossfeed_edges.csv` | Edge list: `producer_mag, compound_id, compound_name, consumer_mag` |
| `crossfeed_summary.json` | Summary stats: n_mags, n_edges, n_metabolites, top metabolites, per-MAG capacity |
| `kegg_cache.json` | Local cache of KEGG API results (reused on subsequent runs) |

### How cross-feeding edges are detected

1. For each KO term in each MAG, the KEGG REST API (`rest.kegg.jp`) is queried
   for the associated reactions and their substrate / product compound IDs.
2. A compound is flagged as a **cross-feeding metabolite** when at least one
   MAG *produces* it (it appears as a product of one of that MAG's reactions)
   and at least one *different* MAG *consumes* it (it appears as a substrate).
3. Each (producer MAG, compound, consumer MAG) triple becomes one edge in the
   bipartite graph.

Example output for gut microbiome MAGs:

```
MAG_001_Bifidobacterium  →  Pyruvate    →  MAG_002_Firmicutes
MAG_001_Bifidobacterium  →  L-Lactate   →  MAG_002_Firmicutes
MAG_004_Bacteroides      →  Pyruvate    →  MAG_002_Firmicutes
MAG_002_Firmicutes       →  NADH        →  MAG_001_Bifidobacterium
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
| NCBI fetch times out / fails | NCBI rate-limits unauthenticated requests; set `--email` and retry |
| `HTTP Error 429` during fetch | Wait a few minutes then re-run with `--skip-existing` |
| Cross-feeding: 0 edges found | KO terms not resolving to reactions — check internet access or supply a `--cache` file |
| Cross-feeding: KEGG rate limit | Reduce `--max-api-kos` or run again (cache will resume from where it stopped) |
| Cross-feeding: `requests` missing | `pip install requests` inside the conda env |
| `mag_annotations.csv` not picked up by `visualize_model.py` | Place the file in the same directory as `visualize_model.py` (project root) |
