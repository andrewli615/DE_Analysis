# Promoter Selection Differential Expression Analysis

## Goal
Identify E. coli promoters activated by aminoglycoside and beta-lactam
antibiotics for use as biosensor components.

## Directory Structure

```text
data/
  aminoglycoside/
    aminoglycoside_candidates.csv   # tracked fallback input for tobramycin
  caz_kan/
    GSE220559_RAW.tar               # tracked raw count archive

scripts/
  aminoglycoside/
    filter.py
    filter_tr_untr.py
    mapping.py
    compare.py
    summarize.py
  caz_kan/
    mapping.py
    deseq2_caz.py
    deseq2_kan.py
    summarize.py
    interactive_plot.py

outputs/
  aminoglycoside/                   # regenerated locally, git-ignored
    final/
    intermediate/
  caz_kan/                          # regenerated locally, git-ignored
    final/
      detailed_lists/
    intermediate/
    plots/
    scratch/
```

`outputs/` is intentionally ignored by git. Delete it any time you want to run
the analysis from scratch.

## Dependencies

```bash
python3 -m pip install pandas numpy pydeseq2 matplotlib plotly openpyxl
```

`matplotlib` and `plotly` are only needed for plots. Use `--no-plots` when you
only want CSV summaries.

## Run CAZ/KAN Analysis

From the repository root:

```bash
python3 scripts/caz_kan/mapping.py
python3 scripts/caz_kan/deseq2_caz.py
python3 scripts/caz_kan/deseq2_kan.py
python3 scripts/caz_kan/summarize.py --no-plots
```

Main readable output:

```text
outputs/caz_kan/final/promoter_summary.csv
```

Additional regenerated outputs include DESeq result CSVs, per-category promoter
lists, and optional volcano PNG/HTML plots:

```text
outputs/caz_kan/intermediate/       # DESeq CSVs and mapping files
outputs/caz_kan/final/              # one-file readable summary
outputs/caz_kan/final/detailed_lists/
outputs/caz_kan/plots/
outputs/caz_kan/scratch/            # extracted raw archive
```

The raw GSE220559 archive stays compressed in `data/caz_kan/GSE220559_RAW.tar`.
The scripts extract it into `outputs/caz_kan/scratch/GSE220559_RAW/` when needed. If
`data/caz_kan/ecoli_k12.gtf.gz` is absent, `mapping.py` writes an ID-only
mapping so the pipeline can still run; add the GTF file for named genes.

## Run Aminoglycoside Summary

From the repository root:

```bash
python3 scripts/aminoglycoside/summarize.py
```

Main readable output:

```text
outputs/aminoglycoside/final/promoter_summary.csv
```

The summary uses `data/aminoglycoside/aminoglycoside_candidates.csv` as a
fallback tobramycin input. To regenerate aminoglycoside intermediates from raw
data, place these optional inputs in `data/aminoglycoside/`:

```text
GSE224240_analysis.xlsx
GSE228373_RAW/
ecoli_annotation.gtf
```

Then run:

```bash
python3 scripts/aminoglycoside/filter.py
python3 scripts/aminoglycoside/filter_tr_untr.py
python3 scripts/aminoglycoside/mapping.py
python3 scripts/aminoglycoside/compare.py
python3 scripts/aminoglycoside/summarize.py
```

## Output Conventions

Promoter summaries classify rows as:

- `upregulated`: `log2FoldChange > 2` and, where available, `padj < 0.05`
- `downregulated`: `log2FoldChange < -2` and, where available, `padj < 0.05`
- `not_regulated`: all other rows

Readable summary CSVs round numeric values to two decimal places and sort by
signal strength.
