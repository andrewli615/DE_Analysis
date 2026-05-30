# Promoter Selection — Differential Expression Analysis

## Goal
Identify E. coli promoters activated by aminoglycoside and beta-lactam antibiotics 
for use as biosensor components. Upregulated genes are candidates for cloning 
upstream of a GFP reporter by wet lab.

## Repository Structure
DE_Analysis/
├── aminoglycoside/
│   ├── filter.py
│   ├── filter_tr_untr.py
│   ├── mapping.py
│   └── compare.py
└── caz_kan_DE/
├── mapping.py
├── deseq2_caz.py
├── deseq2_kan.py
└── plots_and_crossreactivity.py

## Datasets
| Dataset | Antibiotic | Strain | Samples | Notes |
|---------|-----------|--------|---------|-------|
| GSE224240 | Tobramycin (sub-MIC) | E. coli K-12 MG1655 | 12 | Pre-computed DESeq2 results |
| GSE228373 | Gentamicin | E. coli JM109 | 2 | No replicates, fold change only |
| GSE220559 | Ceftazidime + Kanamycin | E. coli K-12 MG1655 | 33 | 3 replicates, full DESeq2 pipeline |

## Analysis 1 — Aminoglycoside (aminoglycoside/)
Analyzed GSE224240 (tobramycin) and GSE228373 (gentamicin). Cross-referencing 
gave 3 overlapping candidates (dppC, astD, puuE) but all were metabolic genes. 
Known stress response promoters (marA, recA, cpxP) were downregulated — consistent 
with known aminoglycoside ribosomal stress mechanism.

### Scripts
- `filter.py` — filters GSE224240 for upregulated genes (log2FC > 2, padj < 0.05)
- `filter_tr_untr.py` — fold change comparison for GSE228373 (no replicates)
- `mapping.py` — maps ER3413 gene IDs to standard names using GTF annotation
- `compare.py` — cross-references tobramycin and gentamicin candidate lists

## Analysis 2 — Beta-Lactam & Aminoglycoside (caz_kan_DE/)
Analyzed GSE220559 using full DESeq2 pipeline with 3 replicates per condition.
Ceftazidime (beta-lactam) and kanamycin (aminoglycoside) vs water control.

### Scripts
Run in this order:
1. `mapping.py` — extracts b-number to gene name mapping from E. coli K-12 GTF
2. `deseq2_caz.py` — DESeq2 analysis for ceftazidime vs control
3. `deseq2_kan.py` — DESeq2 analysis for kanamycin vs control
4. `plots_and_crossreactivity.py` — volcano plots and cross-reactivity check

### Dependencies
pip install pydeseq2 pandas matplotlib numpy

## Key Findings
**Beta-lactam (ceftazidime)**
- 112 primary candidates (log2FC > 2), strong SOS response signal
- Best candidate: recA — specific to beta-lactam, well characterized, commonly used in iGEM biosensors
- dgcZ, lipA — specific to ceftazidime, strong signal

**Aminoglycoside (kanamycin)**
- 896 primary candidates, broad transcriptional response
- Best candidate: ibpB — specific to kanamycin, responds to aminoglycoside ribosomal stress
- 60 genes cross-reactive between ceftazidime and kanamycin — avoided as promoter candidates

## Notes
- Datasets use K-12 MG1655
- Raw data: GSE224240, GSE228373, GSE220559 available on GEO
- GTF annotation: E. coli K-12 MG1655 GCF_000005845.2 from NCBI
