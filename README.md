
# Promoter Selection — Differential Expression Analysis

## Goal
Identify E. coli promoters activated by aminoglycoside and beta-lactam antibiotics 
for use as biosensor components. Upregulated genes are candidates for cloning 
upstream of a GFP reporter by wet lab.

## Datasets
| Dataset | Antibiotic | Strain | Samples | Notes |
|---------|-----------|--------|---------|-------|
| GSE224240 | Tobramycin (sub-MIC) | E. coli K-12 MG1655 | 12 | Pre-computed DESeq2 results |
| GSE228373 | Gentamicin | E. coli JM109 | 2 | No replicates, fold change only, required GTF ID mapping |

## Scripts
- `filter.py` — filters GSE224240 for significantly upregulated genes (log2FC > 2, padj < 0.05)
- `filter_tr_untr.py` — treated vs untreated fold change comparison for GSE228373 (no replicates, no statistical testing)
- `mapping.py` — maps ER3413 gene IDs to standard gene names using GTF annotation
- `compare.py` — cross-references tobramycin and gentamicin candidate lists

## Results
- `aminoglycoside_candidates.csv` — tobramycin upregulated candidates
- `aminoglycoside_overlap.csv` — genes upregulated in both tobramycin and gentamicin (dppC, astD, puuE)

## Notes
- Overlapping candidates are metabolic genes(dppC, astD, puuE), not canonical stress response genes
- Known stress response promoters (marA, recA, cpxP) downregulated under tobramycin
