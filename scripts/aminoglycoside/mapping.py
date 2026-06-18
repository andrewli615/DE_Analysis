import pandas as pd
import re
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DATA_DIR = REPO_ROOT / 'data' / 'aminoglycoside'
OUTPUT_DIR = REPO_ROOT / 'outputs' / 'aminoglycoside'
INTERMEDIATE_DIR = OUTPUT_DIR / 'intermediate'
INTERMEDIATE_DIR.mkdir(parents=True, exist_ok=True)

# Parse GTF to extract ER3413 ID to gene name mapping
mapping = {}
with open(DATA_DIR / 'ecoli_annotation.gtf', 'r') as f:
    for line in f:
        if 'old_locus_tag' in line and 'gene_id' in line:
            old = re.search(r'old_locus_tag "([^"]+)"', line)
            gene = re.search(r'gene "([^"]+)"', line)
            if old and gene:
                # Convert ER3413_1 format to ER3413_1
                old_id = old.group(1)
                mapping[old_id] = gene.group(1)

print(f"Mapped {len(mapping)} genes")

# Load candidates
candidates = pd.read_csv(INTERMEDIATE_DIR / 'gentamicin_candidates.csv')

# Map gene names
candidates['gene'] = candidates['gene_id'].map(mapping)

print(candidates[['gene_id', 'gene', 'log2FC']].head(20))
candidates.to_csv(INTERMEDIATE_DIR / 'gentamicin_candidates_named.csv', index=False)
