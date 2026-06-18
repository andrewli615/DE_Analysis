import pandas as pd
import numpy as np
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
RAW_DIR = REPO_ROOT / 'data' / 'aminoglycoside' / 'GSE228373_RAW'
OUTPUT_DIR = REPO_ROOT / 'outputs' / 'aminoglycoside'
INTERMEDIATE_DIR = OUTPUT_DIR / 'intermediate'
INTERMEDIATE_DIR.mkdir(parents=True, exist_ok=True)

# Load treated and untreated
treated = pd.read_csv(RAW_DIR / 'GSM7119829_476822_S10.txt', 
                      sep='\t', names=['gene_id', 'treated'])
untreated = pd.read_csv(RAW_DIR / 'GSM7119827_385274_DKP.txt', 
                      sep='\t', names=['gene_id', 'untreated'])

# Convert to numeric
treated['treated'] = pd.to_numeric(treated['treated'], errors='coerce')
untreated['untreated'] = pd.to_numeric(untreated['untreated'], errors='coerce')

# Merge
df = treated.merge(untreated, on='gene_id')

# Calculate log2 fold change
df['log2FC'] = np.log2((df['treated'] + 1) / (df['untreated'] + 1))

# Filter upregulated
sig = df[df['log2FC'] > 2].sort_values('log2FC', ascending=False)
print(f"Candidates: {len(sig)}")
print(sig.head(20))

sig.to_csv(INTERMEDIATE_DIR / 'gentamicin_candidates.csv', index=False)
