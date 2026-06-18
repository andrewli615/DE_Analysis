import pandas as pd
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DATA_DIR = REPO_ROOT / 'data' / 'aminoglycoside'
OUTPUT_DIR = REPO_ROOT / 'outputs' / 'aminoglycoside'
INTERMEDIATE_DIR = OUTPUT_DIR / 'intermediate'
INTERMEDIATE_DIR.mkdir(parents=True, exist_ok=True)

# Load both candidate lists
tobramycin_path = INTERMEDIATE_DIR / 'aminoglycoside_candidates.csv'
if not tobramycin_path.exists():
    tobramycin_path = DATA_DIR / 'aminoglycoside_candidates.csv'
tobramycin = pd.read_csv(tobramycin_path)
gentamicin = pd.read_csv(INTERMEDIATE_DIR / 'gentamicin_candidates_named.csv')

# Drop unnamed genes from gentamicin
gentamicin = gentamicin.dropna(subset=['gene'])

# Find overlap
overlap = tobramycin[tobramycin['gene'].isin(gentamicin['gene'])]

print(f"Overlapping genes: {len(overlap)}")
print(overlap[['gene', 'log2FoldChange', 'padj']].sort_values('log2FoldChange', ascending=False).to_string())

# Save
overlap.to_csv(INTERMEDIATE_DIR / 'aminoglycoside_overlap.csv', index=False)
