import pandas as pd

# Load both candidate lists
tobramycin = pd.read_csv('/Users/nichitabulbuc/Desktop/DE_analysis/aminoglycoside_candidates.csv')
gentamicin = pd.read_csv('/Users/nichitabulbuc/Desktop/DE_analysis/gentamicin_candidates_named.csv')

# Drop unnamed genes from gentamicin
gentamicin = gentamicin.dropna(subset=['gene'])

# Find overlap
overlap = tobramycin[tobramycin['gene'].isin(gentamicin['gene'])]

print(f"Overlapping genes: {len(overlap)}")
print(overlap[['gene', 'log2FoldChange', 'padj']].sort_values('log2FoldChange', ascending=False).to_string())

# Save
overlap.to_csv('/Users/nichitabulbuc/Desktop/DE_analysis/aminoglycoside_overlap.csv', index=False)