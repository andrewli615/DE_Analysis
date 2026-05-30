import pandas as pd
import numpy as np

# Load treated and untreated
treated = pd.read_csv('/Users/nichitabulbuc/Desktop/DE_analysis/GSE228373_RAW/GSM7119829_476822_S10.txt', 
                      sep='\t', names=['gene_id', 'treated'])
untreated = pd.read_csv('/Users/nichitabulbuc/Desktop/DE_analysis/GSE228373_RAW/GSM7119827_385274_DKP.txt', 
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

sig.to_csv('/Users/nichitabulbuc/Desktop/DE_analysis/gentamicin_candidates.csv', index=False)