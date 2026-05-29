import pandas as pd

# Load data
df = pd.read_excel('GSE224240_analysis.xlsx')

# Filter for significantly upregulated genes
sig = df[(df['log2FoldChange'] > 2) & (df['padj'] < 0.05)]

# Sort by fold change
sig = sig.sort_values('log2FoldChange', ascending=False)

# Preview
print(sig[['index', 'gene', 'log2FoldChange', 'padj']])

# Save
sig.to_csv('aminoglycoside_candidates.csv', index=False)