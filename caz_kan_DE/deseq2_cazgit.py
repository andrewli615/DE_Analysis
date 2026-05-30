import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import glob

path = 'GSE220559_RAW/'
mapping = pd.read_csv('gene_mapping.csv', index_col=0)['gene']

caz_files = sorted(glob.glob(path + '*CAZ*.txt'))
wu_files = sorted(glob.glob(path + '*Wu*.txt'))

def load_counts(files):
    dfs = []
    for f in files:
        df = pd.read_csv(f, sep='\s+', index_col=0)
        df = df.iloc[:, 0:1]
        df.columns = [f.split('/')[-1].split('.')[0]]
        dfs.append(df)
    return pd.concat(dfs, axis=1)

caz = load_counts(caz_files)
wu = load_counts(wu_files)
counts = pd.concat([caz, wu], axis=1).T

metadata = pd.DataFrame({
    'condition': ['ceftazidime']*3 + ['control']*3
}, index=counts.index)

dds = DeseqDataSet(counts=counts, metadata=metadata, design_factors='condition')
dds.deseq2()
stats = DeseqStats(dds, contrast=['condition', 'ceftazidime', 'control'])
stats.summary()
results = stats.results_df
results['gene'] = results.index.map(mapping)

primary = results[(results['log2FoldChange'] > 2) & (results['padj'] < 0.05)].sort_values('log2FoldChange', ascending=False)
secondary = results[(results['log2FoldChange'] > 1) & (results['log2FoldChange'] <= 2) & (results['padj'] < 0.05)].sort_values('log2FoldChange', ascending=False)

print(f"Primary candidates: {len(primary)}")
print(primary[['gene', 'log2FoldChange', 'padj']].head(20).to_string())
print(f"\nSecondary candidates: {len(secondary)}")
print(secondary[['gene', 'log2FoldChange', 'padj']].head(20).to_string())

results.to_csv('results_caz.csv')
primary.to_csv('betalactam_primary.csv')
secondary.to_csv('betalactam_secondary.csv')