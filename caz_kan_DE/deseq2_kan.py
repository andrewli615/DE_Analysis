import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import glob

path = 'GSE220559_RAW/'
mapping = pd.read_csv('gene_mapping.csv', index_col=0)['gene']

kan_files = sorted(glob.glob(path + '*KAN*.txt'))
wu_files = sorted(glob.glob(path + '*Wu*.txt'))

def load_counts(files):
    dfs = []
    for f in files:
        df = pd.read_csv(f, sep='\s+', index_col=0)
        df = df.iloc[:, 0:1]
        df.columns = [f.split('/')[-1].split('.')[0]]
        dfs.append(df)
    return pd.concat(dfs, axis=1)

kan = load_counts(kan_files)
wu = load_counts(wu_files)
counts = pd.concat([kan, wu], axis=1).T

metadata = pd.DataFrame({
    'condition': ['kanamycin']*3 + ['control']*3
}, index=counts.index)

dds = DeseqDataSet(counts=counts, metadata=metadata, design_factors='condition')
dds.deseq2()
stats = DeseqStats(dds, contrast=['condition', 'kanamycin', 'control'])
stats.summary()
results_kan = stats.results_df
results_kan['gene'] = results_kan.index.map(mapping)

primary_kan = results_kan[(results_kan['log2FoldChange'] > 2) & (results_kan['padj'] < 0.05)].sort_values('log2FoldChange', ascending=False)

print(f"Kanamycin primary candidates: {len(primary_kan)}")
print(primary_kan[['gene', 'log2FoldChange', 'padj']].head(20).to_string())

results_kan.to_csv('results_kan.csv')
primary_kan.to_csv('kanamycin_primary.csv')