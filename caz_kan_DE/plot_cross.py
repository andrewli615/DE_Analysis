import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

results = pd.read_csv('results_caz.csv', index_col=0)
results_kan = pd.read_csv('results_kan.csv', index_col=0)
primary = pd.read_csv('betalactam_primary.csv', index_col=0)
primary_kan = pd.read_csv('kanamycin_primary.csv', index_col=0)

# Cross-reactivity
caz_genes = set(primary['gene'].dropna())
kan_genes = set(primary_kan['gene'].dropna())
overlap = caz_genes.intersection(kan_genes)
print(f"Genes upregulated by both: {len(overlap)}")

key_candidates = ['recA', 'recN', 'emrA', 'pspD', 'lipA', 'dgcZ']
print("\nKey candidate cross-reactivity:")
for gene in key_candidates:
    print(f"{gene}: {'CROSS-REACTIVE' if gene in kan_genes else 'specific to ceftazidime'}")

kan_specific = primary_kan[~primary_kan['gene'].isin(caz_genes)].sort_values('log2FoldChange', ascending=False)
print(f"\nKanamycin specific candidates: {len(kan_specific)}")
kan_specific.to_csv('kanamycin_specific.csv')

def volcano_plot(results, title, output_path, highlight_genes=None):
    fig, ax = plt.subplots(figsize=(10, 8))
    df = results.dropna(subset=['log2FoldChange', 'padj'])
    df['neg_log10_padj'] = -np.log10(df['padj'])
    df['color'] = 'grey'
    df.loc[(df['log2FoldChange'] > 2) & (df['padj'] < 0.05), 'color'] = 'red'
    df.loc[(df['log2FoldChange'] < -2) & (df['padj'] < 0.05), 'color'] = 'blue'
    ax.scatter(df['log2FoldChange'], df['neg_log10_padj'], c=df['color'], alpha=0.5, s=10)
    ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.8)
    ax.axvline(2, color='black', linestyle='--', linewidth=0.8)
    ax.axvline(-2, color='black', linestyle='--', linewidth=0.8)
    if highlight_genes:
        for gene in highlight_genes:
            gene_data = df[df['gene'] == gene]
            if not gene_data.empty:
                ax.annotate(gene, (gene_data['log2FoldChange'].values[0],
                            gene_data['neg_log10_padj'].values[0]), fontsize=8, ha='right')
    ax.set_xlabel('log2 Fold Change')
    ax.set_ylabel('-log10 adjusted p-value')
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.show()

volcano_plot(results, 'Ceftazidime vs Control', 'volcano_ceftazidime.png',
    highlight_genes=['recA', 'recN', 'dgcZ', 'lipA', 'emrA', 'pspD'])

volcano_plot(results_kan, 'Kanamycin vs Control', 'volcano_kanamycin.png',
    highlight_genes=['ibpB', 'pdeL', 'rybB', 'ompG'])