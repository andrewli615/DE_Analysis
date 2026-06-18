import plotly.express as px
import pandas as pd
import numpy as np
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
OUTPUT_DIR = REPO_ROOT / 'outputs' / 'caz_kan'
INTERMEDIATE_DIR = OUTPUT_DIR / 'intermediate'
PLOT_DIR = OUTPUT_DIR / 'plots'
PLOT_DIR.mkdir(parents=True, exist_ok=True)

results = pd.read_csv(INTERMEDIATE_DIR / 'results_caz.csv', index_col=0)
results_kan = pd.read_csv(INTERMEDIATE_DIR / 'results_kan.csv', index_col=0)

def interactive_volcano(results, title, output_path):
    df = results.dropna(subset=['log2FoldChange', 'padj'])
    df['neg_log10_padj'] = -np.log10(df['padj'])
    df['significant'] = 'not significant'
    df.loc[(df['log2FoldChange'] > 2) & (df['padj'] < 0.05), 'significant'] = 'upregulated'
    df.loc[(df['log2FoldChange'] < -2) & (df['padj'] < 0.05), 'significant'] = 'downregulated'

    fig = px.scatter(df, x='log2FoldChange', y='neg_log10_padj',
                     color='significant',
                     color_discrete_map={
                         'upregulated': 'red',
                         'downregulated': 'blue',
                         'not significant': 'grey'
                     },
                     hover_name='gene',
                     hover_data=['log2FoldChange', 'padj'],
                     title=title)
    fig.add_hline(y=-np.log10(0.05), line_dash='dash', line_color='black')
    fig.add_vline(x=2, line_dash='dash', line_color='black')
    fig.add_vline(x=-2, line_dash='dash', line_color='black')
    fig.update_traces(marker=dict(size=5, opacity=0.6))
    fig.update_layout(xaxis_title='log2 Fold Change', yaxis_title='-log10 adjusted p-value')
    fig.write_html(output_path)
    print(f"Saved to {output_path}")

interactive_volcano(results, 'Ceftazidime vs Control',
                    PLOT_DIR / 'volcano_ceftazidime_interactive.html')
interactive_volcano(results_kan, 'Kanamycin vs Control',
                    PLOT_DIR / 'volcano_kanamycin_interactive.html')
