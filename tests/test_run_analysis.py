import numpy as np
import pandas as pd

from scripts.run_analysis import (
    benjamini_hochberg,
    classify,
    finalize_results,
    median_ratio_size_factors,
)


def test_benjamini_hochberg_preserves_index_and_monotonic_adjustment():
    pvalues = pd.Series([0.04, 0.01, 0.03], index=["gene_c", "gene_a", "gene_b"])
    adjusted = benjamini_hochberg(pvalues)

    assert adjusted.index.tolist() == pvalues.index.tolist()
    assert np.allclose(adjusted.to_numpy(), [0.04, 0.03, 0.04])


def test_classify_uses_strict_configured_thresholds():
    assert classify(2.01, 0.049, 2.0, 0.05) == "upregulated"
    assert classify(-2.01, 0.049, 2.0, 0.05) == "downregulated"
    assert classify(2.0, 0.01, 2.0, 0.05) == "not_regulated"
    assert classify(3.0, 0.05, 2.0, 0.05) == "not_regulated"


def test_finalize_results_computes_adjustment_signal_and_order():
    raw = pd.DataFrame(
        {
            "gene": ["geneA", "geneB"],
            "gene_id": ["b0001", "b0002"],
            "log2FoldChange": [3.0, -1.0],
            "pvalue": [0.01, 0.5],
        }
    )
    dataset = {
        "name": "example",
        "antibiotic_class": "beta_lactam",
        "treatment": "example-drug",
        "comparison": "treated_vs_control",
    }

    result = finalize_results(
        raw,
        dataset,
        {"log2_fold_change": 2.0, "padj": 0.05},
    )

    assert result["gene"].tolist() == ["geneA", "geneB"]
    assert result["regulation"].tolist() == ["upregulated", "not_regulated"]
    assert result["padj_source"].tolist() == ["BH_from_pvalue", "BH_from_pvalue"]
    assert result["comparison"].tolist() == ["treated_vs_control"] * 2
    assert result.iloc[0]["signal_strength"] > result.iloc[1]["signal_strength"]


def test_median_ratio_size_factors_are_scale_invariant():
    counts = pd.DataFrame(
        {
            "sample_a": [10, 20, 30],
            "sample_b": [20, 40, 60],
        }
    )
    factors = median_ratio_size_factors(counts)
    normalized = counts.div(factors, axis=1)

    assert np.allclose(normalized["sample_a"], normalized["sample_b"])
