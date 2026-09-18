# Regulatory network and promoter review

The primary analysis uses one consistent RNA-seq study, GSE220559, and a
same-release curated E. coli K-12 reference from RegulonDB. It produces class
proxy gene shortlists and promoter annotations for experimental review.
It does not establish response to the four wet-lab drugs or predict reporter
performance.

## Primary benchmark

`config/benchmark.json` defines:

| Comparison against water control | Role |
| --- | --- |
| Ceftazidime | Beta-lactam discovery proxy |
| Imipenem | Within-study beta-lactam robustness comparison |
| Kanamycin | Aminoglycoside discovery proxy |
| Ciprofloxacin | Cross-reactivity challenge |
| Polymyxin E | Cross-reactivity challenge |

Wet-lab targets are **amoxicillin, cephalexin, gentamicin, and tobramycin**.
None is directly measured in this benchmark. Two beta-lactam comparisons
provide within-study support, not independent validation; the aminoglycoside
arm has only one representative. Shared controls make comparisons dependent.
The original experiment uses IC50 exposures and growth to mid-log phase;
IC50 is not half the MIC, and this is not a fixed-duration acute exposure.

GSE220559 is the only antibiotic-expression study retained in the repository.
The older heterogeneous study inputs, per-dataset exports, and supplementary
configuration have been removed. Benchmark exports are grouped by treatment in
`results/differential_expression/contrasts/`; stable dataset IDs remain in the evidence tables.

## Run the workflow

From the repository root, install the project requirements using
Python >=3.11 in a virtual environment:

```bash
python -m pip install -r requirements.txt
python -m promoter_discovery.run_benchmark
python -m promoter_discovery.setup_data --manifest config/regulondb_assets.json --lock config/regulondb.lock.json --download
python -m promoter_discovery.build_network
python -m promoter_discovery.score_candidates
python -m promoter_discovery.promoter_selection
python -m promoter_discovery.visualize_network
```

The benchmark reads the existing `data/raw/GSE220559/GSE220559_RAW.tar` archive.
The primary reference manifest requests the network, sigma interactions,
gene mapping, TUSet, and PromoterSet from the production RegulonDB endpoint.
It does not require PRECISE. Hashes, remote product identifiers, and observed
release are recorded in `config/regulondb.lock.json`.
Network construction and promoter mapping reject unverified or mismatched
reference files. HTTPS verification stays enabled; configure a trusted CA
bundle with `SSL_CERT_FILE` if the Python installation lacks the required CA.

Run `python -m promoter_discovery.setup_data --manifest config/regulondb_assets.json --lock config/regulondb.lock.json`
to verify cached primary assets. Do not refresh the lock merely to bypass a
checksum mismatch. An intentional reference update requires reviewing the
release and rebuilding downstream outputs.

## Count model and response evidence

`run_benchmark.py` validates finite nonnegative integer raw counts, combines
identical shared samples once, and fits one PyDESeq2 negative-binomial model
with `~condition`. Every selected condition uses the same joint normalization
and dispersion estimation. Genes require >=10 counts in >=3 of the 18 unique
samples, regardless of treatment, before fitting.

Matched drug-versus-water contrasts and all ten pairwise drug contrasts use
Wald tests and BH adjustment across genes separately per contrast. Independent
filtering is disabled after the common prefilter; Cook's filtering remains
enabled. Unavailable tests retain missing p-values and are excluded from the
network background. Approximate 95% Wald intervals accompany fold changes.
The effect cutoff remains the configured absolute log2FC >2 and adjusted
p-value <0.05; significance tests the zero-effect null, not the cutoff itself.
Pairwise p-values are adjusted per contrast, not across all ten comparisons;
exploratory screening of many pairs warrants that additional caution.

`results/differential_expression/model/` contains joint normalized counts, unique sample
metadata, the common filter, direct drug contrasts, and run provenance. Config, DE-export, and direct-contrast hashes prevent
mixing outputs from different joint fits.
`results/differential_expression/contrasts/<treatment>/de_results.csv` preserves full precision,
standard errors, intervals, eligibility, and normalized basal/induced means.
Normalized counts describe relative transcript abundance, not absolute
promoter strength or fluorescence. The optional custom runner (`python -m promoter_discovery.input_data --config <config>`)
rejects this config because it would otherwise use the legacy Welch method.

Network files, enrichment tables, scoring provenance, and the interactive HTML
are written to `results/regulatory_network/`. Class shortlists and full candidate
annotations are in `results/promoter_candidates/`. Promoter mappings and review
entries are in `results/promoter_review/`. See the [results guide](../results/README.md).

## Full network versus displayed candidates

The serialized graph retains all accepted curated regulatory edges and targets,
even without a selected candidate. Candidate selection affects the displayed
subgraph and shortlists, never the statistical universe. Challenge-only hits
do not seed candidates. Use `--candidate-direction either-direction` to
include repression candidates, and `--top-n` for an optional display/shortlist
filter. The HTML defaults to candidates and their regulatory predecessors;
`--full-reference` displays the entire stored network.

Regulatory edge types are `activates`, `represses`, and `dual`. The builder
preserves underlying records when a regulator-target pair has multiple
annotations; opposing effects become `dual` at gene level. This does not
resolve individual promoter or condition contexts. `co-imodulon` relationships
are separate and never count as curated regulatory edges.

## Enrichment and candidate ranking

`regulon_enrichment.csv` tests each full regulon's overlap with upregulated and
downregulated genes separately for every configured contrast. The background
is all measured, eligible genes with verified same-release identity mapping,
including genes without known regulatory edges. Targets outside that
background are excluded from the contingency table but remain in the graph.
Each canonical gene counts once; opposite significant probe estimates cannot
enter either directional response set.

One-sided Fisher exact tests ask whether overlap exceeds chance. BH adjustment
covers the whole family of regulators x contrasts x two directions, including
zero-overlap tests for every regulon with a measured target. Reported effect
sign counts distinguish patterns consistent with activation, repression, and
ambiguous dual regulation. They are descriptive hypotheses, not measurements
of TF activity or proof that the antibiotic acts directly on that regulator.
Small regulons, incomplete annotation, correlated operon genes, and broad
stress responses limit interpretation. Follow-up should assess whether an
apparent enrichment is driven by a single transcription unit.

`regulon_candidate_coverage.csv` describes selected candidates as a fraction
of the full accepted regulon. It is **not antibiotic specificity**. The old
`tf_specificity_scores.csv`, `regulatory_clarity`, and calibrated-sounding
recommendation thresholds are removed from generated outputs.

`annotated_candidates.csv` includes upstream regulators, full response
profiles, approximate uncertainty, regulon evidence, direct drug comparisons,
and explicit validation requirements. A non-significant challenge is labelled
`no_detected_response`, never proof of no response or of class specificity.
Profiles include significant repression even under upregulated-only seeding.

Support labels are:

- `single_comparison_supported`: one qualifying comparison in a class.
- `multiple_comparisons_supported`: more than one qualifying comparison in a
  class, without significant opposite directions within that class.
- `within_class_direction_difference`: significant opposite directions within
  a class; this is visible for interpretation, not labelled technical failure.

Opposite responses across classes are retained as
`cross_class_direction_difference`, not a global conflict. No label asserts
independent validation. Candidates sort by support tier, then their largest
absolute qualifying discovery/robustness fold change, then gene name.
Challenge effects and external annotations cannot improve that ordering. Each
class shortlist uses its own class support tier and class effect size; support
from the other class cannot improve its position.
`beta_lactam_candidates.csv` and `aminoglycoside_candidates.csv` provide
separate class shortlists, with the full response evidence retained.

## Promoter mapping and construct review

`promoter_selection.py` joins canonical candidate names to TUSet, then joins
TU promoter IDs to PromoterSet. `candidate_promoter_mapping.csv` retains
multiple transcription units/promoters and explicitly unmapped genes.
`promoter_review.csv` groups genes sharing a promoter into a single review
entry and includes TSS, strand, sigma factor, annotated sequence, evidence,
confidence, and literature identifiers.

The annotated sequence is **not a ready-to-order fragment**. Operator positions,
all regulatory sites, genomic context, construct boundaries, and the intended
reporter host still need review. This stage does not extract arbitrary upstream
DNA or assume every responsive operon gene has its own promoter. Reporter
measurements with the actual four antibiotics establish transfer and practical
cross-reactivity.

## Optional PRECISE-1K context

The separate `config/precise1k_assets.json` manifest additionally provides a pinned
PRECISE-1K model and same-release companion expression matrix. After acquiring
those optional assets:

```bash
python -m promoter_discovery.setup_data --manifest config/precise1k_assets.json --lock config/precise1k.lock.json --download
python -m promoter_discovery.i_modulon_analysis --precise data/references/precise1k/model.json.gz --expression data/references/precise1k/expression.csv --mapping data/references/regulondb/network_gene_mapping.tsv
python -m promoter_discovery.score_candidates --graph results/regulatory_network/precise1k/regulatory_network_imodulon.pkl --out results/precise1k_context
python -m promoter_discovery.visualize_network --graph results/regulatory_network/precise1k/regulatory_network_imodulon.pkl --out results/regulatory_network/precise1k/regulatory_network.html
```

The default scorer uses the primary graph. Explicitly select the annotated
graph to include external context. iModulon memberships, PRECISE expression,
and activity remain secondary annotations; unmatched conditions do not
validate benchmark responses. Metabolism/translation membership weights are
heuristics, not experimentally measured host or reporter burden, and never
change primary ranking. Do not infer a new ICA regulatory network from this
small antibiotic panel alone.

## Sources and verification

- [GSE220559 original study](https://journals.asm.org/doi/10.1128/spectrum.00317-23)
- [RegulonDB database paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10767902/)
- [DESeq2 statistical framework](https://doi.org/10.1186/s13059-014-0550-8)
- [PyDESeq2 API](https://pydeseq2.readthedocs.io/en/stable/api/docstrings/pydeseq2.dds.DeseqDataSet.html)
- [PRECISE-1K reference](https://doi.org/10.1093/nar/gkad750)

Run `python -m pytest -q` for integration and scientific-contract regressions.
Tests cover full-regulon denominators, fixed enrichment backgrounds under
candidate filtering, global correction, sign preservation, challenge seeding,
shared controls, raw-count validation, optional external context, and operon
promoter grouping.
