# From Antibiotic Exposure to Detectable Fluorescence

## Page purpose

Our dry-lab work asks one connected design question:

> Which *E. coli* promoters can convert clinically relevant antibiotic exposure
> or resistance-associated physiology into a rapid, specific, and measurable
> fluorescent signal in our device?

We address this question at three scales. Transport models estimate the
antibiotic concentration experienced by the cells. Transcriptomic and
regulatory-network analyses identify candidate promoters. A fluorescence model
connects promoter activity to the signal measured by the detector.

This page should emphasize how each model changed an experimental decision. It
should not read as a collection of independent scripts.

## Suggested opening figure

Redraw the following as one horizontal project diagram:

```text
Patient sample
     |
     v
Hydrogel diffusion and degradation
     |
     v
In-chip transport and cell exposure
     |
     v
Antibiotic-responsive regulatory network
     |
     v
Candidate promoter -> GFP transcription -> GFP maturation
     |
     v
Measured fluorescence and detection decision
```

Under the figure, state the contribution of each model in one sentence:

- Hydrogel model: predicts how much antibiotic reaches the cells and when.
- In-chip model: predicts spatial exposure, mixing, and residence time.
- Differential expression: finds genes associated with antibiotic response or
  resistance.
- Regulatory networks: connects responsive genes to actual promoters,
  transcription factors, sigma factors, and operons.
- Lab promoter library: constrains the computational shortlist to constructs
  that can be tested.
- Fluorescence model: predicts response time, dynamic range, and limit of
  detection.

## 1. Design objective

The device is intended to identify resistant *E. coli* using fluorescence.
Promoter selection therefore cannot depend only on a large fold change. A useful
promoter should combine:

- a sufficiently low baseline signal;
- a reproducible induced signal;
- a response that occurs within the device's assay time;
- specificity for an antibiotic, antibiotic class, or resistance phenotype;
- activity across relevant *E. coli* strains and growth conditions;
- a promoter region that can be identified and cloned;
- compatibility with the promoter parts available to the team.

The computational work reduces a genome-scale search to a balanced wet-lab
library. Wet-lab measurements remain the final evidence of promoter performance.

## 2. Differential-expression analysis

### Why we used public transcriptomic data

Testing every *E. coli* promoter experimentally would be impractical. Public
expression datasets let us identify genes that change during antibiotic
exposure or in resistant strains, then prioritize the upstream promoters most
likely to produce a detectable reporter response.

### Datasets and biological comparisons

| Antibiotic | Class | Current comparison | What the comparison tests | Main caveat |
|---|---|---|---|---|
| Amoxicillin | Beta-lactam | Resistant strain 512 vs wild type | Resistance-associated transcriptional state | Not a direct measure of acute induction in a susceptible strain |
| Amoxicillin | Beta-lactam | Resistant strain 512 + amoxicillin vs wild type + amoxicillin | Resistance-associated differences during matched exposure | Resistance genotype and exposure effects remain combined |
| Ceftazidime | Beta-lactam | Ceftazidime-treated vs water control | Acute response to beta-lactam exposure | Study-specific strain, dose, medium, and time point |
| Gentamicin | Aminoglycoside | Gentamicin-treated vs control | Acute aminoglycoside response | Microarray probes may not map one-to-one to genes |
| Tobramycin | Aminoglycoside | Tobramycin-treated vs control | Acute aminoglycoside response | Existing processed DE results are imported rather than re-estimated uniformly |

The direct wild-type amoxicillin-exposure comparison was initially examined but
was nearly flat. We therefore use the two stronger resistance-focused
comparisons above. This change improves candidate discovery for detecting a
resistance-associated state, but it changes the biological claim: these results
do not by themselves establish that a promoter is acutely induced by
amoxicillin.

### Standardization pipeline

The source studies use different platforms and file formats. A configuration-
driven Python pipeline converts them into a common schema:

```text
Public study files
     |
     v
Sample and condition metadata
     |
     v
Platform-appropriate expression table
     |
     v
Fold change and statistical testing
     |
     v
Benjamini-Hochberg false-discovery correction
     |
     v
Standard result columns and matched volcano plots
     |
     v
Upregulated, downregulated, and non-significant candidate tables
```

The standardized result records the antibiotic, class, biological comparison,
gene identifier, effect size, raw p-value, adjusted p-value, regulation class,
and the source of the adjusted p-value. Volcano plots use common axes so visual
differences between studies are not caused by rescaling.

### Current statistical method

The current repository performs differential-expression analysis, but it does
not run DESeq2 for every dataset:

- Count data are normalized by median-ratio size factors, log transformed, and
  tested using Welch's t-test.
- Replicate microarray expression values are tested using Welch's t-test.
- Tobramycin uses the differential-expression statistics supplied in its source
  workbook.
- Missing adjusted p-values are calculated using the Benjamini-Hochberg method.
- The current strict threshold is `|log2 fold change| > 2` and adjusted
  `p < 0.05`.

This workflow is reproducible and useful for exploration, but it is a simplified
analysis. Before presenting the statistics as final, raw RNA-seq counts should
be reanalyzed with DESeq2 or edgeR, while microarray datasets should be analyzed
with limma. These methods estimate gene-wise variance by sharing information
across genes and are better suited to small transcriptomic experiments than
independent t-tests.

### Current result summary

| Analysis | Tested rows | Upregulated | Downregulated |
|---|---:|---:|---:|
| Amoxicillin: resistant vs wild type | 4,237 | 10 | 16 |
| Amoxicillin: resistant + drug vs wild type + drug | 4,237 | 9 | 14 |
| Ceftazidime | 4,498 | 58 | 58 |
| Gentamicin | 10,208 | 283 | 480 |
| Tobramycin | 4,451 | 90 | 25 |

These row counts are not perfectly comparable. In particular, a microarray may
contain multiple probes for one gene. The number of significant rows should not
be interpreted as a direct measure of how strongly an antibiotic affects the
cell.

### Figures to place in this section

1. A dataset-condition table with strain, medium, dose, exposure time, platform,
   replicate count, and accession.
2. A sample QC panel for each dataset: expression distributions, sample
   correlations, and PCA or multidimensional scaling.
3. Four antibiotic volcano plots using the same x and y limits.
4. An UpSet plot showing gene-level overlap within beta-lactams, within
   aminoglycosides, and across all antibiotics.
5. A heatmap of shortlisted genes across all available samples.
6. A candidate table showing effect size, uncertainty, regulatory mechanism,
   specificity, and promoter feasibility.

### From genes to promoter candidates

A differentially expressed gene is not automatically a usable promoter. For
each gene, the next stage should determine:

- whether it is the first gene in a transcription unit;
- which promoter or promoters drive that transcription unit;
- promoter orientation and genomic coordinates;
- known transcription factors and sigma factors;
- whether regulation is activating or repressing;
- the likely upstream sequence to clone;
- whether the promoter is already present in the team's library.

Candidates should be retained in tiers rather than accepted or rejected by one
hard cutoff:

| Tier | Interpretation |
|---|---|
| 1 | Strong DE, coherent regulatory mechanism, usable promoter, and relevant comparison |
| 2 | Strong effect or cross-dataset support with weaker statistical evidence |
| 3 | Literature-supported stress promoter not recovered strongly in these datasets |
| Control | Positive, constitutive, promoterless, or known stress-response control |

## 3. Regulatory-network analysis

Regulatory-network analysis belongs between gene-level DE and final promoter
ranking. Its job is to test whether several targets controlled by the same
regulator respond coherently, and to identify the promoter that could reproduce
that response in a reporter construct.

### Proposed workflow

```text
Standardized DE tables
     |
     +--> map gene identifiers to E. coli K-12 genes
     |
     +--> attach transcription units, promoters, TFs, and sigma factors
     |
     +--> test regulon enrichment and directional coherence
     |
     +--> identify antibiotic-specific and class-conserved modules
     |
     v
Network-supported promoter candidates
```

RegulonDB and EcoCyc can provide curated *E. coli* transcription units,
promoters, regulators, and regulatory interactions. Curated and predicted
relationships should be labelled separately. Strain-specific genes that cannot
be mapped reliably to K-12 should remain flagged rather than silently dropped.

For each regulator, report:

- number of targets represented in the dataset;
- number and proportion changing in the expected direction;
- enrichment p-value and adjusted p-value;
- median target fold change;
- whether the response repeats across antibiotics or datasets;
- whether the network represents a specific drug mechanism or general stress.

This makes interpretation of controls such as `cpxP` and `recA` more rigorous.
`cpxP` is stronger evidence of a Cpx envelope-stress response when other Cpx
targets move coherently. `recA` is stronger evidence of SOS activation when
other LexA-regulated targets respond. Neither should be assumed to be a
beta-lactam-specific promoter without cross-condition testing.

## 4. Matching the lab promoter library

The computational shortlist should be joined to the team's physical inventory.
For each available part, record:

- promoter name and sequence;
- source strain and reference genome;
- transcription unit and regulated gene;
- expected regulator and stress pathway;
- plasmid backbone, copy number, RBS, GFP variant, and terminator;
- antibiotic datasets supporting the candidate;
- cloning status and sequence-verification status.

A balanced 80-promoter test library could contain:

| Candidate group | Suggested number |
|---|---:|
| Antibiotic-specific candidates, 12 per antibiotic | 48 |
| Class-conserved candidates, 6 per class | 12 |
| Broad stress or cross-class candidates | 8 |
| Literature and pathway controls, including `cpxP` and `recA` | 6 |
| Constitutive, promoterless, and non-responsive controls | 6 |
| **Total** | **80** |

The exact allocation should be updated after the available-parts list arrives.
Candidate selection should preserve regulatory and mechanistic diversity rather
than filling all 80 positions with genes from one strongly responding operon.

## 5. Hydrogel and in-chip transport models

These models establish the antibiotic exposure that the reporter actually sees.
They should be presented as connected but distinct components.

### Hydrogel diffusion and degradation

Estimate antibiotic concentration through the hydrogel as a function of
position and time. Relevant parameters include diffusion coefficient, hydrogel
thickness, partitioning, binding, degradation rate, temperature, and initial
sample concentration.

Primary design outputs:

- time required for the antibiotic to reach the cells;
- fraction of the starting concentration that reaches the sensing chamber;
- whether degradation reduces exposure below the reporter activation range;
- hydrogel thickness and composition compatible with the assay time.

### In-chip transport

Model advection, diffusion, mixing, and residence time in the chip geometry.
Relevant parameters include flow rate, channel dimensions, chamber volume,
inlet concentration, cell location, and adsorption to device materials.

Primary design outputs:

- concentration-time profile at the cell chamber;
- spatial uniformity of exposure;
- minimum incubation time before fluorescence modeling begins;
- chip geometries or flow rates that avoid underexposing the cells.

The output of both transport models should be a concentration-versus-time curve
that can be passed directly to the fluorescence model.

## 6. Fluorescence-response model

The fluorescence model should connect promoter activation to what the detector
measures. A minimal mechanistic model can include promoter activation, mRNA
production and degradation, GFP translation, GFP maturation, protein dilution,
cell growth, and instrument background.

```text
antibiotic concentration over time
                |
                v
       promoter activation
                |
                v
       mRNA -> immature GFP -> mature fluorescent GFP
                |
                v
    background-corrected detector signal
```

The wet-lab time course can be used to estimate promoter-specific parameters.
The most useful outputs are:

- fold induction over the unexposed control;
- absolute fluorescence per cell or normalized optical density;
- time to cross a detection threshold;
- dose-response midpoint and dynamic range;
- leakiness and false-positive rate;
- predicted limit of detection under chip conditions.

This model is more informative for device design than ranking promoters by
transcriptomic fold change alone. RNA abundance in the source study does not
directly determine GFP brightness in a new plasmid, host strain, medium, or
device.

## 7. Integrated model and sensitivity analysis

The strongest final dry-lab result would connect the modules:

```text
sample concentration
  -> transport/degradation
  -> cellular exposure
  -> promoter activation
  -> GFP production and maturation
  -> detector signal
  -> resistant/susceptible classification
```

Perform sensitivity analysis to identify which uncertain parameters most affect
time to detection and classification. This can guide experiments: for example,
the model may show that promoter leakiness matters more than maximum fold
induction, or that hydrogel diffusion dominates the assay time.

## 8. Design-Build-Test-Learn iterations

### Iteration 1: public-data discovery

- **Design:** identify antibiotic-responsive genes in public datasets.
- **Build:** implement a configuration-driven standardization and DE pipeline.
- **Test:** compare volcano plots, significant candidates, and overlaps.
- **Learn:** direct wild-type amoxicillin exposure was weak in the selected
  study, so resistance-focused comparisons were added and the claim was
  narrowed accordingly.

### Iteration 2: regulatory and physical feasibility

- **Design:** connect DE genes to transcription units, promoters, and regulators.
- **Build:** create a candidate table integrating network evidence and the lab
  promoter inventory.
- **Test:** select a mechanistically diverse 80-promoter panel.
- **Learn:** remove candidates that are redundant, unavailable, poorly mapped,
  or likely to report nonspecific stress.

### Iteration 3: quantitative reporter validation

- **Design:** test dose, time, antibiotic specificity, and resistant versus
  susceptible strains.
- **Build:** assemble promoter-GFP reporters with common genetic context.
- **Test:** measure background, dynamic range, response time, growth, and
  cross-reactivity.
- **Learn:** fit the fluorescence model and update the candidate ranking using
  measured performance.

## 9. Limitations and assumptions

- The source datasets differ in strain, platform, medium, dose, exposure time,
  controls, and preprocessing.
- The amoxicillin comparisons identify resistance-associated expression and do
  not isolate acute drug induction.
- Statistical significance does not guarantee a large or rapid GFP response.
- A gene may be downstream in an operon and may not have its own promoter.
- Regulatory annotations are most complete for *E. coli* K-12 and may not
  transfer perfectly to patient isolates.
- General stress promoters may be sensitive but produce false positives for
  unrelated stresses.
- Plasmid copy number and sequence context can change promoter behavior.
- Fixed volcano axes improve visual comparison but can hide extreme values;
  tables should preserve the full numeric results.
- Current rounded summary files should not be used when exact p-values are
  required.
- Wet-lab validation across doses, time points, media, and strains is required
  before making a diagnostic claim.

## 10. Software roadmap

### Highest priority before selecting the wet-lab library

1. **Platform-appropriate DE module.** Run DESeq2 or edgeR for raw RNA-seq
   counts and limma for microarrays; preserve study-specific design variables.
2. **Automated QC report.** Generate library-size checks, sample correlations,
   PCA or MDS, replicate variability, MA plots, and outlier flags.
3. **Candidate integration table.** Combine every dataset into one unrounded
   table with within-study percentile ranks, strict and relaxed evidence,
   overlaps, and explicit caveats.
4. **Regulatory-network enrichment.** Attach RegulonDB/EcoCyc interactions and
   score coherent TF and sigma-factor regulons.
5. **Promoter-context mapper.** Map genes to transcription units, extract
   promoter coordinates and sequences, and flag ambiguous or multi-promoter
   operons.
6. **Library selection optimizer.** Choose approximately 80 candidates while
   balancing antibiotic specificity, class overlap, pathway diversity,
   availability, and controls.

### High-value modeling additions

7. **Fluorescence kinetics fitter.** Fit promoter dose-response and time-course
   data and estimate time to detection and dynamic range.
8. **Coupled transport-reporter simulator.** Feed hydrogel and chip exposure
   curves into the fluorescence model.
9. **Sensitivity and uncertainty analysis.** Quantify how parameter uncertainty
   changes device predictions and wet-lab priorities.
10. **Cross-strain robustness module.** Check promoter sequence conservation and
    response consistency across laboratory and clinical *E. coli* backgrounds.

### Useful presentation software

11. **Interactive candidate explorer.** Filter candidates by antibiotic,
    regulator, effect size, adjusted p-value, specificity, and library status.
12. **Automated wiki figure export.** Rebuild publication-ready plots, result
    tables, and method metadata from one versioned analysis run.
13. **Reproducibility report.** Record source accession, file checksum, software
    version, configuration, and analysis date for every result.

An elaborate machine-learning classifier should wait until the team has enough
wet-lab measurements. With approximately 80 constructs, interpretable kinetic
models and held-out validation are likely to be more defensible than a complex
black-box model.

## 11. Recommended page order

Use this order on the final wiki:

1. Project question and integrated pipeline figure
2. How computation changed device and wet-lab decisions
3. Differential-expression and promoter-selection analysis
4. Regulatory-network and promoter-context analysis
5. Lab promoter library design
6. Hydrogel diffusion and degradation model
7. In-chip transport model
8. Fluorescence model
9. Integrated prediction and sensitivity analysis
10. DBTL iterations
11. Limitations and future work
12. Code, data provenance, and references

This follows the useful structure of McGill iGEM's 2025 model page: objective,
method, dataset construction, evaluation, DBTL, limitations, and future
extensions. The content and claims here remain specific to our antibiotic
biosensor project.

## References and implementation resources

- [McGill iGEM 2025 model page](https://2025.igem.wiki/mcgill/model)
- [DESeq2 Bioconductor documentation](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- [limma Bioconductor documentation](https://bioconductor.org/packages/release/bioc/html/limma.html)
- [RegulonDB 11 overview](https://pmc.ncbi.nlm.nih.gov/articles/PMC9465075/)
- [EcoCyc overview](https://pmc.ncbi.nlm.nih.gov/articles/PMC6504970/)
- [Amoxicillin study GSE47221](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47221)
- [Ceftazidime study GSE220559](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220559)
- [Gentamicin study GSE44211](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44211)
- [Tobramycin study GSE224240](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224240)

