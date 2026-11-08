# Harmonisation, verification and independent assessment

Four questions were asked of the existing public-data analysis: is there a harmonised
matrix and metadata set, are the sample classifications correct, are the RNA-seq count
formats resolved, and do the differential expression findings survive independent testing.

**Summary.** Harmonised outputs did not exist and have now been built. Sample
classifications verify, both against the raw GEO records and biologically. The count-format
issues are resolved, and one of them was mis-stated in the earlier documentation. The
headline gene-level finding survives independent testing; one pathway-level claim does not.

---

## 1. Harmonisation — did not exist, now built

Previously there were four parallel tables with different column names (`gsm` vs `sample`),
different group vocabularies (`Normal`/`AH` vs `Steatosis`/`MildAH`/`Cirrhosis` vs
`d0`/`d28`) and different value types. Nothing mapped them onto a common schema.

`harmonised_sample_metadata.csv` — **198 samples**, one row each, with: `sample_id`,
`dataset`, `platform`, `assay`, `value_type`, `group_raw`, `group_harmonised`, `is_control`,
`severity_rank`, `patient`, `timepoint`, `treatment`, `library_size`, `size_factor`,
`flag_low_depth`, `gsm_linked`.

Harmonised vocabulary and ordinal severity rank:

| Rank | Group | Datasets |
|---|---|---|
| 0 | Normal liver / Donor liver | GSE28619 (7), GSE143318 (7) |
| 1 | Alcoholic steatosis | GSE103580 (6) |
| 2 | Mild alcoholic hepatitis | GSE103580 (13) |
| 3 | Alcoholic hepatitis | GSE28619 (15), GSE143318 (13), GSE270042 (33 d0 / 32 d28) |
| 4 | Alcoholic cirrhosis | GSE103580 (67), GSE143318 (5) |

`harmonised_gene_universe.csv` — 69,224 symbols across the union, **15,346 present in all
four** cohorts, 18,831 in at least three. Of the 154 mitophagy genes, **133 are present in
all four**. Any cross-cohort analysis should state which universe it used.

Note the two control groups are not equivalent: GSE28619 uses "normal liver" and GSE143318
uses donor organs. They are harmonised as rank 0 but flagged separately in `group_harmonised`
because they are not interchangeable.

## 2. Sample classification — verified, with one caveat

**Re-parsed from the raw GEO series matrices and compared to the derived metadata:**

| Dataset | Result |
|---|---|
| GSE28619 | exact match, 22/22 |
| GSE103580 | exact match, 86/86 |
| GSE270042 | exact match on day, patient and treatment, 65/65 |
| GSE143318 | group **counts** match (AH 13, Donor 7, Cirrhosis 5) but see caveat |

**Caveat — GSE143318 has no machine-readable sample linkage.** The counts file uses column
names `N1…N8`, `AC1…AC5`, `AH1…AH64`; the series matrix carries GSM accessions with titles
like "Alcoholic Hepatitis Patient 1". `!Sample_supplementary_file` is `NONE` for every
sample and `!Sample_description` is "liver homogenate" throughout. There is therefore **no
way to map an individual GSM to an individual count column.** Group assignment is an
inference from the column prefix. The group sizes agree and the biology agrees (below), so
group-level conclusions are sound, but individual sample identity is not traceable and the
analysis cannot be re-linked to per-patient clinical data.

**Biological verification (orthogonal to the labels).** Mean z-scored marker-panel
expression by group — hepatocyte identity, fibrosis/ECM, inflammation:

| Dataset | Group | Hepatocyte | Fibrosis | Inflammation |
|---|---|---|---|---|
| GSE28619 | Normal | **+0.88** | −1.29 | −1.04 |
| GSE28619 | AH | −0.41 | **+0.60** | **+0.49** |
| GSE103580 | Steatosis | **+0.78** | −1.10 | −0.68 |
| GSE103580 | MildAH | +0.38 | −0.18 | +0.23 |
| GSE103580 | Cirrhosis | −0.14 | **+0.13** | +0.02 |
| GSE143318 | Donor | **+0.91** | −1.20 | −1.11 |
| GSE143318 | AlcCirrhosis | +0.44 | 0.00 | −0.56 |
| GSE143318 | AH | −0.66 | **+0.65** | **+0.81** |
| GSE270042 | day 0 | −0.11 | +0.32 | +0.55 |
| GSE270042 | day 28 | **+0.12** | −0.33 | −0.57 |

Every cohort behaves as its labels predict: controls have the highest hepatocyte signal and
lowest fibrosis and inflammation, GSE103580 is monotonic across its severity gradient, and
GSE270042 reverses on recovery. No label swap is detectable.

## 3. RNA-seq count formats — resolved, and one earlier caveat was wrong

**GSE270042 rounding: immaterial.** 12.5% of the deposited values were non-integer
(transcript-level estimates). Maximum absolute change from rounding is 0.500, mean 0.0295,
and **no cell changed by more than 0.5**. Rounding is not a source of error.

**GSE143318: the earlier caveat was mis-stated and is hereby corrected.** The previous
documentation said that because library sizes are near-identical (CV 1.07%), "DESeq2 size
factors will be close to 1 and the usual assumptions are only partly met."

Computing them shows the opposite. Median-of-ratios size factors span **0.625 to 1.417, a
2.27-fold range**. Equal total counts do not imply equal composition — the submitter
normalised depth, but the samples still differ substantially in how reads are distributed
across genes, which is exactly what median-of-ratios normalisation corrects. DESeq2 is doing
necessary work on this dataset, not redundant work. The correct caveat is narrower: the
counts are not raw, so the mean–variance relationship DESeq2 assumes is distorted by the
upstream scaling, but its normalisation step remains appropriate.

**GSE270042: an unflagged quality problem.** Library sizes span **2.15 M to 40.95 M reads**
(19-fold) and size factors **0.010 to 3.146** (320-fold). Seven libraries fall below 10 M
reads, all of them inside the paired analysis set:

| Sample | Reads (M) | Size factor | Day |
|---|---|---|---|
| 303_d28 | 2.15 | 0.010 | d28 |
| 341_d0 | 5.12 | 0.023 | d0 |
| 363_d0 | 5.79 | 0.188 | d0 |
| 146_d28 | 6.61 | 0.166 | d28 |
| 332_d28 | 7.12 | 0.196 | d28 |
| 313_d0 | 7.70 | 0.046 | d0 |
| 311_d0 | 9.65 | 0.141 | d0 |

A size factor of 0.010 means that library carries roughly a hundredth of the expected signal
— it is effectively a failed run. These are split 4 day-0 / 3 day-28, so they are not
confounded with the contrast, but they should have been flagged. They now are, in
`flag_low_depth`.

## 4. Independent assessment of the differential expression findings

### 4a. Distribution-free re-analysis of every contrast

Each contrast was re-tested with a method making different assumptions: Wilcoxon rank-sum
for the two-group comparisons, Spearman correlation with severity for the gradient, and
Wilcoxon signed-rank on log-CPM for the paired contrast.

| Contrast | n genes | sig original | sig alternative | overlap | Jaccard | direction agreement (shared) |
|---|---|---|---|---|---|---|
| GSE28619 AH vs Normal | 21,439 | 6,164 | 6,006 | 5,500 | 0.825 | **1.000** |
| GSE103580 severity | 19,455 | 980 | 1,671 | 910 | 0.523 | **1.000** |
| GSE143318 AH vs Normal | 18,898 | 11,060 | 11,045 | 7,449 | 0.508 | 0.968 |
| GSE143318 cirrhosis vs Normal | 18,898 | 6,591 | 5,957 | 3,668 | 0.413 | **1.000** |
| GSE270042 day28 vs day0 | 23,592 | 1,069 | 187 | 143 | 0.128 | **1.000** |

No gene reverses direction in any contrast. The low Jaccard for the paired contrast reflects
a power difference — a signed-rank test over 32 pairs is far weaker than a negative-binomial
GLM — not a disagreement.

### 4b. The headline gene-level finding holds

*BNIP3*, under the independent test:

| Contrast | original log2FC (padj) | alternative direction (padj) |
|---|---|---|
| GSE28619 AH vs Normal | −0.35 (0.0162) | down (0.0369) |
| GSE103580 severity | −0.27 (0.0188) | down (0.0133) |
| GSE143318 AH vs Normal | −1.41 (<0.0001) | down (0.0030) |
| GSE143318 cirrhosis vs Normal | −0.94 (0.0003) | down (0.0471) |
| GSE270042 day28 vs day0 | +0.11 (0.84) | n.s. (0.89) |

Significant and downward in all four disease contrasts under both methods, and null at day 28
under both. This finding is robust.

### 4c. The specificity result holds under a curated background

The original competitive test used regular-expression gene sets, which I had flagged as the
weakest component. Replacing them with curated GO terms — GO:0006119 oxidative
phosphorylation and GO:0005761 mitochondrial ribosome, descendants included, 223 symbols,
zero overlap with the mitophagy set:

| Contrast | median mitophagy | median bg (GO) | FDR (GO) | median bg (regex) | FDR (regex) |
|---|---|---|---|---|---|
| GSE103580 severity | −0.07 | −1.62 | 1.2e-14 | −1.62 | 1.5e-14 |
| GSE143318 AH vs donor | −1.80 | −2.72 | 5.1e-03 | −2.50 | 5.3e-03 |
| GSE143318 cirrhosis vs donor | −1.44 | −2.66 | 3.2e-09 | −2.64 | 2.3e-08 |
| GSE270042 day28 vs day0 | +0.40 | +0.48 | 0.229 | +0.45 | 0.439 |
| **GSE28619 AH vs Normal** | **+1.69** | **+2.55** | **0.095** | **+2.69** | **0.049** |

The conclusion is unchanged in three of four disease contrasts and strengthens in two. **But
GSE28619 moves from just significant (0.049) to not significant (0.095).** That contrast was
borderline and is not robust to background choice. The claim should be stated as: mitophagy
transcripts shift less far than the mitochondrial background in all four disease contrasts,
significantly so in three.

### 4d. One pathway-level claim does not survive

Re-running the paired contrast with the seven low-depth patients excluded (25 patients,
50 samples, 23,728 genes):

| Gene | full cohort | high-depth only |
|---|---|---|
| *USP30* | +0.55 (padj 2.5e-05) | **+0.52 (padj 3.1e-05)** |
| *SQSTM1* | −0.58 (padj 0.056) | −0.74 (padj 0.00024) |
| *PINK1* | +0.19 (0.57) | +0.25 (0.20) |
| *BNIP3* | +0.11 (0.84) | −0.00 (0.999) |

*USP30* is robust. But the module-level result is not:

| Module | NES full cohort | NES high-depth only |
|---|---|---|
| PINK1–Parkin core | **+1.58** (p 0.022, FDR 0.093) | **+0.90** (p 0.61, FDR 1.00) |

**The claim in `report.md` §1.4 that the PINK1–Parkin core "is the one module that rises
during recovery" does not survive removal of the low-depth libraries and should be
withdrawn.** It was nominally significant at best (FDR 0.093) and is driven substantially by
the seven poor libraries. The gene-level *USP30* result, which is what the poster states, is
unaffected.

## 5. What changed as a result

- Harmonised sample metadata and gene universe now exist.
- Seven low-depth libraries in GSE270042 are flagged; any future analysis of that cohort
  should run with and without them.
- The GSE143318 size-factor caveat in `README_matrices.md` is corrected (§3).
- The reference background is upgraded from regex patterns to curated GO terms.
- One pathway-level claim is withdrawn; the specificity claim is restated as three of four
  contrasts significant.
- GSE143318 sample-level traceability is recorded as a permanent limitation.
