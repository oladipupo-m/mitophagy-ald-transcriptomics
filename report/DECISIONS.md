# Processing decisions — public-data mitophagy / ALD differential expression

Every analytical decision behind the saved differential expression, enrichment and
specificity outputs, with its value, its justification, and where it is recorded in code.
Reconstructed from the saved script, the saved outputs and the session record. Items that
could **not** be recovered are listed in §10 rather than guessed.

**Outputs this document covers:** `de_all_contrasts.csv.gz`, `mitophagy_gene_results.csv`,
`gsea_mitophagy_modules.csv`, `gsea_reference_sets.csv`,
`mitophagy_vs_mitochondrial_background.csv`, `replicated_mitophagy_genes.csv`,
`key_mitophagy_genes.csv`, and figures 1–4.

Decision classes are marked **[R]** rule (one correct answer), **[J]** judgement call
(defensible alternatives exist, choice must be justified), **[C]** constraint (forced by
the platform, not chosen).

---

## 1. Cohort selection

| Decision | Value | Class | Justification |
|---|---|---|---|
| Cohorts | GSE28619, GSE143318, GSE103580, GSE270042 | **J** | Two case–control cohorts on independent platforms give cross-platform replication; one severity gradient tests stage-dependence; one paired within-patient contrast removes between-subject confounding |
| Species / tissue | Human liver only | **R** | Question is about human ALD; rodent models answer a different question |
| GSE270042 scope | SubSeries of GSE270043 (ISAIAH trial, NCT03775109) | **C** | The parent is a SuperSeries mixing assay types |
| Excluded | GSE94417, GSE95246, GSE155926 | **J** | SuperSeries mixing assay types, or no usable control group |

## 2. Matrix preparation

Full provenance in `README_matrices.md`. Decisions summarised:

| Decision | Value | Class | Justification |
|---|---|---|---|
| GSE28619 scale | used as supplied (log2, range 2.3–15.5) | **R** | Series matrix is already RMA-processed |
| GSE103580 scale | log2(x + 1) applied | **R** | Supplied on a linear scale (range 1.1–16,699) |
| Probe → gene collapse | one probe per gene, **maximum mean expression** | **J** | Common convention; maximum variance or median are defensible alternatives and would give slightly different matrices |
| GSE143318 duplicate rows | deduplicated on ENSEMBL id before summing to symbol | **R** | File contains one row per ENSEMBL–ENTREZ pair; summing without deduplication inflates counts up to 44,521-fold for some genes |
| GSE270042 non-integer counts | rounded to integer | **C** | Estimated counts from transcript quantification; DESeq2 requires integers |
| GSE270042 pairing | 32 patients with both timepoints (64 of 65 samples) | **R** | One patient has day 0 only; a paired model cannot use it |

## 3. Gene set construction

| Decision | Value | Class | Justification |
|---|---|---|---|
| Sources | GO via QuickGO (GO:0000422, GO:1903146, descendants) + Reactome (R-HSA-5205647, R-HSA-5205685, R-HSA-8934903) + a specified mechanistic core | **J** | Ontology-defined sets are reproducible; the manual core guarantees the PINK1/Parkin axis is represented even where annotation is sparse |
| GO relationship types | `is_a`, `part_of`, `regulates` | **J** | Including `regulates` admits regulators that are not pathway components — this widens the set and is the reason modules are tested separately |
| GO:0098779 | requested, returned no convertible results | **C** | QuickGO rejected the identifier; the term contributes nothing to the final set |
| Modules | 5: PINK1–Parkin core (21), receptor-mediated (10), autophagy machinery (14), mitochondrial dynamics/biogenesis (9), GO/Reactome extended (100) | **J** | Lets a mechanistic core be tested separately from broad annotation |
| Cleaning | dropped `PARK2`, `DKFZp762B153`, `HEL-S-7a`; UniProt accessions filtered by regex | **R** | `PARK2` is the retired symbol for *PRKN*, already present — keeping both double-counts one gene. The other two are not gene symbols |
| Final size | 157 → **154 genes** | — | — |

## 4. Identifier reconciliation

Gene-set members absent from a matrix were queried against NCBI Gene (`esearch`/`esummary`
on the `gene` database, `[Preferred Symbol]`, *Homo sapiens*) and accepted only where the
matrix symbol is a **retired** identifier.

**Accepted (14):** `ATPIF1`→*ATP5IF1*, `KIAA0141`→*DELE1*, `GBAS`→*NIPSNAP2*,
`PARK2`→*PRKN*, `C5orf51`→*RIMOC1*, `TOMM70A`→*TOMM70*, `C12orf44`→*ATG101*,
`KIAA0652`→*ATG13*, `KIAA0831`→*ATG14*, `LASS1`→*CERS1*, `SMCR7L`→*MIEF1*,
`MTERFD1`→*MTERF3*, `C16orf55`→*SPATA33*, `WDR45L`→*WDR45B*.

**Rejected (1):** `CIR1`→*UBE2V1*. `CIR1` is the current approved symbol of a different
gene, so the alias is ambiguous and the mapping was not applied. **[R]**

Renaming applied per dataset, and only where the current symbol was absent from that matrix,
so no merge collisions were created. Without this step *PRKN* and 13 others are silently
missing from both microarray cohorts, and an earlier version of this analysis wrongly
reported *PRKN* as undetectable on arrays.

## 5. Differential expression

Genes retained and tested per contrast:

| Dataset | Contrast | Model | Genes tested |
|---|---|---|---|
| GSE28619 | AH vs Normal | two-sample *t*, equal variance | 21,445 |
| GSE103580 | severity trend | OLS slope on ordinal severity | 19,469 |
| GSE143318 | AH vs Normal | DESeq2 Wald | 18,898 |
| GSE143318 | AlcCirrhosis vs Normal | DESeq2 Wald | 18,898 |
| GSE270042 | day28 vs day0, paired | DESeq2 Wald | 23,593 |

| Decision | Value | Class | Justification |
|---|---|---|---|
| Microarray test | `scipy.stats.ttest_ind(..., equal_var=True)` | **C** | limma's empirical-Bayes moderation was unavailable (see §11). Unmoderated statistics are valid but less powerful at n = 15 vs 7 |
| Severity coding | Steatosis 1, MildAH 2, Cirrhosis 3; OLS slope | **J** | A trend test across ordered stages is more powerful than three pairwise contrasts and matches the biological ordering. It assumes equal spacing between stages, which is an assumption, not a fact |
| RNA-seq model | PyDESeq2, negative-binomial GLM, Wald test, `refit_cooks=True` | **R** | Field standard for count data |
| Count filter | `(counts >= 10).sum(axis=1) >= max(5, 0.2 × n_samples)` → threshold 5 for GSE143318 (n=25), 12 for GSE270042 (n=64) | **J** | Independent of the outcome, so it does not bias the test; the 20% rule scales with cohort size |
| Paired design | `~patient + day` | **R** | The patient term absorbs between-subject variation; omitting it inflates residual variance and loses the within-patient contrast |
| Contrast direction | alt vs reference, reference named second | **R** | Positive log2FC = higher in the first-named group |

## 6. Multiplicity

Benjamini–Hochberg FDR (`statsmodels.stats.multitest.multipletests`, `method="fdr_bh"`),
applied **per contrast across all genes with a finite p-value**, recomputed in the harness
rather than taken from DESeq2's own independent-filtering-adjusted values, so that every
contrast is corrected identically. **[J]** — DESeq2's independent filtering would give
slightly more power on the RNA-seq contrasts; uniformity across platforms was preferred.

## 7. Gene-set enrichment

`gseapy.prerank`, ranked by the per-contrast **test statistic** (not fold change).

| Parameter | Modules | Reference sets |
|---|---|---|
| `permutation_num` | 10,000 | 5,000 |
| `min_size` | 8 | 5 |
| `max_size` | 2,000 | 2,000 |
| `seed` | 0 | 0 |
| `threads` | 4 | 4 |

Ranking by test statistic rather than fold change **[R]** — fold change ignores
within-group variance and over-ranks low-expression genes.

## 8. Reference sets and the competitive test

Reference sets, built to have known composition behaviour, used as negative controls:

| Set | Definition | n (GSE103580) |
|---|---|---|
| Hepatocyte identity | explicit list: ALB, APOB, APOA1, APOA2, TTR, SERPINA1, TF, FGA, FGB, FGG, CYP2E1, CYP3A4, CYP1A2, HNF4A, ASGR1, CPS1, TAT, G6PC1, PCK1 | 19 |
| OXPHOS (nuclear-encoded) | regex `^(NDUF[ABCSV]|SDH[ABCD]$|UQCR|COX[45678]|ATP5[FMPO])` | 126 |
| Mitoribosome | regex `^MRP[LS]\d+$` | 75 |
| Fibrosis / ECM | explicit list: COL1A1, COL1A2, COL3A1, ACTA2, TIMP1, LOX, FN1, PDGFRB, TGFB1, SPP1, THBS2 | 11 |
| Immune / inflammatory | explicit list: PTPRC, CD68, CXCL8, IL1B, CCL2, LCN2, S100A8, S100A9, CXCL6, TNF | 10 |

**Competitive test.** Mann–Whitney U (two-sided) comparing the test statistics of mitophagy
genes against a **mitochondrial background** = OXPHOS ∪ Mitoribosome (148–185 genes
depending on cohort). BH-corrected across the ten tests. **[J]** — testing against a
biologically matched background rather than against all genes is what licenses the claim
that the mitophagy signal is not specific; testing against all genes would answer a
different and weaker question.

The regex sets are **approximations of their biology**, not curated gene sets. `^MRP[LS]\d+$`
captures mitochondrial ribosomal proteins well; the OXPHOS pattern captures most
nuclear-encoded subunits of complexes I and III–V but not all, and excludes
mitochondrially-encoded subunits entirely. A curated source (MSigDB `HALLMARK_OXIDATIVE_
PHOSPHORYLATION`, or MitoCarta) would be the stronger choice and is the first thing to
change on a re-run.

## 9. Replication criteria

A gene counts as replicated if it reaches FDR < 0.05 in **at least three of the four disease
contrasts** with a **consistent direction** in all of them. **[J]** — three of four is a
convention chosen here, not a standard; two of four would be permissive, four of four would
discard genes absent from one platform.

## 10. Decisions not recoverable from the saved record

- **Exact package versions at run time.** The kernel has been reset since. Current versions
  in the `mitode2` environment are: pydeseq2 0.5.4, gseapy 1.3.1, scipy 1.17.1,
  statsmodels 0.15.0, pandas 2.3.3, numpy 2.4.6. These are *current*, not verified
  run-time values.
- **QuickGO and Reactome release versions** at the time of retrieval. Retrieval date is
  recoverable from the session record; the database release is not.
- **The exact dataset-screening trail.** Candidate accessions were inspected interactively;
  only the final selection is recorded, so the excluded-with-reason log that a PRISMA-style
  account needs does not exist for the dataset search.

## 11. Known deviations from field standard

1. **No limma.** Microarray contrasts use unmoderated t-statistics because Bioconductor is
   not distributed through conda for Windows and the sandbox refuses to load compiled R
   libraries from a writable path. This costs power at small n and is the first thing to
   re-run on a Linux machine.
2. **No deconvolution.** Cell composition is *controlled for* by reference sets, not
   *adjusted for* by estimating hepatocyte fraction. The conclusion is therefore comparative
   ("mitophagy moves less than the mitochondrial background") rather than absolute.
3. **FDR recomputed outside DESeq2**, forgoing its independent filtering (see §6).
4. **GSE143318 counts were depth-normalised before deposit** (library size CV 1.1%), so
   DESeq2's size factors are near 1 and its normalisation assumptions are only partly met.

## 12. Code locations

| Stage | Where |
|---|---|
| Differential expression, module GSEA | `de_analysis.py` |
| Matrix provenance | `README_matrices.md` |
| Retrieval, gene-set construction, alias reconciliation, reference sets, competitive test, figures | **not in a script** — run interactively; this document is the only record |

Closing that last gap — consolidating everything into one pipeline — is the outstanding
reproducibility task.
