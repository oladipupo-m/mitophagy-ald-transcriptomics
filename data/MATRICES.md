# Expression matrices: provenance

No raw expression matrix is committed to this repository. Each is retrieved from
GEO by `code/download_data.py`. The notes below record exactly which file was used
for each cohort, and how it was processed into the harmonised form in `metadata/`.

---

# Analysis-ready expression matrices and sample metadata

Four public human liver cohorts from NCBI GEO, processed to gene-level matrices for the
mitophagy / alcohol-related liver disease study. These are **derived** files; the raw
downloads are re-fetchable from GEO at the accessions below.

## Files

| File | Contents | Rows x cols |
|---|---|---|
| `GSE28619_expr.csv` | log2 expression, gene-level | 22,880 x 22 |
| `GSE28619_meta.csv` | sample, group (AH 15 / Normal 7) | 22 x 3 |
| `GSE103580_expr.csv` | log2 expression, gene-level | 19,469 x 86 |
| `GSE103580_meta.csv` | sample, disease status (Steatosis 6 / MildAH 13 / Cirrhosis 67) | 86 x 2 |
| `GSE143318_counts.csv` | integer counts, gene-level | 35,934 x 25 |
| `GSE143318_meta.csv` | sample, group (AH 13 / AlcCirrhosis 5 / Normal 7) | 25 x 2 |
| `GSE270042_counts.csv` | integer counts, gene-level | 59,452 x 65 |
| `GSE270042_meta.csv` | sample, patient, day, treatment (d0 33 / d28 32) | 65 x 4 |
| `GPL570_probe2symbol.csv` | probe to gene symbol map | 45,782 x 2 |
| `GPL13667_probe2symbol.csv` | probe to gene symbol map | 49,384 x 2 |

In every matrix, column order matches metadata row order (verified).

## Provenance and processing state

**GSE28619** — Affymetrix HG-U133 Plus 2.0 (GPL570). Source: GEO series matrix, which is
already RMA-processed and log2-scaled (observed range 2.3-15.5). Probes mapped to symbols via
the GPL570 platform table and collapsed one probe per gene by **maximum mean expression**.

**GSE103580** — Affymetrix HG-U219 (GPL13667). Source: GEO series matrix, supplied on a
**linear** scale (range 1.1-16,699); log2(x+1) transformed here before probe collapse. Same
collapse rule.

**GSE143318** — RNA-seq. Source: author supplementary file `GSE143318_Rawcount.txt.gz`. The
file contains multiple rows per gene from one-to-many ENSEMBL-to-ENTREZ mappings; deduplicated
on ENSEMBL id before summing to gene symbol. **Caveat:** library sizes are near-identical
across samples (CV 1.1%), so these counts appear to have been depth-normalised before deposit.
DESeq2 size factors will be close to 1 and the usual assumptions are only partly met.

**GSE270042** — RNA-seq, SubSeries of GSE270043 (ISAIAH trial, NCT03775109). Source: author
supplementary counts matrix. Index was `ENSG_SYMBOL`; split on the first underscore and summed
to symbol. Estimated counts were non-integer and have been **rounded**. All 65 samples are
included here; the paired analysis used only the 32 patients with both timepoints (64 samples).

## Known issues to carry forward

1. **Legacy symbols.** Both microarray platforms annotate *PRKN* as `PARK2` and *ATG101* as
   `C12orf44`. Fourteen gene-set members are only reachable after alias reconciliation against
   NCBI Gene. Do this before any gene-set analysis.
2. The 72 h block of an unrelated file in this project uses uM where mM was meant — unrelated
   to these matrices, noted only to avoid confusion.
3. Probe collapse by maximum mean expression is a defensible convention, not the only one.
   Collapsing by maximum variance or by median gives slightly different matrices.

## Accessions

GSE28619, GSE103580, GSE143318, GSE270042 (all NCBI GEO, public, no access restrictions).
