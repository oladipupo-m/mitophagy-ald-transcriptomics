# Code

| Script | Covers |
|---|---|
| `download_data.py` | Retrieval of every raw file from GEO into `data/raw/` |
| `de_analysis.py` | Differential expression (PyDESeq2 / moderated t) and module GSEA |
| `reference_sets_gsea.py` | Curated GO reference sets, competitive context, Figure 1 |

## What is not yet scripted

Gene-set construction, gene-symbol alias reconciliation, probe-to-symbol collapse and
the Mann–Whitney competitive test were run interactively. `report/DECISIONS.md` records
every parameter and every judgement call from those steps, but they cannot currently be
re-executed from this repository alone.

This is the honest state of the work, and consolidating those steps into the pipeline is
the first item of outstanding work. Until then, the committed tables in `results/` are
the record of what those steps produced.

## Parameters worth knowing

- **RNA-seq gene filter**: `(counts >= 10).sum(axis=1) >= max(5, 0.2 * n_samples)`
- **Designs**: `~group`; `~patient + day` for the paired GSE270042 contrast
- **PyDESeq2**: `refit_cooks=True`
- **GSEA**: `gseapy.prerank`, ranked by test statistic, 10,000 permutations for modules
  and 5,000 for reference sets, `min_size` 8 / 5, `max_size` 2000, `seed=0`
- **FDR**: Benjamini–Hochberg recomputed outside DESeq2 so that microarray and RNA-seq
  contrasts are corrected identically
- **Probe collapse**: maximum mean expression across samples
