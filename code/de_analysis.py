"""Mitophagy gene mining across four public alcohol-related liver disease cohorts.

Array cohorts (GSE28619, GSE103580): ordinary least-squares per gene.
RNA-seq cohorts (GSE143318, GSE270042): PyDESeq2 (negative-binomial GLM, Wald test).
Module-level enrichment: gseapy.prerank on the per-contrast test statistic.
"""
import os, json
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import gseapy

os.makedirs("results", exist_ok=True)
RNG = 0

gs = pd.read_csv("data/R/mitophagy_geneset.csv")
modules = {m: sorted(set(g.symbol)) for m, g in gs.groupby("module")}
modules["Mitophagy (all)"] = sorted(set(gs.symbol))
modules["Mitophagy (GO core)"] = sorted(set(gs.loc[gs.in_GO_mitophagy, "symbol"]))


def bh(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if ok.sum():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def ols_two_group(expr, groups, ref, alt):
    """Per-gene two-sample t-test; returns log2FC (alt - ref), t, p."""
    a = expr.loc[:, groups == alt].to_numpy(float)
    b = expr.loc[:, groups == ref].to_numpy(float)
    t, p = stats.ttest_ind(a, b, axis=1, equal_var=True)
    return pd.DataFrame({"symbol": expr.index, "log2FC": a.mean(1) - b.mean(1),
                         "stat": t, "pvalue": p})


def ols_trend(expr, x):
    """Per-gene OLS slope of expression on an ordinal covariate."""
    Y = expr.to_numpy(float)
    x = np.asarray(x, float)
    xc = x - x.mean()
    Sxx = (xc ** 2).sum()
    beta = (Y * xc).sum(1) / Sxx
    alpha = Y.mean(1) - beta * x.mean()
    resid = Y - (alpha[:, None] + beta[:, None] * x[None, :])
    dfree = len(x) - 2
    se = np.sqrt((resid ** 2).sum(1) / dfree / Sxx)
    t = beta / se
    p = 2 * stats.t.sf(np.abs(t), dfree)
    return pd.DataFrame({"symbol": expr.index, "log2FC": beta, "stat": t, "pvalue": p})


def run_deseq(counts, meta, design, contrast, label, sample_col):
    meta = meta.set_index(sample_col).loc[counts.columns]
    keep = (counts >= 10).sum(1) >= max(5, int(0.2 * counts.shape[1]))
    dds = DeseqDataSet(counts=counts.loc[keep].T, metadata=meta, design=design,
                       refit_cooks=True, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=contrast, quiet=True)
    st.summary()
    r = st.results_df.reset_index().rename(
        columns={"index": "symbol", "log2FoldChange": "log2FC", "padj": "padj_src"})
    r["symbol"] = r.iloc[:, 0]
    return dds, r[["symbol", "log2FC", "stat", "pvalue"]].copy()


de_tables, gsea_rows = [], []


def finalise(df, dataset, contrast):
    df = df.dropna(subset=["pvalue"]).copy()
    df["padj"] = bh(df["pvalue"])
    df.insert(0, "contrast", contrast)
    df.insert(0, "dataset", dataset)
    de_tables.append(df)
    rank = df.dropna(subset=["stat"]).drop_duplicates("symbol").set_index("symbol")["stat"].sort_values(ascending=False)
    pre = gseapy.prerank(rnk=rank, gene_sets=modules, permutation_num=10000,
                         min_size=8, max_size=2000, seed=RNG, no_plot=True,
                         outdir=None, threads=4)
    g = pre.res2d.copy()
    g["dataset"], g["contrast"] = dataset, contrast
    gsea_rows.append(g)
    return df


# ---------- 1. GSE28619: AH vs normal liver (Affymetrix GPL570, log2) ----------
e1 = pd.read_csv("data/R/GSE28619_expr.csv", index_col=0)
m1 = pd.read_csv("data/R/GSE28619_meta.csv").set_index("gsm").loc[e1.columns]
finalise(ols_two_group(e1, m1["group"].values, "Normal", "AH"),
         "GSE28619", "AH vs Normal")

# ---------- 2. GSE103580: ALD severity trend (Affymetrix GPL13667, log2) ----------
e2 = pd.read_csv("data/R/GSE103580_expr.csv", index_col=0)
m2 = pd.read_csv("data/R/GSE103580_meta.csv").set_index("gsm").loc[e2.columns]
sev = m2["status"].map({"Steatosis": 1, "MildAH": 2, "Cirrhosis": 3}).values
finalise(ols_trend(e2, sev), "GSE103580", "ALD severity trend")

# ---------- 3. GSE143318: AH / alcoholic cirrhosis vs donor liver (RNA-seq) ----------
c3 = pd.read_csv("data/R/GSE143318_counts.csv", index_col=0)
m3 = pd.read_csv("data/R/GSE143318_meta.csv")
for alt in ["AH", "AlcCirrhosis"]:
    dds3, r3 = run_deseq(c3, m3.copy(), "~group", ["group", alt, "Normal"],
                         f"{alt} vs Normal", "sample")
    finalise(r3, "GSE143318", f"{alt} vs Normal")

# ---------- 4. GSE270042: paired AH day 28 vs day 0 (RNA-seq, ISAIAH trial) ----------
c4 = pd.read_csv("data/R/GSE270042_counts.csv", index_col=0)
m4 = pd.read_csv("data/R/GSE270042_meta.csv")
paired = m4.groupby("patient").size()
paired = paired[paired == 2].index
m4p = m4[m4.patient.isin(paired)].copy()
c4p = c4[m4p["sample"].tolist()]
dds4, r4 = run_deseq(c4p, m4p, "~patient + day", ["day", "d28", "d0"],
                     "day28 vs day0 (paired)", "sample")
finalise(r4, "GSE270042", "day28 vs day0 (paired)")

# ---------- combine ----------
DE = pd.concat(de_tables, ignore_index=True)
GSEA = pd.concat(gsea_rows, ignore_index=True)
GSEA = GSEA.rename(columns={"Term": "module", "NES": "NES", "NOM p-val": "pval",
                            "FDR q-val": "fdr", "Gene %": "gene_pct",
                            "Lead_genes": "leading_edge", "Tag %": "tag_pct"})
keep_cols = [c for c in ["dataset", "contrast", "module", "ES", "NES", "pval", "fdr",
                         "tag_pct", "leading_edge"] if c in GSEA.columns]
GSEA = GSEA[keep_cols]

DE.to_csv("results/de_all_contrasts.csv.gz", index=False)
GSEA.to_csv("results/gsea_mitophagy_modules.csv", index=False)

MITO = DE[DE.symbol.isin(set(gs.symbol))].merge(
    gs[["symbol", "module", "in_GO_mitophagy", "in_Reactome_mitophagy"]], on="symbol", how="left")
MITO = MITO.sort_values(["contrast", "pvalue"])
MITO.to_csv("results/mitophagy_gene_results.csv", index=False)

# normalised expression for plotting
np.save("results/.placeholder.npy", np.zeros(1))
e1.to_pickle("results/expr_GSE28619.pkl"); e2.to_pickle("results/expr_GSE103580.pkl")
pd.DataFrame(dds3.layers["normed_counts"].T, index=dds3.var_names,
             columns=dds3.obs_names).to_pickle("results/expr_GSE143318.pkl")
pd.DataFrame(dds4.layers["normed_counts"].T, index=dds4.var_names,
             columns=dds4.obs_names).to_pickle("results/expr_GSE270042.pkl")
m1.to_pickle("results/meta_GSE28619.pkl"); m2.to_pickle("results/meta_GSE103580.pkl")
m3.to_pickle("results/meta_GSE143318.pkl"); m4p.to_pickle("results/meta_GSE270042.pkl")

print("contrasts:", DE.contrast.nunique(), "| genes tested per contrast:")
print(DE.groupby(["dataset", "contrast"]).size().to_string())
print("\nmitophagy genes with padj<0.05 per contrast:")
print(MITO[MITO.padj < 0.05].groupby(["dataset", "contrast"]).agg(
    n_sig=("symbol", "size"), n_up=("log2FC", lambda x: int((x > 0).sum())),
    n_down=("log2FC", lambda x: int((x < 0).sum()))).to_string())
print("\nmodule GSEA:")
print(GSEA.assign(NES=lambda d: pd.to_numeric(d.NES).round(2),
                  fdr=lambda d: pd.to_numeric(d.fdr))
      [["dataset", "contrast", "module", "NES", "pval", "fdr"]].to_string(index=False))
