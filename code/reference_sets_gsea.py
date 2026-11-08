"""Reference-set GSEA and Figure 1.

Reads results/de_all_contrasts.csv.gz (produced by de_analysis.py), fetches the two
curated GO reference sets from QuickGO, runs gseapy prerank for seven gene sets across
five contrasts, and renders figures/fig1_pathway_enrichment.png.

Reference sets are deliberately chosen to have KNOWN composition behaviour: they are
negative controls, not hypotheses.  Hepatocyte identity, fibrosis/ECM and
immune/inflammatory are explicit marker lists; OXPHOS and mitoribosome come from GO.
"""
import io, pathlib, re
import gseapy, matplotlib.pyplot as plt, numpy as np, pandas as pd, requests

ROOT = pathlib.Path(__file__).resolve().parents[1]
UNIPROT = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")

HEPATOCYTE = ["ALB","APOB","APOA1","APOA2","TTR","SERPINA1","TF","FGA","FGB","FGG",
              "CYP2E1","CYP3A4","CYP1A2","HNF4A","ASGR1","CPS1","TAT","G6PC1","PCK1"]
FIBROSIS   = ["COL1A1","COL1A2","COL3A1","ACTA2","TIMP1","LOX","FN1","PDGFRB","TGFB1","SPP1","THBS2"]
IMMUNE     = ["PTPRC","CD68","CXCL8","IL1B","CCL2","LCN2","S100A8","S100A9","CXCL6","TNF"]

CONTRASTS = [("GSE28619","AH vs Normal","AH vs\nnormal"),
             ("GSE143318","AH vs Normal","AH vs\ndonor"),
             ("GSE143318","AlcCirrhosis vs Normal","Cirrhosis\nvs donor"),
             ("GSE103580","ALD severity trend","Severity\ntrend"),
             ("GSE270042","day28 vs day0 (paired)","Day 28 vs\nday 0")]
ROWS = ["Hepatocyte identity","OXPHOS (GO:0006119)","Mitoribosome (GO:0005761)",
        "Fibrosis / ECM","Immune / inflammatory","Mitophagy (all genes)","PINK1-Parkin core"]


def quickgo_symbols(go_id):
    """Human protein-coding gene symbols annotated to a GO term and its descendants."""
    r = requests.get("https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch",
                     params={"goId": go_id, "taxonId": 9606, "geneProductType": "protein",
                             "goUsage": "descendants", "goUsageRelationships": "is_a,part_of",
                             "downloadLimit": 50000},
                     headers={"Accept": "text/tsv"}, timeout=240)
    r.raise_for_status()
    t = pd.read_csv(io.StringIO(r.text), sep="\t")
    col = [c for c in t.columns if c.strip().upper() == "SYMBOL"][0]
    # QuickGO returns UniProt accessions alongside symbols; keep symbols only.
    return sorted({s for s in t[col].dropna().astype(str) if not UNIPROT.match(s)})


def build_sets():
    mg = pd.read_csv(ROOT / "data/genesets/mitophagy_geneset.csv")
    return {"Hepatocyte identity": HEPATOCYTE,
            "OXPHOS (GO:0006119)": quickgo_symbols("GO:0006119"),
            "Mitoribosome (GO:0005761)": quickgo_symbols("GO:0005761"),
            "Fibrosis / ECM": FIBROSIS,
            "Immune / inflammatory": IMMUNE,
            "Mitophagy (all genes)": sorted(mg["symbol"].unique()),
            "PINK1-Parkin core": sorted(mg.loc[mg.module == "PINK1-Parkin core", "symbol"].unique())}


def run_gsea(de, sets):
    out = []
    for ds, ct, _ in CONTRASTS:
        sub = de[(de.dataset == ds) & (de.contrast == ct)].dropna(subset=["stat"])
        rnk = sub.groupby("symbol")["stat"].mean().sort_values(ascending=False)
        res = gseapy.prerank(rnk=rnk, gene_sets=sets, permutation_num=5000, min_size=5,
                             max_size=2000, seed=0, threads=4, no_plot=True,
                             outdir=None, verbose=False).res2d
        for _, r in res.iterrows():
            out.append({"dataset": ds, "contrast": ct, "set": r["Term"],
                        "NES": float(r["NES"]), "pval": float(r["NOM p-val"]),
                        "fdr": float(r["FDR q-val"]),
                        "n_matched": len(set(sets[r["Term"]]) & set(rnk.index))})
    return pd.DataFrame(out)


def figure1(ref, path):
    idx = ref.set_index(["dataset", "contrast", "set"])
    M = np.array([[idx.loc[(d, c, r), "NES"] for d, c, _ in CONTRASTS] for r in ROWS], float)
    Q = np.array([[idx.loc[(d, c, r), "fdr"] for d, c, _ in CONTRASTS] for r in ROWS], float)

    # The title is a claim; test it before drawing it.
    mito = M[ROWS.index("Mitophagy (all genes)")]
    bg = np.vstack([M[ROWS.index("OXPHOS (GO:0006119)")],
                    M[ROWS.index("Mitoribosome (GO:0005761)")]]).mean(0)
    assert np.all(np.sign(mito) == np.sign(bg)) and np.all(np.abs(mito) < np.abs(bg))

    fig, ax = plt.subplots(figsize=(7.4, 4.0))
    lim = np.abs(M).max()
    im = ax.imshow(M, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    for i in range(len(ROWS)):
        for j in range(len(CONTRASTS)):
            v = M[i, j]
            ax.text(j, i, f"{v:.2f}" + ("*" if Q[i, j] < 0.05 else ""), ha="center",
                    va="center", fontsize=8, color="white" if abs(v) > 0.60 * lim else "black")
    ax.set_xticks(range(len(CONTRASTS)), [c[2] for c in CONTRASTS])
    ax.set_yticks(range(len(ROWS)), ROWS)
    ax.set_xticks(np.arange(len(CONTRASTS)) - .5, minor=True)
    ax.set_yticks(np.arange(len(ROWS)) - .5, minor=True)
    ax.grid(which="minor", color="white", lw=1.4)
    ax.tick_params(which="minor", length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.axhline(4.5, color="black", lw=1.6)
    for j, (ds, _, _) in enumerate(CONTRASTS):
        ax.text(j, -0.75, ds, ha="center", va="bottom", fontsize=7, style="italic", color="0.45")
    ax.set_title("Mitophagy gene sets move with the bulk mitochondrial programme,\n"
                 "not independently of it", loc="left", pad=22)
    cb = fig.colorbar(im, ax=ax, fraction=0.030, pad=0.07)
    cb.set_label("GSEA normalised enrichment score\n(positive = up in the first-named group)", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    fig.text(0.013, 0.015, "* FDR < 0.05 (gseapy prerank, 5,000 permutations)", fontsize=6.5, color="0.45")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(path, dpi=300, bbox_inches="tight")


def main():
    de = pd.read_csv(ROOT / "results/de_all_contrasts.csv.gz")
    ref = run_gsea(de, build_sets())
    ref.to_csv(ROOT / "results/gsea_reference_sets_curatedGO.csv", index=False)
    figure1(ref, ROOT / "figures/fig1_pathway_enrichment.png")


if __name__ == "__main__":
    main()
