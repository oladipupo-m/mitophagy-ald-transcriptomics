"""Fetch the raw GEO files used by this analysis.

Nothing downloaded here is committed to the repository: provenance stays with NCBI.
Files land in data/raw/.  Roughly 200 MB total.

    python code/download_data.py
"""
import pathlib, re, sys, urllib.request

RAW = pathlib.Path(__file__).resolve().parents[1] / "data" / "raw"
SERIES = ["GSE28619", "GSE103580", "GSE143318", "GSE270042"]
PLATFORMS = ["GPL570", "GPL13667"]          # annotation tables for the two arrays
FTP = "https://ftp.ncbi.nlm.nih.gov/geo"
ACC = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"


def stem(acc):
    """GSE143318 -> GSE143nnn (GEO's directory convention)."""
    return acc[:-3] + "nnn"


def get(url, dest):
    if dest.exists():
        print(f"  have   {dest.name}")
        return
    print(f"  fetch  {dest.name}")
    urllib.request.urlretrieve(url, dest)


def suppl_filenames(gse):
    """List the supplementary files GEO holds for a series."""
    url = f"{FTP}/series/{stem(gse)}/{gse}/suppl/"
    with urllib.request.urlopen(url) as r:
        html = r.read().decode("utf-8", "replace")
    return sorted(set(re.findall(rf'href="({gse}_[^"?]+?\.(?:gz|tsv|txt|csv))"', html)))


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    for gse in SERIES:
        print(gse)
        get(f"{FTP}/series/{stem(gse)}/{gse}/matrix/{gse}_series_matrix.txt.gz",
            RAW / f"{gse}_series_matrix.txt.gz")
        for fn in suppl_filenames(gse):
            get(f"{FTP}/series/{stem(gse)}/{gse}/suppl/{fn}", RAW / fn)
    for gpl in PLATFORMS:
        print(gpl)
        get(f"{ACC}?acc={gpl}&targ=self&form=text&view=data", RAW / f"{gpl}_table.txt")
    print(f"\ndone -> {RAW}")


if __name__ == "__main__":
    sys.exit(main())
