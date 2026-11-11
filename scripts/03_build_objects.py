"""Build per-library expression objects with QC metrics, clone IDs and cluster labels.

CLAUDE.md "Next data tasks" item 3. Runs in the cloud workspace (Linux, scanpy).

- One AnnData per 10x library (22), raw integer counts, cells x 33,694 genes.
  EGFP (lineage-vector reporter, present only in replicates 2-3) is removed from the
  gene matrix and kept as obs['egfp_counts'].
- Replicate ID comes from the viral barcode via results/02_library_map (GEO sample
  names swap replicates 2/3 in Exp3).
- Clone-table cells are joined on (library, core barcode), so stage-IV barcode
  collisions across libraries cannot mis-assign cells.
- QC uses the paper's published thresholds (You et al. 2023, STAR Methods), not tuned:
  >= 1,000 genes, <= 30,000 UMIs, <= 15% mitochondrial UMIs. The paper additionally
  dropped clusters not reproducible between replicates; that step is not reproduced.

Outputs
  data/processed/libs/<library>.h5ad      (git-ignored; large)
  results/03_objects/cell_obs.csv.gz      every barcode: library, stage, exp, rep, QC,
                                          qc_pass, clone, Cluster.Id
  results/03_objects/qc_by_library.csv
  results/03_objects/clone_counts_after_qc.csv
  results/03_objects/summary.csv
  results/03_objects/provenance.txt
"""
import hashlib
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import scipy.io
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
MTX = ROOT / "data/interim/10x"
MAP = ROOT / "results/02_library_map"
LIBOUT = ROOT / "data/processed/libs"
OUT = ROOT / "results/03_objects"
TARS = [ROOT / "data/raw/GSE207921/GSE207921_RAW.tar",
        ROOT / "data/raw/GSE221592/GSE221592_RAW.tar"]
MIN_GENES, MAX_UMI, MAX_MT = 1000, 30000, 15.0
PUBLISHED = {"I->II": (768, 2762), "II->III": (1213, 5856), "III->IV": (559, 2816)}
PAIRS = {"I->II": ("I", "II", 1), "II->III": ("II", "III", 2), "III->IV": ("III", "IV", 3)}


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def load_library(lib):
    feat = next(MTX.glob(f"{lib}_features.tsv.gz"), None) or next(MTX.glob(f"{lib}_genes.tsv.gz"))
    genes = pd.read_csv(feat, sep="\t", header=None).iloc[:, :2]
    genes.columns = ["gene_id", "symbol"]
    bc = pd.read_csv(MTX / f"{lib}_barcodes.tsv.gz", header=None)[0]
    X = sp.csr_matrix(scipy.io.mmread(MTX / f"{lib}_matrix.mtx.gz").T, dtype=np.float32)
    is_egfp = (genes.gene_id == "EGFP").to_numpy()
    egfp = np.asarray(X[:, is_egfp].sum(axis=1)).ravel() if is_egfp.any() else np.zeros(X.shape[0])
    X = X[:, ~is_egfp]
    genes = genes[~is_egfp].reset_index(drop=True)
    a = ad.AnnData(X=X, var=pd.DataFrame(index=genes.gene_id.values))
    a.var["symbol"] = genes.symbol.values
    a.obs_names = (lib + ":" + bc.str.split("-").str[0]).values
    a.obs["core"] = bc.str.split("-").str[0].values
    a.obs["egfp_counts"] = egfp
    return a


def main():
    libmap = pd.read_csv(MAP / "library_coverage.csv")
    lib_meta = libmap.groupby("library").agg(stage=("stage", "first"), exp=("exp", "first"),
                                             rep=("rep", "first"), gsm=("gsm", "first"))
    assert lib_meta.groupby("library").size().eq(1).all()
    cells = pd.read_csv(MAP / "cell_library_map.csv.gz")
    cells = cells.rename(columns={"Cluster.Id": "cluster"})
    LIBOUT.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)

    obs_all, gene_ids = [], None
    for lib, m in lib_meta.iterrows():
        a = load_library(lib)
        if gene_ids is None:
            gene_ids = a.var_names.copy()
        assert a.var_names.equals(gene_ids), f"gene order differs in {lib}"
        a.var["mt"] = a.var.symbol.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
        a.obs["library"], a.obs["gsm"] = lib, m.gsm
        a.obs["stage"], a.obs["exp"], a.obs["rep"] = m.stage, int(m.exp), int(m.rep)
        a.obs["qc_pass"] = ((a.obs.n_genes_by_counts >= MIN_GENES)
                            & (a.obs.total_counts <= MAX_UMI)
                            & (a.obs.pct_counts_mt <= MAX_MT))
        c = cells[cells.library == lib].set_index("core")[["clone", "cluster", "Cell.Barcode"]]
        assert c.index.is_unique, f"duplicate clone-table barcodes in {lib}"
        a.obs = a.obs.join(c, on="core")
        a.obs["in_clone_table"] = a.obs.clone.notna()
        a.X = a.X.astype(np.int32) if np.allclose(a.X.data, np.round(a.X.data)) else a.X
        a.write_h5ad(LIBOUT / f"{lib}.h5ad", compression="gzip")
        obs_all.append(a.obs.reset_index(names="cell_id"))
        print(f"{lib}: {a.n_obs} barcodes, {int(a.obs.qc_pass.sum())} pass, "
              f"{int(a.obs.in_clone_table.sum())} in clone table", flush=True)
        del a

    obs = pd.concat(obs_all, ignore_index=True)
    keep = ["cell_id", "library", "gsm", "stage", "exp", "rep", "core", "n_genes_by_counts",
            "total_counts", "pct_counts_mt", "egfp_counts", "qc_pass", "in_clone_table",
            "clone", "cluster", "Cell.Barcode"]
    obs[keep].to_csv(OUT / "cell_obs.csv.gz", index=False)

    qc = obs.groupby(["library", "stage", "exp", "rep"]).agg(
        n_barcodes=("cell_id", "size"), n_qc_pass=("qc_pass", "sum"),
        n_clone_table=("in_clone_table", "sum"),
        n_clone_table_qc_pass=("in_clone_table", lambda s: int((s & obs.loc[s.index, "qc_pass"]).sum())),
        median_genes=("n_genes_by_counts", "median"), median_umi=("total_counts", "median"),
        median_pct_mt=("pct_counts_mt", "median")).reset_index()
    qc.to_csv(OUT / "qc_by_library.csv", index=False)

    rows = []
    for name, (s0, s1, e) in PAIRS.items():
        for label, sub in [("all_clone_table_cells", obs[obs.in_clone_table]),
                           ("qc_pass_only", obs[obs.in_clone_table & obs.qc_pass])]:
            d = sub[(sub.exp == e) & sub.stage.isin([s0, s1])]
            st = d.groupby("clone")["stage"].agg(set)
            shared = st[st.apply(lambda x: {s0, s1} <= x)].index
            rows.append(dict(transition=name, cells_used=label, clones=len(shared),
                             cells=int(d.clone.isin(shared).sum()),
                             published_clones=PUBLISHED[name][0],
                             published_cells=PUBLISHED[name][1]))
    cc = pd.DataFrame(rows)
    cc.to_csv(OUT / "clone_counts_after_qc.csv", index=False)

    summary = pd.Series(dict(
        libraries=obs.library.nunique(), genes=len(gene_ids), barcodes=len(obs),
        qc_pass=int(obs.qc_pass.sum()), published_hq_cells=173564,
        clone_table_cells_matched=int(obs.in_clone_table.sum()),
        clone_table_cells_total=len(cells),
        clone_table_cells_failing_paper_qc=int((obs.in_clone_table & ~obs.qc_pass).sum()),
        thresholds=f"genes>={MIN_GENES}; UMI<={MAX_UMI}; mt%<={MAX_MT}"))
    summary.to_csv(OUT / "summary.csv", header=["value"])

    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    prov = [f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
            f"python {sys.version.split()[0]} scanpy {sc.__version__} anndata {ad.__version__} "
            f"pandas {pd.__version__} numpy {np.__version__} scipy {scipy.__version__}",
            "environment: Claude cloud workspace (Linux)"]
    prov += [f"input {t.name} sha256 {sha256(t)}" for t in TARS]
    prov += [f"input {p.relative_to(ROOT)} sha256 {sha256(p)}"
             for p in [MAP / "cell_library_map.csv.gz", MAP / "library_coverage.csv"]]
    (OUT / "provenance.txt").write_text("\n".join(prov) + "\n")

    pd.set_option("display.width", 200)
    print(summary.to_string())
    print(qc.to_string(index=False))
    print(cc.to_string(index=False))


if __name__ == "__main__":
    main()
