"""A3.7 (ANALYSIS_PLAN amendment A3, post-review): RQ1 restricted to single-barcode clones.

Recomputes the pooled forward RQ1 summary (clone-weighted means, 1,000-resample clone
bootstrap, seed 0, as in 06_rq1_coupling.summarise) from the saved per-cell metrics,
for (a) all cross-stage clones and (b) clones whose viral-barcode set is one barcode
(clone string without ';'). Primary cell set, capped pools, couplings: OT balanced
eps = 0.05, kNN, independent. Nothing is refitted.
"""
import hashlib
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INP = ROOT / "results/06_rq1/cells_forward.csv.gz"
OUT = ROOT / "results/11_single_barcode"
OUT.mkdir(parents=True, exist_ok=True)
COUPLINGS = ["OT_bal_eps0.05", "kNN", "independent"]
COLS = ["auroc", "enrichment", "m", "m0", "top1", "top10"]
N_BOOT = 1000


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    f = pd.read_csv(INP, low_memory=False)
    f = f[(f.cellset == "primary") & (~f.uncapped.astype(bool)) & f.coupling.isin(COUPLINGS)].copy()
    f["multi"] = f.clone.astype(str).str.contains(";")
    f["cl"] = f.rep.astype(str) + "|" + f.clone.astype(str)
    f["enrichment"] = f.m / f.m0
    rng = np.random.default_rng(0)
    rows = []
    for (t, cp), g in f.groupby(["transition", "coupling"]):
        for subset, h in (("all", g), ("single_barcode", g[~g.multi])):
            c = h.groupby("cl")[COLS].mean()
            arr = c.to_numpy()
            idx = rng.integers(0, len(c), size=(N_BOOT, len(c)))
            boots = np.stack([arr[i].mean(0) for i in idx])
            r = dict(transition=t, coupling=cp, subset=subset, n_clones=len(c), n_cells=len(h),
                     pct_cells_multi=100 * g.multi.mean())
            for k, col in enumerate(COLS):
                r[col] = arr[:, k].mean()
                r[col + "_lo"], r[col + "_hi"] = np.percentile(boots[:, k], [2.5, 97.5])
            rows.append(r)
    pd.DataFrame(rows).to_csv(OUT / "summary.csv", index=False)
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    (OUT / "provenance.txt").write_text("\n".join([
        f"git_commit {commit}{' (DIRTY scripts/)' if dirty else ''}",
        "plan: ANALYSIS_PLAN.md amendment A3.7 (post-review; preliminary check disclosed there)",
        f"python {sys.version.split()[0]} pandas {pd.__version__} numpy {np.__version__}",
        f"input {INP.relative_to(ROOT).as_posix()} sha256 {sha256(INP)}"]) + "\n")
    print(pd.DataFrame(rows)[["transition", "coupling", "subset", "n_clones", "auroc", "enrichment"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
