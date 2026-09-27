#!/usr/bin/env python
"""GTEx tissue-specific gene co-expression (co-occurrence) edges.

Adapted from the hackathon's gtex_coexpr_snapshot.py; no Cosmos checkout is needed.
Requires Python 3.10+, numpy, pandas and scipy >= 1.11:
    python -m pip install numpy pandas 'scipy>=1.11'

Inputs: sample-level GTEx TPM GCT.gz files (not median-expression summaries), and
CSV/TSV with columns symbol,ensg (unversioned Ensembl gene IDs). The original
Cosmos JSON with a genes list containing human_symbol/human_id is also accepted.
GTEx input files must be obtained separately; this script does not download data.

    python gtex_coexpression.py --gtex-dir /path/to/gtex_v10 --genes genes.tsv
    python gtex_coexpression.py --gtex-dir /path/to/gtex_v10 --genes genes.tsv \
        --tissue muscle=gene_tpm_v10_muscle_skeletal.gct.gz --out results/gtex

Default tissues: muscle, subcutaneous adipose, whole blood. A gene is in the
background universe if TPM > 0.1 in at least 20% of samples and it varies across
samples. Query genes additionally require median TPM >= 1. Regress expression
PCs out of log2(TPM+1), then correlate residual ranks. For each PC count, apply
BH across all eligible query pairs and tissues. The t-based p values use
n - 2 - K degrees of freedom and are approximate, because PCs are estimated
from the same expression data. PC adjustment does not guarantee removal of all
confounding; these associations are reference evidence, not causal effects.

An edge must pass FDR < 0.05, exceed the background absolute-r 95th percentile,
and have the same sign for all PC counts (default 5,10,20). The default weight
is r at K=10. The final matrix averages supported tissue weights; opposite signs
can cancel, so inspect the per-tissue outputs and snapshot_pairs_long.csv.
Zero in an adjacency matrix means no retained edge, not proof of independence.
Outputs: pair evidence CSV, per-tissue adjacency/correlation CSVs, final adjacency
and supporting-tissue matrices, and summary.json with input and analysis settings.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

TISSUES = {"muscle": "gene_tpm_v10_muscle_skeletal.gct.gz",
           "adipose": "gene_tpm_v10_adipose_subcutaneous.gct.gz",
           "blood": "gene_tpm_v10_whole_blood.gct.gz"}


def snapshot_genes(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".json":
        g = json.loads(path.read_text())["genes"]
        df = pd.DataFrame([(x["human_symbol"], x["human_id"]) for x in g], columns=["symbol", "ensg"])
    else:
        df = pd.read_csv(path, sep="\t" if path.suffix.lower() in (".tsv", ".txt") else ",", dtype=str)
    if not {"symbol", "ensg"}.issubset(df.columns) or len(df) < 2:
        raise ValueError("Gene file requires symbol and ensg columns and at least two rows")
    df = df[["symbol", "ensg"]].copy()
    if df.isna().any().any():
        raise ValueError("Gene symbols and Ensembl IDs must not be missing")
    for col in df:
        df[col] = df[col].str.strip()
    df["ensg"] = df.ensg.str.split(".").str[0]
    if (df == "").any().any() or any(df[col].duplicated().any() for col in df):
        raise ValueError("Gene symbols and Ensembl IDs must be nonempty and unique")
    return df


def load(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", skiprows=2, index_col=0)
    df = df[~df.index.str.endswith("_PAR_Y")]
    df.index = df.index.str.split(".").str[0]
    df = df.drop(columns="Description").loc[lambda d: ~d.index.duplicated()]
    values = df.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError(f"{path}: TPM values must be finite and nonnegative")
    return df


def residualize(x: np.ndarray, pcs: np.ndarray) -> np.ndarray:
    beta, *_ = np.linalg.lstsq(pcs, x, rcond=None)
    return x - pcs @ beta


def zranks(x: np.ndarray) -> np.ndarray:
    r = np.apply_along_axis(stats.rankdata, 0, x)
    sd = r.std(0)
    return np.divide(r - r.mean(0), sd, out=np.full_like(r, np.nan), where=sd > 1e-12)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gtex-dir", type=Path, required=True)
    parser.add_argument("--genes", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("results/gtex_coexpr"))
    parser.add_argument("--tissue", action="append", metavar="NAME=FILE", help="Repeat for each TPM GCT file; paths relative to --gtex-dir")
    parser.add_argument("--pcs", type=int, nargs="+", default=[5, 10, 20])
    parser.add_argument("--main-pc", type=int, default=10)
    parser.add_argument("--background-pairs", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    KS, K_MAIN, FDR = tuple(sorted(set(args.pcs))), args.main_pc, 0.05
    if min(KS) < 0 or K_MAIN not in KS or args.background_pairs < 1:
        parser.error("PC counts must be nonnegative, --main-pc must occur in --pcs, and --background-pairs must be positive")
    tissues = dict(TISSUES)
    if args.tissue:
        tissues = {}
        for item in args.tissue:
            name, sep, filename = item.partition("=")
            if not sep or not filename or not re.fullmatch(r"[A-Za-z0-9_-]+", name) or name in tissues:
                parser.error("Each --tissue must be a unique NAME=FILE; NAME uses letters, digits, underscores or hyphens")
            tissues[name] = filename
    OUT = args.out
    genes = snapshot_genes(args.genes)
    OUT.mkdir(parents=True, exist_ok=True)
    sym = genes.symbol.tolist()
    rng = np.random.default_rng(args.seed)
    rows, summary = [], {}
    for short, tissue in tissues.items():
        tpm = load(args.gtex_dir / tissue)
        n = tpm.shape[1]
        lg = np.log2(tpm + 1)
        keep = ((tpm > 0.1).mean(axis=1) >= 0.2) & (lg.std(axis=1) > 0)
        universe = tpm.index[keep.to_numpy()]
        u = lg.loc[universe].to_numpy().T
        if len(universe) < 2 or max(KS) >= min(n - 2, len(universe)):
            raise ValueError(f"{short}: too few samples or variable genes for PC counts {KS}")
        U = np.linalg.svd((u - u.mean(0)) / u.std(0), full_matrices=False)[0]
        med = tpm.reindex(genes.ensg).median(axis=1).to_numpy()
        expressed = genes.ensg.isin(tpm.index).to_numpy() & (med >= 1) & genes.ensg.isin(universe).to_numpy()
        x = lg.reindex(genes.ensg).fillna(0).to_numpy().T
        zr_raw = zranks(x)
        r_raw = zr_raw.T @ zr_raw / n
        pairs = rng.choice(len(universe), size=(args.background_pairs, 2))
        pairs = pairs[pairs[:, 0] != pairs[:, 1]]
        if not len(pairs):
            raise ValueError("No distinct background pairs sampled; increase --background-pairs")
        per_k = {}
        for k in KS:
            pcs = np.column_stack([np.ones(n), U[:, :k]])
            zx = zranks(residualize(x, pcs))
            zu = zranks(residualize(u, pcs))
            bg = np.abs((zu[:, pairs[:, 0]] * zu[:, pairs[:, 1]]).mean(0))
            bg = bg[np.isfinite(bg)]
            if not len(bg):
                raise ValueError(f"{short}: no finite background correlations at K={k}")
            per_k[k] = (zx.T @ zx / n, float(np.quantile(bg, 0.95)), float(np.median(bg)))
        summary[short] = {"samples": n, "expressed_universe": int(len(universe)),
                          "background_abs_r_p95": {k: round(v[1], 3) for k, v in per_k.items()},
                          "median_tpm": {symbol: round(float(value), 2) if np.isfinite(value) else None for symbol, value in zip(sym, med)}}
        for i in range(len(sym)):
            for j in range(i + 1, len(sym)):
                row = {"tissue": short, "gene_a": sym[i], "gene_b": sym[j], "expressed": bool(expressed[i] and expressed[j]),
                       "r_raw": round(float(r_raw[i, j]), 3)}
                for k, (r, p95, _) in per_k.items():
                    rv = float(np.clip(r[i, j], -1, 1)); df = n - 2 - k
                    t = rv * np.sqrt(df / max(1e-12, 1 - rv * rv))
                    row[f"r_k{k}"] = rv
                    row[f"p_k{k}"] = float(2 * stats.t.sf(abs(t), df))
                    row[f"above_bg_k{k}"] = abs(rv) > p95
                rows.append(row)
    L = pd.DataFrame(rows)
    m = L.expressed
    for k in KS:
        L[f"fdr_k{k}"] = np.nan
        valid = m & np.isfinite(L[f"p_k{k}"])
        if valid.any():
            L.loc[valid, f"fdr_k{k}"] = stats.false_discovery_control(L.loc[valid, f"p_k{k}"].to_numpy(), method="bh")
    signs = np.sign(L[[f"r_k{k}" for k in KS]]).nunique(axis=1) == 1
    L["edge"] = m & signs & np.logical_and.reduce([(L[f"fdr_k{k}"] < FDR) & L[f"above_bg_k{k}"] for k in KS])
    L["weight"] = L[f"r_k{K_MAIN}"].where(L.edge, 0.0)
    L.to_csv(OUT / "snapshot_pairs_long.csv", index=False)

    def square(values: dict, fill=0.0) -> pd.DataFrame:
        M = pd.DataFrame(fill, index=sym, columns=sym, dtype=object if isinstance(fill, str) else float)
        for (ga, gb), v in values.items():
            M.loc[ga, gb] = M.loc[gb, ga] = v
        return M

    for short in tissues:
        s = L[L.tissue == short]
        square({(a, b): w for a, b, w in s[["gene_a", "gene_b", "weight"]].itertuples(index=False)}).to_csv(OUT / f"adjacency_{short}.csv")
        square({(a, b): (r if e else np.nan) for a, b, r, e in s[["gene_a", "gene_b", f"r_k{K_MAIN}", "expressed"]].itertuples(index=False)},
               fill=np.nan).to_csv(OUT / f"correlation_{short}_k{K_MAIN}.csv")
    E = L[L.edge]
    final = {(a, b): round(float(g.weight.mean()), 3) for (a, b), g in E.groupby(["gene_a", "gene_b"])}
    where = {(a, b): "+".join(g.tissue) for (a, b), g in E.groupby(["gene_a", "gene_b"])}
    square(final).to_csv(OUT / "adjacency_final.csv")
    square(where, fill="").to_csv(OUT / "adjacency_final_tissues.csv")
    (OUT / "summary.json").write_text(json.dumps({"settings": {"gtex_dir": str(args.gtex_dir), "genes": str(args.genes), "tissues": tissues, "pcs": KS, "main_pc": K_MAIN, "fdr": FDR, "background_pairs": args.background_pairs, "seed": args.seed}, "tissues": summary}, indent=2, allow_nan=False))
    print(E[["tissue", "gene_a", "gene_b", "r_raw", *[f"r_k{k}" for k in KS], f"fdr_k{K_MAIN}"]].to_string(index=False))
    print(json.dumps({t: d["background_abs_r_p95"] for t, d in summary.items()}))
    print(square(final).round(2).replace(0.0, "").to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
