#!/usr/bin/env python3
"""Adjacency matrices from the significant-response graph, in the hackathon repo's format
(rows x columns CSV with a row-label column, plus row/column metadata and a long edge list).

  python3 scripts/build_adjacency.py kg/exports_sig kg/adjacency

Cell value for response matrices = sign(direction) x -log10(adj p)  (0 = not significant).
  human adj p = BH within the control-adjusted contrast [D8, M1]; rat adj p = MoTrPAC training q [D1].
When several molecules of one gene fall in the same cell, the one with the smallest adj p is kept.
"""
import os, sys
import numpy as np
import pandas as pd

EX = sys.argv[1] if len(sys.argv) > 1 else "kg/exports_sig"
OUT = sys.argv[2] if len(sys.argv) > 2 else "kg/adjacency"
os.makedirs(OUT, exist_ok=True)


def rd(name):
    d = pd.read_csv(f"{EX}/{name}.csv.gz", low_memory=False)
    d.columns = [c if c.startswith(":") else c.split(":")[0] for c in d.columns]
    return d


GROUP_CODE = {"human:endurance_acute": "EE-CON", "human:resistance_acute": "RE-CON", "rat:endurance_training": "TRN-SED"}
HUMAN_TP = ["during_20_min", "during_40_min", "post_10_min", "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr"]
RAT_TP = ["1w", "2w", "4w", "8w"]
LAYER = {"transcript-rna-seq": "transcript", "prot-pr": "protein", "prot-ph": "phosphosite", "prot-ac": "acetylsite",
         "prot-ub-protein-corrected": "ubiquitylsite", "prot-ol": "olink", "immunoassay": "immunoassay",
         "epigen-atac-seq": "atac", "epigen-rrbs": "methylation"}

# ---------------------------------------------------------------- load
reg = pd.concat([rd("edges_UPREGULATED_IN"), rd("edges_DOWNREGULATED_IN")])
tis = rd("nodes_Tissue")
reg = reg[reg[":END_ID"].isin(set(tis["id"]))].copy()          # organ edges carry group/timepoint as properties
reg["dir"] = np.where(reg[":TYPE"] == "UPREGULATED_IN", 1, -1)
reg["nl"] = -np.log10(reg["adj_p_value"].astype(float).clip(lower=1e-300))
reg["value"] = (reg["dir"] * reg["nl"]).round(4)
reg["gcode"] = reg["group"].map(GROUP_CODE)
reg["layer"] = reg["assay"].map(LAYER).fillna("metabolite")
reg["is_metab"] = reg[":START_ID"].str.startswith("REFMET:") | reg["ome"].str.lower().str.startswith("metab")

genes = rd("nodes_Gene").set_index("id")
metab = rd("nodes_Metabolite").set_index("id")
m2g = rd("edges_MAPS_TO_GENE").sort_values(["weight", ":END_ID"], ascending=[False, True])
best_gene = m2g.drop_duplicates(":START_ID").set_index(":START_ID")[":END_ID"]          # one gene per molecule
orth = rd("edges_ORTHOLOG_OF")
rat2hum = orth.groupby(":START_ID")[":END_ID"].apply(lambda s: ";".join(genes.symbol.reindex(s).fillna("").astype(str))).to_dict()


def metab_name(fid):
    n = metab.refmet_name.get(fid)
    return n if isinstance(n, str) and n else str(metab.source_id.get(fid, fid))


def col_label(d, species):
    if species == "human":
        return d["tissue"] + "|" + d["gcode"] + "|" + d["timepoint"]
    return d["tissue"] + "|" + d["sex"].str[0].str.upper() + "|" + d["timepoint"]


def col_order(species, tissues):
    if species == "human":
        cols = []
        for t in tissues:
            for g, tps in [("EE-CON", HUMAN_TP), ("RE-CON", HUMAN_TP[2:])]:
                cols += [f"{t}|{g}|{tp}" for tp in tps]
        return cols
    return [f"{t}|{s}|{w}" for t in tissues for s in "FM" for w in RAT_TP]


def matrix(long, row_key, col, species, row_order):
    """long: one row per (row_key, col) after picking the smallest adj p; returns wide matrix."""
    w = long.pivot(index=row_key, columns=col, values="value")
    cols = [c for c in col_order(species, sorted(long.tissue.unique())) if c in w.columns]
    return w.reindex(index=row_order, columns=cols).fillna(0.0)


def pick_best(d, keys):
    d = d.sort_values(["nl"], ascending=False)
    agg = d.groupby(keys).agg(n_molecules=(":START_ID", "nunique"), n_up=("dir", lambda s: int((s > 0).sum())),
                              n_down=("dir", lambda s: int((s < 0).sum())))
    top = d.drop_duplicates(keys).set_index(keys)
    return top.join(agg).reset_index()


def write(name, M, rows_meta, cols_meta, edges, note):
    R, C = M.shape
    base = f"{name}_{R}x{C}"
    M.index.name = "row"; M.columns.name = None
    M.to_csv(f"{OUT}/{base}.csv", float_format="%.4g")
    rm = rows_meta.reindex(M.index); rm.index.name = "row"; rm.to_csv(f"{OUT}/{base}_rows.csv")
    cm = cols_meta.reindex(M.columns); cm.index.name = "column"; cm.to_csv(f"{OUT}/{base}_cols.csv")
    edges.to_csv(f"{OUT}/{base}_edges.csv", index=False, float_format="%.4g")
    nz = int((M.values != 0).sum())
    print(f"{base}: {nz} non-zero cells ({nz / M.size:.2%}); {note}")
    return dict(file=f"{base}.csv", rows=R, cols=C, nonzero=nz, density=nz / M.size, note=note)


def cols_meta_for(cols, species):
    parts = pd.Series(cols).str.split("|", expand=True)
    if species == "human":
        cm = pd.DataFrame({"column": cols, "tissue": parts[0].values, "contrast": parts[1].values, "timepoint": parts[2].values,
                           "group": parts[1].map({"EE-CON": "acute endurance vs control", "RE-CON": "acute resistance vs control"}).values})
    else:
        cm = pd.DataFrame({"column": cols, "tissue": parts[0].values, "sex": parts[1].map({"F": "female", "M": "male"}).values,
                           "week": parts[2].values, "contrast": "8-wk endurance training vs sedentary control"})
    return cm.set_index("column")


def dedupe(label):
    """make row/column labels unique by appending the ID tail to duplicated names"""
    label = pd.Series([str(x) for x in label.values], index=label.index, dtype=object)
    dup = label.duplicated(keep=False)
    tail = pd.Series([str(i).split(":")[-1] for i in label.index], index=label.index, dtype=object)
    label[dup] = label[dup] + "|" + tail[dup]
    return label


summary = []

# ---------------------------------------------------------------- 1-2. gene x contrast (per species)
for sp in ["human", "rat"]:
    d = reg[(reg.species == sp) & ~reg.is_metab].copy()
    d["gene"] = d[":START_ID"].map(best_gene)
    d = d.dropna(subset=["gene"])
    d["col"] = col_label(d, sp)
    best = pick_best(d, ["gene", "col"])
    order = best.groupby("gene").nl.max().sort_values(ascending=False).index          # strongest genes first
    sym = genes.symbol.reindex(order).fillna("").astype(str)
    label = np.where(sym != "", sym, order)
    label = pd.Series(label, index=order)
    label = dedupe(label)                                                               # keep symbols unique
    long = best.assign(row=best.gene.map(label))
    M = matrix(long, "row", "col", sp, list(label.values))
    rm = pd.DataFrame({"row": label.values, "gene_id": order, "symbol": sym.values,
                       "entrez_id": genes.entrez_id.reindex(order).values,
                       "human_ortholog": [rat2hum.get(g, "") for g in order] if sp == "rat" else "",
                       "n_cells": long.groupby("gene").col.nunique().reindex(order).values,
                       "layers": d.groupby("gene").layer.apply(lambda s: ";".join(sorted(set(s)))).reindex(order).values,
                       "best_adj_p": d.groupby("gene").adj_p_value.min().reindex(order).values}).set_index("row")
    ed = long[["row", "gene", "col", "value", "adj_p_value", "dir", "layer", ":START_ID", "n_molecules", "n_up", "n_down",
               "z", "significance_basis", "source_refs"]].rename(columns={":START_ID": "best_molecule", "col": "column", "dir": "direction"})
    summary.append(write(f"{sp}_gene_contrast_adjacency", M, rm, cols_meta_for(M.columns, sp), ed,
                         "genes x (tissue|contrast|timepoint), signed -log10 adj p of the gene's strongest molecule"))

# ---------------------------------------------------------------- 3-4. metabolite x contrast (per species)
for sp in ["human", "rat"]:
    d = reg[(reg.species == sp) & reg.is_metab].copy()
    d["col"] = col_label(d, sp)
    best = pick_best(d, [":START_ID", "col"])
    order = best.groupby(":START_ID").nl.max().sort_values(ascending=False).index
    label = pd.Series([metab_name(f) for f in order], index=order)
    label = dedupe(label)
    long = best.assign(row=best[":START_ID"].map(label))
    M = matrix(long, "row", "col", sp, list(label.values))
    rm = pd.DataFrame({"row": label.values, "metabolite_id": order,
                       "refmet_name": metab.refmet_name.reindex(order).values, "platform": metab.platform.reindex(order).values,
                       "shared_with_other_species": metab.species.reindex(order).eq("rat;human").values,
                       "n_cells": long.groupby(":START_ID").col.nunique().reindex(order).values,
                       "best_adj_p": d.groupby(":START_ID").adj_p_value.min().reindex(order).values}).set_index("row")
    ed = long[["row", ":START_ID", "col", "value", "adj_p_value", "dir", "assay", "z", "significance_basis", "source_refs"]] \
        .rename(columns={":START_ID": "metabolite_id", "col": "column", "dir": "direction"})
    summary.append(write(f"{sp}_metab_contrast_adjacency", M, rm, cols_meta_for(M.columns, sp), ed,
                         "metabolites x (tissue|contrast|timepoint), signed -log10 adj p"))

# ---------------------------------------------------------------- 5. gene x tissue (both species, human-ortholog rows)
d = reg[~reg.is_metab].copy()
d["gene"] = d[":START_ID"].map(best_gene)
d = d.dropna(subset=["gene"])
hum_of = orth.drop_duplicates(":START_ID").set_index(":START_ID")[":END_ID"]   # first human ortholog
d["hgene"] = np.where(d.species == "human", d.gene, d.gene.map(hum_of))
d = d.dropna(subset=["hgene"])
d["col"] = d.species + ":" + d.tissue
best = pick_best(d, ["hgene", "col"])
order = best.groupby("hgene").nl.max().sort_values(ascending=False).index
label = pd.Series(genes.symbol.reindex(order).fillna("").astype(str).values, index=order)
label = pd.Series(np.where(label.values == "", label.index, label.values), index=label.index)
label = dedupe(label)
long = best.assign(row=best.hgene.map(label))
tcols = [f"human:{t}" for t in ["muscle", "adipose", "blood"]] + sorted(c for c in long.col.unique() if c.startswith("rat:"))
M = long.pivot(index="row", columns="col", values="value").reindex(index=list(label.values), columns=[c for c in tcols if c in set(long.col)]).fillna(0.0)
rm = pd.DataFrame({"row": label.values, "human_gene_id": order, "symbol": genes.symbol.reindex(order).values,
                   "responds_in_human": d[d.species == "human"].groupby("hgene").size().reindex(order).fillna(0).gt(0).values,
                   "responds_in_rat": d[d.species == "rat"].groupby("hgene").size().reindex(order).fillna(0).gt(0).values,
                   "n_tissues": long.groupby("hgene").col.nunique().reindex(order).values}).set_index("row")
cm = pd.DataFrame({"column": M.columns, "species": [c.split(":")[0] for c in M.columns], "tissue": [c.split(":")[1] for c in M.columns]}).set_index("column")
ed = long[["row", "hgene", "col", "value", "adj_p_value", "dir", "species", "layer", ":START_ID", "n_molecules", "n_up", "n_down", "source_refs"]] \
    .rename(columns={"hgene": "human_gene_id", ":START_ID": "best_molecule", "col": "column", "dir": "direction"})
summary.append(write("gene_tissue_adjacency", M, rm, cm, ed,
                     "human-ortholog genes x species:tissue, signed -log10 adj p of the strongest response at any timepoint"))

# ---------------------------------------------------------------- 6-7. metabolite x protein co-regulation (per species)
co = rd("edges_CO_REGULATED_WITH")
co = co[co.partner_type != "metabolite"] if "partner_type" in co else co
feat_names = {}
for lab in ["ProteinFeature", "AffinityProtein"]:
    n = rd(f"nodes_{lab}")
    for fid in n["id"]:
        g = best_gene.get(fid)
        s = genes.symbol.get(g) if isinstance(g, str) else None
        feat_names[fid] = f"{s}" if isinstance(s, str) and s else fid.split(":")[-1]
for sp in ["human", "rat"]:
    c = co[co.species == sp].copy()
    ism = c[":START_ID"].str.startswith("REFMET:") | c[":START_ID"].isin(metab.index)
    c["met"] = np.where(ism, c[":START_ID"], c[":END_ID"])
    c["prot"] = np.where(ism, c[":END_ID"], c[":START_ID"])
    c = c[c.met.isin(metab.index) & ~c.prot.isin(metab.index)]
    if c.empty:
        continue
    c["abs"] = c.pearson_r.abs()
    c = c.sort_values("abs", ascending=False).drop_duplicates(["met", "prot"])   # strongest tissue per pair
    ml = pd.Series({m: metab_name(m) for m in c.met.unique()})
    pl = pd.Series({p: feat_names.get(p, p.split(":")[-1]) + "|" + p.split(":")[1] for p in c.prot.unique()})
    ml, pl = dedupe(ml), dedupe(pl)
    c["row"] = c.met.map(ml); c["colname"] = c.prot.map(pl)
    rows = c.groupby("row")["abs"].max().sort_values(ascending=False).index
    cols = c.groupby("colname")["abs"].max().sort_values(ascending=False).index
    M = c.pivot(index="row", columns="colname", values="pearson_r").reindex(index=rows, columns=cols).fillna(0.0).round(3)
    rm = pd.DataFrame({"row": rows}).set_index("row").join(c.drop_duplicates("row").set_index("row")[["met"]].rename(columns={"met": "metabolite_id"}))
    rm["n_proteins"] = c.groupby("row").size().reindex(rows).values
    cm = pd.DataFrame({"column": cols}).set_index("column").join(c.drop_duplicates("colname").set_index("colname")[["prot"]].rename(columns={"prot": "protein_feature_id"}))
    cm["gene"] = [best_gene.get(p, "") for p in cm.protein_feature_id]
    ed = c[["row", "met", "colname", "prot", "pearson_r", "tissue", "n_points", "sign", "significance_basis", "source_refs"]] \
        .rename(columns={"met": "metabolite_id", "colname": "column", "prot": "protein_feature_id"})
    summary.append(write(f"{sp}_metab_protein_coreg_adjacency", M, rm, cm, ed,
                         "metabolites x proteins, Pearson r of response profiles (|r| >= 0.9, top 10 per metabolite; strongest tissue kept)"))

pd.DataFrame(summary).to_csv(f"{OUT}/_index.csv", index=False)
print("wrote", OUT)
