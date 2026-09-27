#!/usr/bin/env python3
"""Disease and drug layer from Hetionet v1.0 [R2], added to the significant-response graph.

  python3 scripts/add_disease_layer.py --exports kg/exports_sig --hetionet raw/ppi \
          --rat processed/rat_pkg --human processed/human_pkg

Adds to kg/exports_sig/:
  nodes_Disease.csv.gz      137 Disease Ontology diseases [R19]
  nodes_Compound.csv.gz     DrugBank compounds that bind a responding gene or treat a disease [R20]
  edges_ASSOCIATED_WITH     Gene -> Disease   (Hetionet DaG: GWAS Catalog, DISEASES, DisGeNET, DOAF)     responding human genes only
  edges_UP_IN_DISEASE /     Gene -> Disease   (Hetionet DuG / DdG: STARGEO disease-vs-control expression)   responding human genes only
        DOWN_IN_DISEASE
  edges_BINDS               Compound -> Gene  (Hetionet CbG)
  edges_TREATS / PALLIATES  Compound -> Disease (Hetionet CtD / CpD)
  edges_DISEASE_GENES_ENRICHED_IN   Disease -> Tissue: disease genes over-represented among the genes that respond
                            to exercise in that species x tissue (hypergeometric test [M6], BH [M1]) [P1]
  edges_EXERCISE_OPPOSES / EXERCISE_MIMICS   Disease -> Tissue: exercise moves disease-dysregulated genes the opposite /
                            same way as the disease (exact binomial test on direction agreement, BH [M1]) [P1]
Rat genes enter through their human ortholog [D2]. Metabolites are not touched (metabolite rules unchanged).
"""
import argparse, glob, json, os
import numpy as np
import pandas as pd
from scipy import stats

ap = argparse.ArgumentParser()
ap.add_argument("--exports", default="kg/exports_sig")
ap.add_argument("--hetionet", default="raw/ppi")
ap.add_argument("--rat", default="processed/rat_pkg")
ap.add_argument("--human", default="processed/human_pkg")
a = ap.parse_args()
EX = a.exports

TISSUE_REFS = {"muscle": "D3;D4;D8", "adipose": "D3;D5;D8", "blood": "D3;D6;D8"}
EXCLUDE_OMES = {"ATAC", "METHYL", "METAB", "metab"}           # no full tested-feature universe (ATAC/METHYL) or not genes


def rd(name):
    d = pd.read_csv(f"{EX}/{name}.csv.gz", low_memory=False)
    d.columns = [c if c.startswith(":") else c.split(":")[0] for c in d.columns]
    return d


def write_nodes(label, df):
    df = df.rename(columns={"id": "id:ID"}).copy()
    df.insert(1, ":LABEL", label)
    df.to_csv(f"{EX}/nodes_{label}.csv.gz", index=False)
    return len(df)


def write_edges(etype, df, float_cols=(), int_cols=()):
    df = df.copy()
    df = df.dropna(subset=["src", "dst"]).drop_duplicates(["src", "dst"])
    # per-node ranks for the 5-edges-per-node view (by weight, within this edge type)
    df["rank_from_src"] = df.groupby("src")["weight"].rank(ascending=False, method="first").astype(int)
    df["rank_at_dst"] = df.groupby("dst")["weight"].rank(ascending=False, method="first").astype(int)
    ren = {"src": ":START_ID", "dst": ":END_ID", "weight": "weight:float", "rank_from_src": "rank_from_src:int", "rank_at_dst": "rank_at_dst:int"}
    ren.update({c: f"{c}:float" for c in float_cols}); ren.update({c: f"{c}:int" for c in int_cols})
    df = df.rename(columns=ren)
    df.insert(2, ":TYPE", etype)
    cols = [":START_ID", ":END_ID", ":TYPE"] + [c for c in df.columns if c not in (":START_ID", ":END_ID", ":TYPE")]
    df[cols].to_csv(f"{EX}/edges_{etype}.csv.gz", index=False)
    return len(df)


# ---------------------------------------------------------------- graph genes and responses
genes = rd("nodes_Gene")
hum = genes[genes.species == "human"].dropna(subset=["entrez_id"])
ez2gene = dict(zip(hum.entrez_id.astype(int), hum["id"]))
gene2ez = {v: k for k, v in ez2gene.items()}
orth = rd("edges_ORTHOLOG_OF")                                   # rat -> human
r2h = orth.groupby(":START_ID")[":END_ID"].apply(list).to_dict()
m2g = rd("edges_MAPS_TO_GENE").sort_values("weight", ascending=False).drop_duplicates(":START_ID")
f2g = dict(zip(m2g[":START_ID"], m2g[":END_ID"]))
tis = rd("nodes_Tissue")
reg = pd.concat([rd("edges_UPREGULATED_IN"), rd("edges_DOWNREGULATED_IN")])
reg = reg[reg[":END_ID"].isin(set(tis["id"]))]
reg = reg[~reg.ome.isin(EXCLUDE_OMES)]
reg["dir"] = np.where(reg[":TYPE"] == "UPREGULATED_IN", 1, -1)
reg["gene"] = reg[":START_ID"].map(f2g)
reg = reg.dropna(subset=["gene"])
rows = []
for r in reg[["species", ":END_ID", "gene", "dir"]].itertuples(index=False):
    hs = [r.gene] if r.species == "human" else r2h.get(r.gene, [])
    for h in hs:
        if h in gene2ez:
            rows.append((r.species, r[1], gene2ez[h], r.dir))
resp = pd.DataFrame(rows, columns=["species", "tissue", "ez", "dir"])
# one direction per gene x tissue (majority; ties = 0)
resp = resp.groupby(["species", "tissue", "ez"]).dir.sum().apply(np.sign).reset_index()
responding_ez = set(resp.ez)
print("responding human-ortholog genes:", len(responding_ez), "| species x tissue sets:", resp.groupby(["species", "tissue"]).ngroups)

# ---------------------------------------------------------------- universes (all tested features per species x tissue)
hf = pd.read_csv(f"{a.human}/HUMAN_FEATURE_TO_GENE.tsv.gz", sep="\t", low_memory=False, usecols=["assay", "feature_id", "entrez_gene"])
hmap = hf.dropna(subset=["entrez_gene"]).drop_duplicates(["assay", "feature_id"]).set_index(["assay", "feature_id"]).entrez_gene.astype(int)
universe = {}
for f in glob.glob(f"{a.human}/*_DA.tsv.gz"):
    if "METAB" in f or "CLINICAL" in f:
        continue
    d = pd.read_csv(f, sep="\t", usecols=["tissue", "assay", "feature_id"], low_memory=False).drop_duplicates()
    for (t, asy), s in d.groupby(["tissue", "assay"]):
        ez = hmap.reindex(list(zip([asy] * len(s), s.feature_id))).dropna().astype(int)
        universe.setdefault(("human", f"human:{t}"), set()).update(ez)
rf = pd.read_csv(f"{a.rat}/FEATURE_TO_GENE.tsv.gz", sep="\t", low_memory=False, usecols=["feature_ID", "entrez_gene"]).dropna()
rt = pd.read_csv(f"{a.rat}/RAT_TO_HUMAN_GENE.tsv.gz", sep="\t", low_memory=False).dropna(subset=["RAT_NCBI_GENE_ID", "HUMAN_ORTHOLOG_NCBI_GENE_ID"])
rat2hez = {}
for r_ez, h in zip(rt.RAT_NCBI_GENE_ID, rt.HUMAN_ORTHOLOG_NCBI_GENE_ID.astype(str)):
    for x in str(r_ez).split("|"):
        for y in h.split("|"):
            if x.strip().isdigit() and y.strip().split(".")[0].isdigit():
                rat2hez.setdefault(int(x), set()).add(int(y.strip().split(".")[0]))
rfeat = rf.groupby("feature_ID").entrez_gene.apply(lambda s: {h for x in s for h in rat2hez.get(int(float(x)), ())}).to_dict()
for f in glob.glob(f"{a.rat}/*_DA.tsv.gz"):
    if os.path.basename(f).split("_")[0] in ("METAB",):
        continue
    d = pd.read_csv(f, sep="\t", usecols=["tissue", "feature_ID"], low_memory=False).drop_duplicates()
    for t, s in d.groupby("tissue"):
        u = universe.setdefault(("rat", f"rat:{t}"), set())
        for fid in s.feature_ID:
            u.update(rfeat.get(fid, ()))
print("universes:", {k: len(v) for k, v in universe.items()})

# ---------------------------------------------------------------- Hetionet
sif = pd.read_csv(f"{a.hetionet}/hetionet_edges.sif.gz", sep="\t")
nodes = pd.read_csv(f"{a.hetionet}/hetionet_nodes.tsv", sep="\t")
name = dict(zip(nodes["id"], nodes["name"]))


def sub(meta):
    return sif[sif.metaedge == meta]


gid = lambda s: s.str.replace("Gene::", "", regex=False).astype(int)
dis = nodes[nodes.kind == "Disease"].copy()
dis["id2"] = dis["id"].str.replace("Disease::", "", regex=False)
dis_id = lambda s: s.str.replace("Disease::", "", regex=False)
cmp_id = lambda s: "DRUGBANK:" + s.str.replace("Compound::", "", regex=False)

DaG = sub("DaG").assign(dz=lambda d: dis_id(d.source), ez=lambda d: gid(d.target))
DuG = sub("DuG").assign(dz=lambda d: dis_id(d.source), ez=lambda d: gid(d.target), ddir=1)
DdG = sub("DdG").assign(dz=lambda d: dis_id(d.source), ez=lambda d: gid(d.target), ddir=-1)
CbG = sub("CbG").assign(cp=lambda d: cmp_id(d.source), ez=lambda d: gid(d.target))
CtD = sub("CtD").assign(cp=lambda d: cmp_id(d.source), dz=lambda d: dis_id(d.target))
CpD = sub("CpD").assign(cp=lambda d: cmp_id(d.source), dz=lambda d: dis_id(d.target))

# ---------------------------------------------------------------- nodes
nD = write_nodes("Disease", pd.DataFrame({"id": dis.id2, "name": dis.name, "doid": dis.id2,
                                          "url": "https://disease-ontology.org/?id=" + dis.id2, "source_refs": "R2;R19"}))
cb = CbG[CbG.ez.isin(responding_ez)]
keep_c = set(cb.cp) | set(CtD.cp) | set(CpD.cp)
cn = nodes[nodes.kind == "Compound"].assign(cid=lambda d: cmp_id(d["id"]))
cn = cn[cn.cid.isin(keep_c)]
db = cn.cid.str.replace("DRUGBANK:", "", regex=False)
nC = write_nodes("Compound", pd.DataFrame({"id": cn.cid, "name": cn.name, "drugbank_id": db,
                                           "url": "https://go.drugbank.com/drugs/" + db, "source_refs": "R2;R20"}))

# ---------------------------------------------------------------- gene-level edges (responding human genes)
counts = {}
d = DaG[DaG.ez.isin(responding_ez)]
counts["ASSOCIATED_WITH"] = write_edges("ASSOCIATED_WITH", pd.DataFrame({"src": d.ez.map(ez2gene), "dst": d.dz, "weight": 1.0,
    "evidence": "Hetionet DaG (GWAS Catalog, DISEASES, DisGeNET, DOAF)", "significance_basis": "curated association", "source_refs": "R2;R19"}))
for et, dd in [("UP_IN_DISEASE", DuG), ("DOWN_IN_DISEASE", DdG)]:
    d = dd[dd.ez.isin(responding_ez)]
    counts[et] = write_edges(et, pd.DataFrame({"src": d.ez.map(ez2gene), "dst": d.dz, "weight": 1.0,
        "evidence": "Hetionet D%sG (STARGEO disease vs control expression)" % ("u" if et.startswith("UP") else "d"),
        "significance_basis": "Hetionet signature", "source_refs": "R2;R19"}))
counts["BINDS"] = write_edges("BINDS", pd.DataFrame({"src": cb.cp, "dst": cb.ez.map(ez2gene), "weight": 1.0,
    "evidence": "Hetionet CbG (DrugBank, ChEMBL, BindingDB, DrugCentral)", "significance_basis": "curated binding", "source_refs": "R2;R20"}))
for et, dd, ev in [("TREATS", CtD, "Hetionet CtD (indications, expert-curated)"), ("PALLIATES", CpD, "Hetionet CpD (symptomatic indications)")]:
    counts[et] = write_edges(et, pd.DataFrame({"src": dd.cp, "dst": dd.dz, "weight": 1.0, "evidence": ev,
                                               "significance_basis": "curated indication", "source_refs": "R2;R19;R20"}))

# ---------------------------------------------------------------- disease enrichment among responders (hypergeometric)
dz_genes = DaG.groupby("dz").ez.apply(set).to_dict()
sig = {}
for (sp, t), g in resp.groupby(["species", "tissue"]):
    sig[(sp, t)] = dict(zip(g.ez, g.dir))
tests = []
for (sp, t), R in sig.items():
    U = universe.get((sp, t))
    if not U:
        continue
    Rset = set(R) & U
    N, n = len(U), len(Rset)
    for dz, G in dz_genes.items():
        GU = G & U
        K = len(GU)
        hits = GU & Rset
        k = len(hits)
        if K < 5 or k < 3:
            continue
        p = stats.hypergeom.sf(k - 1, N, K, n)
        up = sum(1 for e in hits if R[e] > 0); down = sum(1 for e in hits if R[e] < 0)
        tests.append((dz, t, sp, k, K, n, N, (k / K) / (n / N), p, up, down,
                      ";".join(sorted(genes.set_index("id").symbol.reindex([ez2gene.get(e) for e in hits]).dropna().astype(str))[:40])))
te = pd.DataFrame(tests, columns=["dz", "tissue_id", "species", "k", "K", "n", "N", "fold", "p", "n_up", "n_down", "genes"])
te["adj_p"] = stats.false_discovery_control(te.p) if len(te) else []
te = te[(te.adj_p < 0.05) & (te.fold > 1)].copy()
te["weight"] = (-np.log10(te.adj_p.clip(lower=1e-300))).rank(pct=True).round(4)
te["tissue"] = te.tissue_id.str.split(":").str[1]
te["refs"] = np.where(te.species == "rat", "D1;D2", te.tissue.map(TISSUE_REFS)) + ";R2;R19;M6;M1;P1"
counts["DISEASE_GENES_ENRICHED_IN"] = write_edges("DISEASE_GENES_ENRICHED_IN", pd.DataFrame({
    "src": te.dz, "dst": te.tissue_id, "species": te.species, "tissue": te.tissue, "overlap": te.k, "disease_genes_measured": te.K,
    "responding_genes": te.n, "measured_genes": te.N, "fold_enrichment": te.fold.round(3), "p_value": te.p, "adj_p_value": te.adj_p,
    "n_up": te.n_up, "n_down": te.n_down, "genes": te.genes, "weight": te.weight,
    "significance_basis": "hypergeometric, BH adj_p<0.05 across all disease x tissue tests", "source_refs": te.refs}),
    float_cols=["fold_enrichment", "p_value", "adj_p_value"], int_cols=["overlap", "disease_genes_measured", "responding_genes", "measured_genes", "n_up", "n_down"])

# ---------------------------------------------------------------- direction: does exercise oppose or mimic the disease signature?
sigd = pd.concat([DuG, DdG]).groupby(["dz", "ez"]).ddir.sum().apply(np.sign)
sigd = sigd[sigd != 0]
dsig = {dz: s.droplevel(0).to_dict() for dz, s in sigd.groupby(level=0)}
rows = []
for (sp, t), R in sig.items():
    for dz, S in dsig.items():
        both = [(S[e], R[e]) for e in set(S) & set(R) if R[e] != 0]
        n = len(both)
        if n < 10:
            continue
        opp = sum(1 for s, r in both if s != r)
        p = stats.binomtest(opp, n, 0.5).pvalue
        rows.append((dz, t, sp, n, opp, n - opp, opp / n, p))
co = pd.DataFrame(rows, columns=["dz", "tissue_id", "species", "n", "opposite", "same", "frac_opposite", "p"])
co["adj_p"] = stats.false_discovery_control(co.p) if len(co) else []
co = co[co.adj_p < 0.05].copy()
co["tissue"] = co.tissue_id.str.split(":").str[1]
co["refs"] = np.where(co.species == "rat", "D1;D2", co.tissue.map(TISSUE_REFS)) + ";R2;R19;M1;P1"
for et, mask, w in [("EXERCISE_OPPOSES", co.frac_opposite > 0.5, co.frac_opposite), ("EXERCISE_MIMICS", co.frac_opposite < 0.5, 1 - co.frac_opposite)]:
    s = co[mask]
    counts[et] = write_edges(et, pd.DataFrame({"src": s.dz, "dst": s.tissue_id, "species": s.species, "tissue": s.tissue,
        "genes_compared": s.n, "opposite": s.opposite, "same": s.same, "fraction_opposite": s.frac_opposite.round(4),
        "p_value": s.p, "adj_p_value": s.adj_p, "weight": w[mask].round(4),
        "significance_basis": "exact binomial test on direction agreement vs 0.5, BH adj_p<0.05", "source_refs": s.refs}),
        float_cols=["fraction_opposite", "p_value", "adj_p_value"], int_cols=["genes_compared", "opposite", "same"])

# ---------------------------------------------------------------- manifest
m = json.load(open(f"{EX}/manifest.json"))
m["nodes"].update({"Disease": nD, "Compound": nC})
for k, v in counts.items():
    m["edges"][k] = {"n": v, "dropped_dangling_or_duplicate": 0}
m["disease_layer"] = {"source": "Hetionet v1.0", "tests_enrichment": int(len(tests)), "tests_direction": int(len(rows))}
json.dump(m, open(f"{EX}/manifest.json", "w"), indent=2)
print("nodes: Disease", nD, "Compound", nC); print("edges:", counts)
