#!/usr/bin/env python3
"""Compact columnar JSON of kg/exports_sig for the web explorer (GO omitted to keep it light)."""
import json, sys, os
import numpy as np
import pandas as pd

EX = sys.argv[1] if len(sys.argv) > 1 else "kg/exports_sig"
OUT = sys.argv[2] if len(sys.argv) > 2 else "explorer/data.json"


def N(label):
    d = pd.read_csv(f"{EX}/nodes_{label}.csv.gz", low_memory=False)
    d.columns = [c.split(":")[0] if not c.startswith(":") else c for c in d.columns]
    return d


def E(t):
    d = pd.read_csv(f"{EX}/edges_{t}.csv.gz", low_memory=False)
    d.columns = [c.split(":")[0] if not c.startswith(":") else c for c in d.columns]
    return d


tis = N("Tissue"); grp = N("ExerciseGroup"); gene = N("Gene"); pw = N("Pathway"); ph = N("Phenotype")
pw = pw[~pw.database.isin(["GOBP", "GOCC", "GOMF"])].reset_index(drop=True)
FEAT = ["TranscriptFeature", "ProteinFeature", "PTMSite", "AffinityProtein", "ChromatinRegion",
        "MethylationRegion", "Metabolite"]
feats = pd.concat([N(l).assign(ftype=i) for i, l in enumerate(FEAT)], ignore_index=True)

ti = {k: i for i, k in enumerate(tis["id"])}
gi = {k: i for i, k in enumerate(grp["id"])}
gni = {k: i for i, k in enumerate(gene["id"])}
pwi = {k: i for i, k in enumerate(pw["id"])}
fi = {k: i for i, k in enumerate(feats["id"])}

# feature -> genes (best weight first)
m2g = E("MAPS_TO_GENE").sort_values("weight", ascending=False)
m2g = m2g[m2g[":START_ID"].isin(fi) & m2g[":END_ID"].isin(gni)]
fgenes = m2g.groupby(":START_ID")[":END_ID"].apply(lambda s: [gni[x] for x in s.unique()[:3]]).to_dict()
fg = [fgenes.get(f, []) for f in feats["id"]]

sym = gene.set_index("id").symbol.fillna("").to_dict()


def fname(r):
    sid = str(r.source_id)
    if r.ftype == 6:
        return r.refmet_name if isinstance(r.refmet_name, str) and r.refmet_name else sid
    g = fgenes.get(r.id)
    s = sym.get(gene["id"][g[0]], "") if g else ""
    if r.ftype == 2:  # PTM: gene + site
        site = sid.rsplit("_", 1)[-1]
        return f"{s or sid.rsplit('_', 1)[0]} {site}"
    if r.ftype in (4, 5):
        return f"{s} · {sid}" if s else sid
    return s or sid


feats["name"] = [fname(r) for r in feats.itertuples(index=False)]

TP = ["1w", "2w", "4w", "8w", "during_20_min", "during_40_min", "post_10_min", "post_15_30_45_min",
      "post_3.5_4_hr", "post_24_hr"]
tpi = {k: i for i, k in enumerate(TP)}
SEX = {"female": 0, "male": 1, "both": 2}

reg = pd.concat([E("UPREGULATED_IN"), E("DOWNREGULATED_IN")])
reg = reg[reg[":END_ID"].isin(ti)]           # organ edges only (they carry the group as a property)
reg = reg[reg[":START_ID"].isin(fi)]
R = {"f": reg[":START_ID"].map(fi).tolist(), "t": reg[":END_ID"].map(ti).tolist(),
     "g": reg["group"].map(gi).tolist(), "s": reg["sex"].map(SEX).tolist(),
     "tp": reg["timepoint"].map(tpi).tolist(),
     "l": np.where(reg[":TYPE"] == "UPREGULATED_IN", 1, -1).tolist(),          # direction only
     "q": (reg["neg_log10_adj_p"].clip(upper=300) * 100).round().astype(int).tolist(),   # -log10(adj p) x 100
     "z": (reg["z"] * 10).round().astype(int).tolist(),
     "w": (reg["weight"] * 100).round().astype(int).tolist(),
     "a": reg["assay"].astype(str).tolist()}
assays = sorted(set(R["a"])); ai = {a: i for i, a in enumerate(assays)}
R["a"] = [ai[a] for a in R["a"]]

mem = E("IN_PATHWAY")
gm = mem[mem[":START_ID"].isin(gni) & mem[":END_ID"].isin(pwi)]
members = [[] for _ in range(len(pw))]
for p, gsrc in zip(gm[":END_ID"].map(pwi), gm[":START_ID"].map(gni)):
    members[p].append(int(gsrc))
sm = mem[mem[":START_ID"].isin(fi) & mem[":END_ID"].isin(pwi)]
site_members = [[] for _ in range(len(pw))]
for p, f in zip(sm[":END_ID"].map(pwi), sm[":START_ID"].map(fi)):
    site_members[p].append(int(f))

en = pd.concat([E("ENRICHED_UP_IN").assign(d=1), E("ENRICHED_DOWN_IN").assign(d=-1)])
en = en[en[":END_ID"].isin(ti) & en[":START_ID"].isin(pwi)]
en["tpk"] = en.timepoint.astype(str)
PE = {"p": en[":START_ID"].map(pwi).tolist(), "t": en[":END_ID"].map(ti).tolist(),
      "g": en["group"].map(gi).tolist(), "s": en["sex"].map(SEX).tolist(), "tp": en["tpk"].tolist(),
      "d": en["d"].tolist(), "q": (-np.log10(en["adj_p_value"].clip(1e-300)) * 10).round().astype(int).tolist(),
      "a": en["assay"].astype(str).tolist()}

co = E("CO_REGULATED_WITH")
co = co[co[":START_ID"].isin(fi) & co[":END_ID"].isin(fi)]
CO = {"a": co[":START_ID"].map(fi).tolist(), "b": co[":END_ID"].map(fi).tolist(),
      "r": (co["pearson_r"] * 100).round().astype(int).tolist(), "t": (co["species"] + ":" + co["tissue"]).map(ti).tolist()}

orth = E("ORTHOLOG_OF")
orth = orth[orth[":START_ID"].isin(gni) & orth[":END_ID"].isin(gni)]
OS = E("ORTHOLOGOUS_SITE")
OS = OS[OS[":START_ID"].isin(fi) & OS[":END_ID"].isin(fi)]

def pairs(t, a_map, b_map):
    d = E(t)
    d = d[d[":START_ID"].isin(a_map) & d[":END_ID"].isin(b_map)]
    return [d[":START_ID"].map(a_map).tolist(), d[":END_ID"].map(b_map).tolist()]


pp = E("INTERACTS_WITH"); pp = pp[pp[":START_ID"].isin(gni) & pp[":END_ID"].isin(gni)]
PPI = [pp[":START_ID"].map(gni).tolist(), pp[":END_ID"].map(gni).tolist(),
       pp["source"].map({"Hetionet": 0, "STRING": 1, "Hetionet;STRING": 2}).tolist(),
       pp["string_score"].fillna(0).astype(int).tolist()]
LR = pairs("LIGAND_OF", gni, gni)
KIN = pairs("PHOSPHORYLATES", gni, fi)

phe = pd.concat([E("INCREASED_IN"), E("DECREASED_IN")])
phi = {k: i for i, k in enumerate(ph["id"])}
PH = {"p": phe[":START_ID"].map(phi).tolist(), "g": phe[":END_ID"].map(gi).tolist(),
      "s": phe["sex"].map(SEX).tolist(), "tp": phe["timepoint"].map(tpi).tolist(),
      "d": np.where(phe[":TYPE"] == "INCREASED_IN", 1, -1).tolist(),
      "v": phe["difference"].round(3).tolist() if "difference" in phe else [None]*len(phe),
      "ap": phe["adj_p_value"].map(lambda x: float(f"{x:.2g}")).tolist(),
      "p_value": phe["p_value"].map(lambda x: float(f"{x:.2g}")).tolist()}
phm = E("MEASURED_IN")

data = {
    "tissues": tis[["id", "name", "species"]].values.tolist(),
    "groups": grp[["id", "species", "modality", "regimen", "description"]].values.tolist(),
    "timepoints": TP, "assays": assays,
    "genes": {"id": gene["id"].tolist(), "sym": gene["symbol"].fillna("").tolist(),
              "sp": gene["species"].map({"rat": 0, "human": 1}).tolist()},
    "pathways": {"id": pw["id"].tolist(), "name": pw["name"].tolist(), "db": pw["database"].tolist(),
                 "members": members, "sites": site_members},
    "ftypes": FEAT,
    "features": {"id": feats["id"].tolist(), "name": feats["name"].tolist(), "type": feats["ftype"].tolist(),
                 "sp": feats["species"].map({"rat": 0, "human": 1, "rat;human": 2}).tolist(),
                 "genes": fg},
    "reg": R, "penr": PE, "coreg": CO,
    "ortho": [orth[":START_ID"].map(gni).tolist(), orth[":END_ID"].map(gni).tolist()],
    "ppi": PPI, "lr": LR, "kin": KIN,
    "orthosite": [OS[":START_ID"].map(fi).tolist(), OS[":END_ID"].map(fi).tolist()],
    "pheno": {"id": ph["id"].tolist(), "name": ph["name"].tolist(), "species": ph["species"].tolist(),
              "unit": ph["unit"].tolist(), "edges": PH,
              "measured": [phm[":START_ID"].map(phi).tolist(), phm[":END_ID"].map(ti).tolist()]},
}
def clean(x):
    if isinstance(x, dict):
        return {k: clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, float) and (x != x):
        return None
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return None if np.isnan(x) else float(x)
    return x


data["refs"] = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "refs", "references.json")))
data = clean(data)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
json.dump(data, open(OUT, "w"), separators=(",", ":"), allow_nan=False)
print(OUT, round(os.path.getsize(OUT) / 1e6, 2), "MB", "features", len(feats), "reg", len(R["f"]),
      "pathways", len(pw), "memberships", sum(map(len, members)))
