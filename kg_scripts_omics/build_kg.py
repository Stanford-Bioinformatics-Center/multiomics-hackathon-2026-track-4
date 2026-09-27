#!/usr/bin/env python3
"""
build_kg.py — Exercise-as-medicine multi-omics knowledge graph (schema v0.1)

Reads the tables converted from the MoTrPAC R packages
  rat   : MotrpacRatTraining6moData      -> data/processed/rat_pkg/*.tsv.gz
  human : MotrpacHumanPreSuspensionAnalysis (c2.0) -> data/processed/human_pkg/*.tsv.gz
and writes one CSV per node label and per edge type in Neo4j bulk-import
format (:ID / :LABEL / :START_ID / :END_ID / :TYPE), plus a manifest.

Usage:
  python build_kg.py --rat DIR --human DIR --out DIR [--human-padj 0.05] [--human-contrasts controls|all]

Edge weights (all in [0, 1]; raw statistics are kept alongside):
  RESPONDS_IN            weight = min(|z|, Z_CAP) / Z_CAP        (sign kept in `direction` and signed `z`)
  MAPS_TO_GENE           1.0 for gene products; genomic regions weighted by annotation (ANNOT_WEIGHT)
  ORTHOLOG_OF            1 / (number of human orthologs of the rat gene)
  ORTHOLOGOUS_SITE       1.0
  EQUIVALENT_TISSUE      1.0 exact anatomical match, 0.5 proxy
  MEMBER_OF_CLUSTER      rat 1.0 (hard state); human = fuzzy c-means membership
  *_ENRICHED_FOR         weight = min(-log10(adj_p), P_CAP) / P_CAP
"""
import argparse, json, math, os, re
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.stats import norm

Z_CAP = 10.0          # |z| at which a response edge saturates to weight 1
P_CAP = 10.0          # -log10(adj p) at which an enrichment edge saturates
ANNOT_WEIGHT = {      # ChIPseeker-style annotation of ATAC / methylation regions
    "Promoter (<=1kb)": 1.0, "Promoter (1-2kb)": 0.8, "5' UTR": 0.8, "Exon": 0.6,
    "Intron": 0.5, "3' UTR": 0.5, "Upstream (<5kb)": 0.5, "Downstream (<5kb)": 0.4,
    "Distal Intergenic": 0.25, "Overlaps Gene": 0.6,
}

RAT_WEEKS = ["1w", "2w", "4w", "8w"]
HUMAN_TP_ORDER = ["pre_exercise", "during_20_min", "during_40_min", "post_10_min",
                  "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr"]
HUMAN_GROUP = {"ADUEndur": "endurance", "ADUResist": "resistance", "ADUControl": "control"}
TISSUE_EQUIV = [  # rat, human, match
    ("SKM-VL", "muscle", "exact"), ("SKM-GN", "muscle", "proxy"),
    ("WAT-SC", "adipose", "exact"), ("BLOOD", "blood", "proxy"), ("PLASMA", "blood", "proxy"),
]
UBERON = {"SKM-GN": "UBERON:0001388", "SKM-VL": "UBERON:0001379", "muscle": "UBERON:0001379",
          "WAT-SC": "UBERON:0002190", "adipose": "UBERON:0002190", "HEART": "UBERON:0000948",
          "LIVER": "UBERON:0002107", "LUNG": "UBERON:0002048", "KIDNEY": "UBERON:0002113",
          "BAT": "UBERON:0001348", "BLOOD": "UBERON:0000178", "blood": "UBERON:0000178",
          "PLASMA": "UBERON:0001969", "CORTEX": "UBERON:0001851", "HIPPOC": "UBERON:0002421",
          "HYPOTH": "UBERON:0001898", "ADRNL": "UBERON:0002369", "COLON": "UBERON:0001155",
          "SMLINT": "UBERON:0002108", "SPLEEN": "UBERON:0002106", "OVARY": "UBERON:0000992",
          "TESTES": "UBERON:0000473", "VENACV": "UBERON:0001585"}

# rat ome -> (node label, human-comparable ome name)
RAT_OME_LABEL = {"TRNSCRPT": "TranscriptFeature", "PROT": "ProteinFeature", "PHOSPHO": "PTMSite",
                 "ACETYL": "PTMSite", "UBIQ": "PTMSite", "IMMUNO": "AffinityProtein",
                 "METAB": "Metabolite", "ATAC": "ChromatinRegion", "METHYL": "MethylationRegion"}
HUMAN_ASSAY_LABEL = {"transcript-rna-seq": "TranscriptFeature", "prot-pr": "ProteinFeature",
                     "prot-ph": "PTMSite", "prot-ol": "AffinityProtein", "metab": "Metabolite",
                     "epigen-atac-seq": "ChromatinRegion", "epigen-methylcap-seq": "MethylationRegion",
                     "prot-clinical": "ClinicalAnalyte", "metab-t-clinical": "ClinicalAnalyte"}


def rd(path, **kw):
    return pd.read_csv(path, sep="\t", low_memory=False, **kw)


def strip_ver(x):
    return x.split(".")[0] if isinstance(x, str) and x.startswith("ENS") else x


def w_from_z(z):
    return np.minimum(np.abs(z), Z_CAP) / Z_CAP


def w_from_p(p):
    return np.minimum(-np.log10(np.clip(p, 1e-300, 1)), P_CAP) / P_CAP


class Graph:
    """Accumulates node/edge rows per label/type and de-duplicates on write."""
    def __init__(self):
        self.nodes = defaultdict(list)
        self.edges = defaultdict(list)

    def node(self, label, df):
        self.nodes[label].append(df)

    def edge(self, etype, df):
        self.edges[etype].append(df)

    def node_ids(self, label=None):
        labels = [label] if label else self.nodes
        out = set()
        for l in labels:
            for d in self.nodes.get(l, []):
                out.update(d["id"])
        return out

    def write(self, out):
        os.makedirs(out, exist_ok=True)
        manifest = {"nodes": {}, "edges": {}}
        for label, parts in self.nodes.items():
            d = pd.concat(parts, ignore_index=True).drop_duplicates("id")
            d = d.rename(columns={"id": "id:ID"})
            d.insert(1, ":LABEL", "MolecularFeature;" + label if label in FEATURE_LABELS else label)
            d.to_csv(f"{out}/nodes_{label}.csv.gz", index=False)
            manifest["nodes"][label] = len(d)
        all_ids = self.node_ids()
        for etype, parts in self.edges.items():
            d = pd.concat(parts, ignore_index=True)
            before = len(d)
            d = d[d["src"].isin(all_ids) & d["dst"].isin(all_ids)].drop_duplicates()
            dropped = before - len(d)
            d = d.rename(columns={"src": ":START_ID", "dst": ":END_ID"})
            d.insert(2, ":TYPE", etype)
            d.to_csv(f"{out}/edges_{etype}.csv.gz", index=False)
            manifest["edges"][etype] = {"n": len(d), "dropped_dangling_or_duplicate": dropped}
        json.dump(manifest, open(f"{out}/manifest.json", "w"), indent=2)
        return manifest


FEATURE_LABELS = set(RAT_OME_LABEL.values()) | set(HUMAN_ASSAY_LABEL.values())


# ----------------------------------------------------------------------------- static nodes
def add_static(g):
    g.node("Species", pd.DataFrame({"id": ["NCBITaxon:10116", "NCBITaxon:9606"],
                                    "name": ["Rattus norvegicus", "Homo sapiens"]}))
    g.node("Study", pd.DataFrame({
        "id": ["motrpac:rat-training-06", "motrpac:human-precovid-sed-adu"],
        "design": ["8-week progressive treadmill endurance training, 6-month-old F344 rats",
                   "single acute endurance or resistance bout vs non-exercise control, sedentary adults (pre-COVID)"],
        "source": ["MotrpacRatTraining6moData", "MotrpacHumanPreSuspensionAnalysis c2.0"],
        "species": ["rat", "human"]}))
    g.edge("IN_SPECIES", pd.DataFrame({"src": ["motrpac:rat-training-06", "motrpac:human-precovid-sed-adu"],
                                       "dst": ["NCBITaxon:10116", "NCBITaxon:9606"], "weight": 1.0}))
    g.node("ExerciseModality", pd.DataFrame({"id": ["modality:endurance", "modality:resistance", "modality:control"],
                                             "name": ["endurance", "resistance", "control (no exercise)"]}))
    tissues = pd.DataFrame(
        [(f"rat:{t}", t, "rat", UBERON.get(t, "")) for t in UBERON if t.isupper() or "-" in t] +
        [(f"human:{t}", t, "human", UBERON[t]) for t in ["muscle", "adipose", "blood"]],
        columns=["id", "name", "species", "uberon_id"])
    g.node("Tissue", tissues)
    g.edge("EQUIVALENT_TISSUE", pd.DataFrame(
        [(f"rat:{r}", f"human:{h}", m, 1.0 if m == "exact" else 0.5) for r, h, m in TISSUE_EQUIV],
        columns=["src", "dst", "match", "weight"]))


# ----------------------------------------------------------------------------- rat
def add_rat(g, rat):
    trf = rd(f"{rat}/TRAINING_REGULATED_FEATURES.tsv.gz")
    f2g = rd(f"{rat}/FEATURE_TO_GENE_FILT.tsv.gz")
    mmap = rd(f"{rat}/METAB_FEATURE_ID_MAP.tsv.gz")

    # metabolites -> species-agnostic RefMet node where possible
    refmet = (mmap.dropna(subset=["metabolite_refmet"])
                  .drop_duplicates(["tissue", "feature_ID_metareg"])
                  .set_index(["tissue", "feature_ID_metareg"])["metabolite_refmet"].to_dict())

    def fid(row_assay, row_tissue, row_id):
        if row_assay == "METAB":
            rm = refmet.get((row_tissue, row_id))
            if isinstance(rm, str):
                return f"REFMET:{rm}"
        return f"rat:{row_assay}:{row_id}"

    trf["fid"] = [fid(a, t, i) for a, t, i in zip(trf.assay, trf.tissue, trf.feature_ID)]

    # feature nodes
    feats = trf.drop_duplicates("fid")[["fid", "assay", "assay_code", "feature_ID", "platform"]]
    for ome, sub in feats.groupby("assay"):
        label = RAT_OME_LABEL[ome]
        d = pd.DataFrame({"id": sub.fid, "species": np.where(sub.fid.str.startswith("REFMET:"), "shared", "rat"),
                          "ome": ome, "assay_code": sub.assay_code, "source_id": sub.feature_ID,
                          "platform": sub.platform.fillna("")})
        if label == "PTMSite":
            d["ptm_type"] = ome
            d["protein_id"] = sub.feature_ID.str.rsplit("_", n=1).str[0].values
        if label in ("ChromatinRegion",):
            m = sub.feature_ID.str.extract(r"^(chr[^:]+):(\d+)-(\d+)$")
            d["chrom"], d["start"], d["end"] = m[0].values, m[1].values, m[2].values
        if label == "Metabolite":
            d["refmet_name"] = np.where(sub.fid.str.startswith("REFMET:"), sub.fid.str[7:], "")
        g.node(label, d)

    # contrast nodes: tissue x sex x week (trained vs sedentary control)
    con = trf[["tissue", "sex", "training_group"]].drop_duplicates()
    con["id"] = "rat:" + con.tissue + ":" + con.sex + ":" + con.training_group
    g.node("Contrast", pd.DataFrame({
        "id": con.id, "species": "rat", "tissue": con.tissue, "sex": con.sex,
        "group": "endurance", "timepoint": con.training_group,
        "timepoint_order": con.training_group.map({w: i for i, w in enumerate(RAT_WEEKS)}),
        "regimen": "chronic_training", "comparison": "trained vs sedentary control"}))
    g.edge("IN_TISSUE", pd.DataFrame({"src": con.id, "dst": "rat:" + con.tissue, "weight": 1.0}))
    g.edge("PART_OF_STUDY", pd.DataFrame({"src": con.id, "dst": "motrpac:rat-training-06", "weight": 1.0}))
    g.edge("OF_MODALITY", pd.DataFrame({"src": con.id, "dst": "modality:endurance", "weight": 1.0}))

    # RESPONDS_IN (skip the handful of rows with no estimate for that sex/week)
    trf = trf.dropna(subset=["timewise_logFC"])
    z = trf.timewise_zscore.astype(float)
    z = z.fillna(trf.timewise_logFC / trf.timewise_logFC_se)          # a few rows lack a z-score
    z = z.fillna(pd.Series(np.sign(trf.timewise_logFC) * norm.isf(trf.timewise_p_value / 2), index=trf.index))
    g.edge("RESPONDS_IN", pd.DataFrame({
        "src": trf.fid, "dst": "rat:" + trf.tissue + ":" + trf.sex + ":" + trf.training_group,
        "assay": trf.assay, "logFC:float": trf.timewise_logFC, "logFC_se:float": trf.timewise_logFC_se,
        "p_value:float": trf.timewise_p_value, "z:float": z, "training_q:float": trf.training_q,
        "direction": np.where(trf.timewise_logFC > 0, "up", "down"), "weight:float": w_from_z(z)}))

    # MAPS_TO_GENE
    keep = f2g[f2g.feature_ID.isin(trf.feature_ID)].dropna(subset=["ensembl_gene"]).copy()
    id2fids = trf.groupby("feature_ID").fid.unique().to_dict()
    rows = []
    for r in keep.itertuples(index=False):
        w = ANNOT_WEIGHT.get(r.custom_annotation, 1.0) if isinstance(r.custom_annotation, str) else 1.0
        for f in id2fids.get(r.feature_ID, []):
            rows.append((f, f"ENSEMBL:{r.ensembl_gene}", r.custom_annotation if isinstance(r.custom_annotation, str) else "",
                         r.relationship_to_gene, w))
    m2g = pd.DataFrame(rows, columns=["src", "dst", "annotation", "distance_to_gene", "weight:float"])
    g.edge("MAPS_TO_GENE", m2g)
    genes = keep.drop_duplicates("ensembl_gene")
    g.node("Gene", pd.DataFrame({"id": "ENSEMBL:" + genes.ensembl_gene, "symbol": genes.gene_symbol,
                                 "species": "rat", "entrez_id": genes.entrez_gene.astype("Int64").astype(str),
                                 "rgd_id": genes.rgd_gene.astype("Int64").astype(str)}))

    # SITE_ON (rat PTM site -> rat protein, when the protein is itself a node)
    ptm = feats[feats.assay.isin(["PHOSPHO", "ACETYL", "UBIQ"])]
    g.edge("SITE_ON", pd.DataFrame({"src": ptm.fid,
                                    "dst": "rat:PROT:" + ptm.feature_ID.str.rsplit("_", n=1).str[0],
                                    "weight": 1.0}))

    # response clusters: rat graphical states (tissue:week_state)
    gs = rd(f"{rat}/GRAPH_STATES.tsv.gz")
    gs["fid"] = [fid(o, t, i) for o, t, i in zip(gs.ome, gs.tissue, gs.feature_ID)]
    long = gs.melt(id_vars=["fid", "tissue", "ome"], value_vars=[f"state_{w}" for w in RAT_WEEKS],
                   var_name="week", value_name="state")
    long["week"] = long.week.str.replace("state_", "")
    long = long.dropna(subset=["state"])   # F0_M0 kept: MoTrPAC enriches those clusters too
    long["cid"] = "rat:graph:" + long.tissue + ":" + long.week + "_" + long.state
    cl = long.groupby(["cid", "tissue", "week", "state"]).size().reset_index(name="size")
    g.node("ResponseCluster", pd.DataFrame({"id": cl.cid, "species": "rat", "method": "graphical_state",
                                            "tissue": cl.tissue, "timepoint": cl.week, "state": cl.state,
                                            "size:int": cl["size"]}))
    g.edge("MEMBER_OF_CLUSTER", pd.DataFrame({"src": long.fid, "dst": long.cid, "ome": long.ome, "weight:float": 1.0}))
    # full 1w->8w trajectories (MoTrPAC's main "graphical clusters")
    gs["pid"] = "rat:graph:" + gs.tissue_path
    pz = gs.groupby(["pid", "tissue", "path"]).size().reset_index(name="size")
    g.node("ResponseCluster", pd.DataFrame({"id": pz.pid, "species": "rat", "method": "graphical_path",
                                            "tissue": pz.tissue, "timepoint": "1w-8w", "state": pz.path,
                                            "size:int": pz["size"]}))
    g.edge("MEMBER_OF_CLUSTER", pd.DataFrame({"src": gs.fid, "dst": gs.pid, "ome": gs.ome, "weight:float": 1.0}))

    pw = rd(f"{rat}/GRAPH_PW_ENRICH.tsv.gz")
    pw = pw[(pw.adj_p_value < 0.05) & ~pw.cluster.str.contains("---", na=False)]
    g.node("GeneSet", pd.DataFrame({"id": pw.term_id, "name": pw.term_name, "database": pw.source,
                                    "size:int": pw.term_size}))
    g.edge("CLUSTER_ENRICHED_FOR", pd.DataFrame({
        "src": "rat:graph:" + pw.cluster, "dst": pw.term_id, "ome": pw.ome, "method": "gprofiler",
        "intersection_size:int": pw.intersection_size, "adj_p_value:float": pw.adj_p_value,
        "weight:float": w_from_p(pw.adj_p_value)}))
    return trf


def add_orthology(g, rat):
    o = rd(f"{rat}/RAT_TO_HUMAN_GENE.tsv.gz").dropna(subset=["RAT_ENSEMBL_ID", "HUMAN_ORTHOLOG_ENSEMBL_ID"])
    o = o.assign(RAT_ENSEMBL_ID=o.RAT_ENSEMBL_ID.str.split(";")).explode("RAT_ENSEMBL_ID")
    o["h"] = o.HUMAN_ORTHOLOG_ENSEMBL_ID.map(strip_ver)
    o["src"] = "ENSEMBL:" + o.RAT_ENSEMBL_ID
    o["dst"] = "ENSEMBL:" + o.h
    rat_genes = g.node_ids("Gene")
    o = o[o.src.isin(rat_genes)].drop_duplicates(["src", "dst"])
    n = o.groupby("src").dst.transform("count")
    g.edge("ORTHOLOG_OF", pd.DataFrame({"src": o.src, "dst": o.dst, "source": o.HUMAN_ORTHOLOG_SOURCE,
                                        "weight:float": 1.0 / n}))
    # make sure the human ortholog exists as a Gene node even if it has no human response edge
    g.node("Gene", pd.DataFrame({"id": o.dst, "symbol": o.HUMAN_ORTHOLOG_SYMBOL, "species": "human",
                                 "entrez_id": o.HUMAN_ORTHOLOG_NCBI_GENE_ID.astype(str), "rgd_id": ""}))

    ps = rd(f"{rat}/RAT_TO_HUMAN_PHOSPHO.tsv.gz").dropna()
    g.edge("ORTHOLOGOUS_SITE", pd.DataFrame({"src": "rat:PHOSPHO:" + ps.ptm_id_rat_refseq,
                                             "dst": "human:prot-ph:" + ps.ptm_id_human_uniprot, "weight": 1.0}))


# ----------------------------------------------------------------------------- human
def add_human(g, hum, padj, contrasts):
    f2g = rd(f"{hum}/HUMAN_FEATURE_TO_GENE.tsv.gz")
    metab_refmet = (f2g[f2g.assay == "metab"].dropna(subset=["refmet_name"])
                    .drop_duplicates("feature_id").set_index("feature_id").refmet_name.to_dict())

    def fid(assay, feature_id):
        if assay == "metab":
            return f"REFMET:{metab_refmet.get(feature_id, feature_id)}"
        return f"human:{assay}:{feature_id}"

    keep_types = {"controls": ["exercise_with_controls"],
                  "all": ["exercise_with_controls", "exercise_no_controls"]}[contrasts]
    parts = []
    for f in sorted(os.listdir(hum)):
        if not f.endswith("_DA.tsv.gz"):
            continue
        d = rd(f"{hum}/{f}", usecols=lambda c: c in {
            "tissue", "assay", "contrast_type", "contrast_short", "randomGroupCode", "Timepoint",
            "feature_id", "logFC", "CI.L_calculated", "CI.R_calculated", "z.std", "p_value", "adj_p_value"})
        d = d[d.contrast_type.isin(keep_types) & (d.adj_p_value < padj)]
        parts.append(d)
    da = pd.concat(parts, ignore_index=True)
    # clinical assays share the generic 'metab'/'prot' names in some tables -> keep the file's assay
    da["fid"] = [fid(a, i) for a, i in zip(da.assay, da.feature_id)]
    da["cid"] = ("human:" + da.tissue + ":" + da.randomGroupCode.map(HUMAN_GROUP) + ":" + da.Timepoint +
                 np.where(da.contrast_type == "exercise_with_controls", ":vs_control", ":vs_pre"))

    # feature nodes
    feats = da.drop_duplicates("fid")[["fid", "assay", "feature_id"]]
    for assay, sub in feats.groupby("assay"):
        label = HUMAN_ASSAY_LABEL.get(assay, "ClinicalAnalyte")
        d = pd.DataFrame({"id": sub.fid, "species": np.where(sub.fid.str.startswith("REFMET:"), "shared", "human"),
                          "ome": assay, "assay_code": assay, "source_id": sub.feature_id})
        if label == "PTMSite":
            d["ptm_type"] = "PHOSPHO"
            d["protein_id"] = sub.feature_id.str.rsplit("_", n=1).str[0].values
        if label == "Metabolite":
            d["refmet_name"] = sub.fid.str[7:].values
        g.node(label, d)

    # contrast nodes
    con = da.drop_duplicates("cid")
    g.node("Contrast", pd.DataFrame({
        "id": con.cid, "species": "human", "tissue": con.tissue, "sex": "both",
        "group": con.randomGroupCode.map(HUMAN_GROUP), "timepoint": con.Timepoint,
        "timepoint_order": con.Timepoint.map({t: i for i, t in enumerate(HUMAN_TP_ORDER)}),
        "regimen": "acute", "comparison": con.contrast_short}))
    g.edge("IN_TISSUE", pd.DataFrame({"src": con.cid, "dst": "human:" + con.tissue, "weight": 1.0}))
    g.edge("PART_OF_STUDY", pd.DataFrame({"src": con.cid, "dst": "motrpac:human-precovid-sed-adu", "weight": 1.0}))
    g.edge("OF_MODALITY", pd.DataFrame({"src": con.cid, "dst": "modality:" + con.randomGroupCode.map(HUMAN_GROUP),
                                        "weight": 1.0}))

    # RESPONDS_IN
    z = da["z.std"].astype(float)
    g.edge("RESPONDS_IN", pd.DataFrame({
        "src": da.fid, "dst": da.cid, "assay": da.assay, "logFC:float": da.logFC,
        "ci_low:float": da["CI.L_calculated"], "ci_high:float": da["CI.R_calculated"],
        "p_value:float": da.p_value, "adj_p_value:float": da.adj_p_value, "z:float": z,
        "direction": np.where(da.logFC > 0, "up", "down"), "weight:float": w_from_z(z)}))

    # MAPS_TO_GENE
    m = f2g.merge(feats.rename(columns={"feature_id": "fid_src"}), left_on=["assay", "feature_id"],
                  right_on=["assay", "fid_src"]).dropna(subset=["ensembl_gene"])
    m["gene"] = "ENSEMBL:" + m.ensembl_gene.map(strip_ver)
    m["w"] = [ANNOT_WEIGHT.get(a, 1.0) if isinstance(a, str) else 1.0 for a in m.custom_annotation]
    g.edge("MAPS_TO_GENE", pd.DataFrame({"src": m.fid, "dst": m.gene, "annotation": m.custom_annotation.fillna(""),
                                         "distance_to_gene": m.relationship_to_gene, "weight:float": m.w}))
    tfs = set(rd(f"{hum}/UTORONTO_TFs.tsv.gz").gene_symbol)
    gn = m.drop_duplicates("gene")
    g.node("Gene", pd.DataFrame({"id": gn.gene, "symbol": gn.gene_symbol, "species": "human",
                                 "entrez_id": gn.entrez_gene.astype("Int64").astype(str), "rgd_id": "",
                                 "is_TF:boolean": gn.gene_symbol.isin(tfs)}))

    # SITE_ON (human phosphosite -> protein, if the protein is a node)
    ph = feats[feats.assay == "prot-ph"]
    g.edge("SITE_ON", pd.DataFrame({"src": ph.fid, "dst": "human:prot-pr:" + ph.feature_id.str.rsplit("_", n=1).str[0],
                                    "weight": 1.0}))

    # pathway / kinase enrichment per contrast
    cam = rd(f"{hum}/CAMERA_RESULTS.tsv.gz")
    ptm = rd(f"{hum}/PTMSEA_RESULTS.tsv.gz")
    enr = []
    for d, meth, zc in [(cam, "CAMERA-PR", "z.std"), (ptm, "PTM-SEA", "NES")]:
        d = d[d.contrast_type.isin(keep_types) & (d.adj_p_value < 0.05)].merge(
            rd(f"{hum}/CONTRAST_CONVERTER.tsv.gz")[["contrast", "randomGroupCode", "Timepoint"]], on="contrast")
        d["cid"] = ("human:" + d.tissue + ":" + d.randomGroupCode.map(HUMAN_GROUP) + ":" + d.Timepoint +
                    np.where(d.contrast_type == "exercise_with_controls", ":vs_control", ":vs_pre"))
        d["gsid"] = d.database + ":" + d.set_id.astype(str)
        g.node("GeneSet", pd.DataFrame({"id": d.gsid, "name": d.set_short, "database": d.database,
                                        "collection": d.collection, "size:int": d.set_size_DB}))
        enr.append(pd.DataFrame({"src": d.cid, "dst": d.gsid, "assay": d.assay, "method": meth,
                                 "direction": d.direction, "score:float": d[zc], "adj_p_value:float": d.adj_p_value,
                                 "weight:float": w_from_p(d.adj_p_value)}))
    g.edge("CONTRAST_ENRICHED_FOR", pd.concat(enr, ignore_index=True))

    # fuzzy c-means clusters (per tissue, across omes)
    fm = rd(f"{hum}/FCM_MEMBERSHIP.tsv.gz")
    fm["fid"] = [fid(a, i) for a, i in zip(fm.assay, fm.feature_id)]
    fm["cid"] = "human:fcm:" + fm.tissue + ":" + fm.cluster.astype(str)
    sz = fm.groupby(["cid", "tissue", "cluster"]).size().reset_index(name="size")
    g.node("ResponseCluster", pd.DataFrame({"id": sz.cid, "species": "human", "method": "fuzzy_c_means",
                                            "tissue": sz.tissue, "timepoint": "", "state": "cluster " + sz.cluster.astype(str),
                                            "size:int": sz["size"]}))
    g.edge("MEMBER_OF_CLUSTER", pd.DataFrame({"src": fm.fid, "dst": fm.cid, "ome": fm.assay,
                                              "weight:float": fm.membership}))
    ora = rd(f"{hum}/FCM_ORA.tsv.gz")
    ora = ora[ora.adj_p_value < 0.05]
    ora["gsid"] = ora.database + ":" + ora.set_id.astype(str)
    g.node("GeneSet", pd.DataFrame({"id": ora.gsid, "name": ora.set_short, "database": ora.database,
                                    "collection": ora.collection, "size:int": ora.set_size_DB}))
    g.edge("CLUSTER_ENRICHED_FOR", pd.DataFrame({
        "src": "human:fcm:" + ora.tissue + ":" + ora.cluster.astype(str), "dst": ora.gsid, "ome": ora.assay,
        "method": "ORA", "intersection_size:int": ora.set_size_in_cluster,
        "adj_p_value:float": ora.adj_p_value, "weight:float": w_from_p(ora.adj_p_value)}))
    return da


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rat", required=True)
    ap.add_argument("--human", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--human-padj", type=float, default=0.05)
    ap.add_argument("--human-contrasts", choices=["controls", "all"], default="controls")
    a = ap.parse_args()

    g = Graph()
    add_static(g)
    add_rat(g, a.rat)
    add_human(g, a.human, a.human_padj, a.human_contrasts)
    add_orthology(g, a.rat)          # after human so human Gene rows (with is_TF) win de-duplication
    man = g.write(a.out)
    print(json.dumps(man, indent=2))


if __name__ == "__main__":
    main()
