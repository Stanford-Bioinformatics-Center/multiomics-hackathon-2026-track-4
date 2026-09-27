#!/usr/bin/env python3
"""
build_sig_graph.py — significant-response, pathway-anchored exercise KG (schema v0.2)

Only molecules that are SIGNIFICANTLY up- or down-regulated are included:
  rat   : per sex x week call from MoTrPAC graphical states (GRAPH_STATES, F/M = -1/0/1);
          features without a state (3%, mostly ovary/testes) fall back to timewise p < 0.05
  human : control-adjusted (delta-delta) contrasts, BH adj p < 0.05

Every regulation edge carries: species, tissue, group, sex, timepoint, timepoint_order,
timescale, assay, ome, logFC, z, p, (adj_p), weight.

Metabolite rules (per project decision):
  * NO metabolite -> ExerciseGroup edges
  * NO metabolite -> Gene edges, NO metabolite -> Pathway edges
  * Metabolites connect only to: Tissue (organ; UP/DOWNREGULATED_IN, per timepoint),
    other Metabolites and Proteins (CO_REGULATED_WITH, response-profile correlation).
  * Human clinical-chemistry analytes (lactate, glucose, cortisol, NEFA, glycerol, ketones,
    insulin, glucagon, CK) are modelled as Phenotypes, not metabolites/proteins.

Usage:
  python build_sig_graph.py --rat DIR --human DIR --out DIR [--ppi-dir DIR] [--r-min 0.9] [--top-k 10] [--go enriched|all|none]

Interaction layers (with --ppi-dir): INTERACTS_WITH (human physical PPI, Hetionet), LIGAND_OF (LIANA consensus),
PHOSPHORYLATES (PhosphoSitePlus kinase -> significant phosphosite). All human evidence; rat reaches it via orthologs.
Requires: pandas, numpy, scipy; R (only to read MOLECULAR_SIGNATURES.rds once -> cached TSV)
"""
import argparse, json, os, subprocess
import numpy as np
import pandas as pd
from scipy import stats

from build_kg import (Graph, rd, strip_ver, w_from_z, w_from_p, ANNOT_WEIGHT, UBERON,
                      TISSUE_EQUIV, HUMAN_GROUP, HUMAN_TP_ORDER, RAT_WEEKS)

RAT_ASSAY = {"TRNSCRPT": ("TranscriptFeature", "transcriptome"), "PROT": ("ProteinFeature", "proteome"),
             "PHOSPHO": ("PTMSite", "phosphoproteome"), "ACETYL": ("PTMSite", "acetylome"),
             "UBIQ": ("PTMSite", "ubiquitylome"), "IMMUNO": ("AffinityProtein", "proteome-targeted"),
             "METAB": ("Metabolite", "metabolome"), "ATAC": ("ChromatinRegion", "epigenome-chromatin"),
             "METHYL": ("MethylationRegion", "epigenome-methylation")}
HUM_ASSAY = {"transcript-rna-seq": ("TranscriptFeature", "transcriptome"),
             "prot-pr": ("ProteinFeature", "proteome"), "prot-ph": ("PTMSite", "phosphoproteome"),
             "prot-ol": ("AffinityProtein", "proteome-targeted"), "metab": ("Metabolite", "metabolome")}
PROTEIN_LABELS = {"ProteinFeature", "AffinityProtein"}
CURATED = ["KEGG_MEDICUS", "REACTOME", "WP", "BIOCARTA", "PID", "MITOCARTA"]
GO = ["GOBP", "GOCC", "GOMF"]
PTM_SETS = ["PSP", "PTMSIGDB"]
GROUPS = pd.DataFrame([
    ("rat:endurance_training", "rat", "endurance", "chronic_training", "8-wk progressive treadmill training vs sedentary control"),
    ("human:endurance_acute", "human", "endurance", "acute", "single cycling bout vs non-exercise control (delta-delta)"),
    ("human:resistance_acute", "human", "resistance", "acute", "single resistance bout vs non-exercise control (delta-delta)")],
    columns=["id", "species", "modality", "regimen", "description"])
HUMAN_GROUP_NODE = {"ADUEndur": "human:endurance_acute", "ADUResist": "human:resistance_acute"}
HUMAN_CLINICAL = {"metab-t-clinical", "prot-clinical"}


def reg_type(logfc):
    return np.where(np.asarray(logfc) > 0, "UPREGULATED_IN", "DOWNREGULATED_IN")


def add_reg_edges(g, df, to_group=True):
    """df columns: fid, species, tissue, group_id, sex, timepoint, timepoint_order, timescale,
    assay, ome, logFC, z, p_value, adj_p_value, is_metab"""
    az = df.z.astype(float).abs()
    nl = -np.log10(df.adj_p_value.astype(float).clip(lower=1e-300))
    # tie-break equal adj p (e.g. feature-level q) by |z| so ranks stay informative
    pct = (nl + az * 1e-6).groupby([df.species, df.ome, df.tissue]).rank(pct=True, method="average")
    base = pd.DataFrame({
        "src": df.fid, "species": df.species, "tissue": df.tissue, "group": df.group_id, "sex": df.sex,
        "timepoint": df.timepoint, "timepoint_order:int": df.timepoint_order, "timescale": df.timescale,
        "assay": df.assay, "ome": df.ome, "adj_p_value:float": df.adj_p_value, "neg_log10_adj_p:float": nl.round(3),
        "adj_p_type": df.adj_p_type, "z:float": df.z, "p_value:float": df.p_value,
        "training_q:float": df.get("training_q", np.nan), "timewise_adj_p:float": df.get("timewise_adj_p", np.nan),
        "significance_basis": df.significance_basis,
        # weight = percentile of -log10(adj p) among significant changes of the same species x assay x organ
        "weight:float": pct.round(4), "_type": reg_type(df.logFC)})
    # feature -> organ (all omics, incl. metabolites)
    t = base.assign(dst=df.species + ":" + df.tissue)
    for et, sub in t.groupby("_type"):
        g.edge(et, sub.drop(columns="_type"))
    # feature -> exercise group (never for metabolites)
    if to_group:
        grp = base[~df.is_metab.values].assign(dst=df.group_id[~df.is_metab.values].values)
        for et, sub in grp.groupby("_type"):
            g.edge(et, sub.drop(columns="_type"))


# ------------------------------------------------------------------------------------ static
def add_static(g):
    g.node("ExerciseGroup", GROUPS.copy())
    tis = [(f"rat:{t}", t, "rat", UBERON.get(t, "")) for t in UBERON if t.isupper() or "-" in t]
    tis += [(f"human:{t}", t, "human", UBERON[t]) for t in ["muscle", "adipose", "blood"]]
    g.node("Tissue", pd.DataFrame(tis, columns=["id", "name", "species", "uberon_id"]))
    g.edge("EQUIVALENT_TISSUE", pd.DataFrame(
        [(f"rat:{r}", f"human:{h}", m, 1.0 if m == "exact" else 0.5) for r, h, m in TISSUE_EQUIV],
        columns=["src", "dst", "match", "weight:float"]))


# ------------------------------------------------------------------------------------ rat
def rat_significant(rat):
    trf = rd(f"{rat}/TRAINING_REGULATED_FEATURES.tsv.gz").dropna(subset=["timewise_logFC"])
    gs = rd(f"{rat}/GRAPH_STATES.tsv.gz")
    long = gs.melt(id_vars=["feature"], value_vars=[f"state_{w}" for w in RAT_WEEKS],
                   var_name="week", value_name="state").dropna()
    long["week"] = long.week.str[6:]
    s = long.state.str.extract(r"F(-?\d)_M(-?\d)").astype(int)
    calls = pd.concat([long.assign(sex="female", call=s[0].values), long.assign(sex="male", call=s[1].values)])
    calls = calls[calls.call != 0][["feature", "week", "sex", "call"]]
    m = trf.merge(calls, left_on=["feature", "training_group", "sex"], right_on=["feature", "week", "sex"], how="left")
    has_state = trf.feature.isin(gs.feature)
    m["has_state"] = m.feature.isin(gs.feature)
    sig = m[(m.call.notna()) | (~m.has_state & (m.timewise_p_value < 0.05))].copy()
    sig["significance_basis"] = np.where(sig.call.notna(), "graphical_state", "timewise_p<0.05")
    z = sig.timewise_zscore.astype(float).fillna(sig.timewise_logFC / sig.timewise_logFC_se)
    sig["z"] = z
    # timewise BH-adjusted p (per tissue x assay x sex x week) from MoTrPAC's per-tissue DA tables;
    # metabolites use the meta-regression tables (TRF metabolite IDs are meta-regression IDs).
    import glob
    parts = []
    for f in glob.glob(f"{rat}/*_DA.tsv.gz") + glob.glob(f"{rat}/*_DA_METAREG.tsv.gz"):
        name = os.path.basename(f)
        if name.startswith("METAB_") and not name.endswith("_METAREG.tsv.gz"):
            continue
        d = rd(f, usecols=["tissue", "feature_ID", "sex", "comparison_group", "adj_p_value"])
        parts.append(d)
    da = pd.concat(parts).groupby(["tissue", "feature_ID", "sex", "comparison_group"], as_index=False).adj_p_value.min()
    sig = sig.merge(da.rename(columns={"comparison_group": "training_group", "adj_p_value": "timewise_adj_p"}),
                    on=["tissue", "feature_ID", "sex", "training_group"], how="left")
    # MoTrPAC's significance for rat is the feature-level training q (one test across all sex x week groups);
    # per-week BH adj p is mostly ~1 at n~5/group and is kept only as an extra property.
    sig["adj_p_type"] = "training q (MoTrPAC feature-level FDR across sexes and weeks)"
    sig["adj_p_used"] = sig.training_q
    return trf, sig


def add_rat(g, rat):
    trf, sig = rat_significant(rat)
    f2g = rd(f"{rat}/FEATURE_TO_GENE_FILT.tsv.gz")
    mmap = rd(f"{rat}/METAB_FEATURE_ID_MAP.tsv.gz").dropna(subset=["metabolite_refmet"])
    refmet = mmap.drop_duplicates(["tissue", "feature_ID_metareg"]).set_index(
        ["tissue", "feature_ID_metareg"]).metabolite_refmet.to_dict()

    def fid(a, t, i):
        if a == "METAB" and isinstance(refmet.get((t, i)), str):
            return f"REFMET:{refmet[(t, i)]}"
        return f"rat:{a}:{i}"

    for d in (trf, sig):
        d["fid"] = [fid(a, t, i) for a, t, i in zip(d.assay, d.tissue, d.feature_ID)]

    # nodes
    fe = sig.drop_duplicates("fid")
    for ome, sub in fe.groupby("assay"):
        label, omics = RAT_ASSAY[ome]
        d = pd.DataFrame({"id": sub.fid, "species": np.where(sub.fid.str.startswith("REFMET:"), "rat;human", "rat"),
                          "omics_layer": omics, "assay": sub.assay_code, "source_id": sub.feature_ID,
                          "platform": sub.platform.fillna("")})
        if label == "PTMSite":
            d["protein_id"] = sub.feature_ID.str.rsplit("_", n=1).str[0].values
            d["ptm_type"] = ome
        if label == "Metabolite":
            d["refmet_name"] = np.where(sub.fid.str.startswith("REFMET:"), sub.fid.str[7:], "")
        g.node(label, d)

    ed = pd.DataFrame({
        "fid": sig.fid, "species": "rat", "tissue": sig.tissue, "group_id": "rat:endurance_training",
        "sex": sig.sex, "timepoint": sig.training_group,
        "timepoint_order": sig.training_group.map({w: i for i, w in enumerate(RAT_WEEKS)}),
        "timescale": "weeks of training", "assay": sig.assay_code, "ome": sig.assay,
        "logFC": sig.timewise_logFC, "z": sig.z, "p_value": sig.timewise_p_value, "adj_p_value": sig.adj_p_used,
        "adj_p_type": sig.adj_p_type, "timewise_adj_p": sig.timewise_adj_p,
        "training_q": sig.training_q, "is_metab": sig.assay == "METAB",
        "significance_basis": np.where(sig.significance_basis == "graphical_state",
                                       "training_q<0.05 + MoTrPAC per-week state call",
                                       "training_q<0.05 + timewise p<0.05 (no state call)")})
    add_reg_edges(g, ed)

    # feature -> gene (not for metabolites)
    nm = sig[sig.assay != "METAB"].drop_duplicates("fid")
    k = f2g[f2g.feature_ID.isin(nm.feature_ID)].dropna(subset=["ensembl_gene"])
    id2fid = nm.set_index("feature_ID").fid.to_dict()
    k = k.assign(fid=k.feature_ID.map(id2fid)).dropna(subset=["fid"])
    g.edge("MAPS_TO_GENE", pd.DataFrame({
        "src": k.fid, "dst": "ENSEMBL:" + k.ensembl_gene, "species": "rat",
        "annotation": k.custom_annotation.fillna(""), "distance_to_gene": k.relationship_to_gene,
        "weight:float": [ANNOT_WEIGHT.get(a, 1.0) if isinstance(a, str) else 1.0 for a in k.custom_annotation]}))
    gn = k.drop_duplicates("ensembl_gene")
    g.node("Gene", pd.DataFrame({"id": "ENSEMBL:" + gn.ensembl_gene, "symbol": gn.gene_symbol, "species": "rat",
                                 "entrez_id": gn.entrez_gene.astype("Int64").astype(str)}))
    ptm = nm[nm.assay.isin(["PHOSPHO", "ACETYL", "UBIQ"])]
    g.edge("SITE_ON", pd.DataFrame({"src": ptm.fid, "dst": "rat:PROT:" + ptm.feature_ID.str.rsplit("_", n=1).str[0],
                                    "species": "rat", "weight:float": 1.0}))
    return trf, sig


# ------------------------------------------------------------------------------------ human
def load_human_da(hum):
    parts = []
    for f in sorted(os.listdir(hum)):
        if f.endswith("_DA.tsv.gz"):
            d = rd(f"{hum}/{f}")
            if "platform" not in d:
                d["platform"] = d.assay
            parts.append(d[d.contrast_type == "exercise_with_controls"])
    da = pd.concat(parts, ignore_index=True)
    da["platform"] = da.platform.fillna(da.assay)
    return da


def add_human(g, hum):
    f2g = rd(f"{hum}/HUMAN_FEATURE_TO_GENE.tsv.gz")
    refmet = f2g[f2g.assay == "metab"].dropna(subset=["refmet_name"]).drop_duplicates("feature_id") \
        .set_index("feature_id").refmet_name.to_dict()
    da = load_human_da(hum)
    da["is_clinical"] = da.platform.isin(HUMAN_CLINICAL) | da.assay.isin(HUMAN_CLINICAL)
    da["fid"] = np.where(da.assay == "metab", "REFMET:" + da.feature_id.map(lambda x: refmet.get(x, x)),
                         "human:" + da.assay + ":" + da.feature_id)
    sig = da[(da.adj_p_value < 0.05) & ~da.is_clinical].copy()

    fe = sig.drop_duplicates("fid")
    for assay, sub in fe.groupby("assay"):
        label, omics = HUM_ASSAY[assay]
        d = pd.DataFrame({"id": sub.fid, "species": np.where(sub.fid.str.startswith("REFMET:"), "rat;human", "human"),
                          "omics_layer": omics, "assay": assay, "source_id": sub.feature_id,
                          "platform": sub.platform})
        if label == "PTMSite":
            d["protein_id"] = sub.feature_id.str.rsplit("_", n=1).str[0].values
            d["ptm_type"] = "PHOSPHO"
        if label == "Metabolite":
            d["refmet_name"] = sub.fid.str[7:].values
        g.node(label, d)

    ed = pd.DataFrame({
        "fid": sig.fid, "species": "human", "tissue": sig.tissue, "group_id": sig.randomGroupCode.map(HUMAN_GROUP_NODE),
        "sex": "both", "timepoint": sig.Timepoint,
        "timepoint_order": sig.Timepoint.map({t: i for i, t in enumerate(HUMAN_TP_ORDER)}),
        "timescale": "time relative to acute bout", "assay": sig.platform, "ome": sig.assay,
        "logFC": sig.logFC, "z": sig["z.std"], "p_value": sig.p_value, "adj_p_value": sig.adj_p_value,
        "is_metab": sig.assay == "metab", "significance_basis": "adj_p<0.05 (BH within contrast)",
        "adj_p_type": "adj p (BH within contrast)"})
    add_reg_edges(g, ed)

    nm = sig[sig.assay != "metab"].drop_duplicates("fid")
    m = f2g.merge(nm[["assay", "feature_id", "fid"]], on=["assay", "feature_id"]).dropna(subset=["ensembl_gene"])
    m["gene"] = "ENSEMBL:" + m.ensembl_gene.map(strip_ver)
    g.edge("MAPS_TO_GENE", pd.DataFrame({
        "src": m.fid, "dst": m.gene, "species": "human", "annotation": m.custom_annotation.fillna(""),
        "distance_to_gene": m.relationship_to_gene,
        "weight:float": [ANNOT_WEIGHT.get(a, 1.0) if isinstance(a, str) else 1.0 for a in m.custom_annotation]}))
    tfs = set(rd(f"{hum}/UTORONTO_TFs.tsv.gz").gene_symbol)
    gn = m.drop_duplicates("gene")
    g.node("Gene", pd.DataFrame({"id": gn.gene, "symbol": gn.gene_symbol, "species": "human",
                                 "entrez_id": gn.entrez_gene.astype("Int64").astype(str),
                                 "is_TF:boolean": gn.gene_symbol.isin(tfs)}))
    ph = nm[nm.assay == "prot-ph"]
    g.edge("SITE_ON", pd.DataFrame({"src": ph.fid, "dst": "human:prot-pr:" + ph.feature_id.str.rsplit("_", n=1).str[0],
                                    "species": "human", "weight:float": 1.0}))
    flank = f2g[f2g.assay == "prot-ph"].dropna(subset=["flanking_sequence"])[["feature_id", "flanking_sequence"]]
    return da, sig, flank


# ------------------------------------------------------------------------------------ orthology
def add_orthology(g, rat):
    o = rd(f"{rat}/RAT_TO_HUMAN_GENE.tsv.gz").dropna(subset=["RAT_ENSEMBL_ID", "HUMAN_ORTHOLOG_ENSEMBL_ID"])
    o = o.assign(RAT_ENSEMBL_ID=o.RAT_ENSEMBL_ID.str.split(";")).explode("RAT_ENSEMBL_ID")
    o["src"], o["dst"] = "ENSEMBL:" + o.RAT_ENSEMBL_ID, "ENSEMBL:" + o.HUMAN_ORTHOLOG_ENSEMBL_ID.map(strip_ver)
    genes = g.node_ids("Gene")
    o = o[o.src.isin(genes)].drop_duplicates(["src", "dst"])
    g.edge("ORTHOLOG_OF", pd.DataFrame({"src": o.src, "dst": o.dst, "source": o.HUMAN_ORTHOLOG_SOURCE,
                                        "weight:float": 1.0 / o.groupby("src").dst.transform("count")}))
    # human ortholog nodes (so rat genes reach human-symbol pathways even if not significant in human)
    g.node("Gene", pd.DataFrame({"id": o.dst, "symbol": o.HUMAN_ORTHOLOG_SYMBOL, "species": "human",
                                 "entrez_id": o.HUMAN_ORTHOLOG_NCBI_GENE_ID.astype(str),
                                 "significant_in_human:boolean": o.dst.isin(genes)}))
    ps = rd(f"{rat}/RAT_TO_HUMAN_PHOSPHO.tsv.gz").dropna()
    g.edge("ORTHOLOGOUS_SITE", pd.DataFrame({"src": "rat:PHOSPHO:" + ps.ptm_id_rat_refseq,
                                             "dst": "human:prot-ph:" + ps.ptm_id_human_uniprot, "weight:float": 1.0}))
    return o, ps


# ------------------------------------------------------------------------------------ pathways
def molecular_signatures(hum):
    cache = f"{hum}/MOLECULAR_SIGNATURES_long.tsv.gz"
    if not os.path.exists(cache):
        r = ('x<-readRDS("%s/MOLECULAR_SIGNATURES.rds"); out<-do.call(rbind, lapply(names(x), function(db) '
             'do.call(rbind, lapply(names(x[[db]]), function(s) data.frame(database=db, set=s, member=x[[db]][[s]]))))); '
             'con<-gzfile("%s","w"); write.table(out,con,sep="\\t",quote=FALSE,row.names=FALSE); close(con)') % (hum, cache)
        subprocess.run(["Rscript", "-e", r], check=True)
    return rd(cache)


def add_pathways(g, rat, hum, ortho, rat_ps, flank, rat_sig, hum_sig, go="enriched"):
    ms = molecular_signatures(hum)
    s2id = rd(f"{hum}/SET_TO_ID.tsv.gz")
    cam = rd(f"{hum}/CAMERA_RESULTS.tsv.gz")
    ptmsea = rd(f"{hum}/PTMSEA_RESULTS.tsv.gz")
    conv = rd(f"{hum}/CONTRAST_CONVERTER.tsv.gz")[["contrast", "randomGroupCode", "Timepoint"]]
    enriched_sets = set(cam[(cam.adj_p_value < 0.05) & (cam.contrast_type == "exercise_with_controls")].set) | \
        set(ptmsea[(ptmsea.adj_p_value < 0.05) & (ptmsea.contrast_type == "exercise_with_controls")].set)
    genes = pd.concat(g.nodes["Gene"]).drop_duplicates("id")
    hsym = genes[genes.species == "human"].set_index("symbol").id.to_dict()

    # --- human-curated gene sets (MSigDB-style); GO only if enriched in some contrast
    go_ok = ms.database.isin(GO) & (ms.set.isin(enriched_sets) if go == "enriched" else (go == "all"))
    gene_sets = ms[ms.database.isin(CURATED) | go_ok]
    gm = gene_sets[gene_sets.member.isin(hsym)].copy()
    gm["gene"] = gm.member.map(hsym)
    # human genes that are significant in human, plus human orthologs of significant rat genes
    rat_to_h = ortho.groupby("dst").src.apply(list).to_dict()
    sig_h = set(e for d in g.edges["MAPS_TO_GENE"] for e in d.dst if e.startswith("ENSEMBL:ENSG"))
    rows = []
    for r in gm.itertuples(index=False):
        if r.gene in sig_h:
            rows.append((r.gene, r.database + ":" + r.set, "human", "direct"))
        for rg in rat_to_h.get(r.gene, []):
            rows.append((rg, r.database + ":" + r.set, "rat", "via_human_ortholog"))
    mem = pd.DataFrame(rows, columns=["src", "dst", "species", "evidence"]).drop_duplicates()
    mem["weight:float"] = np.where(mem.evidence == "direct", 1.0, 0.8)
    g.edge("IN_PATHWAY", mem)
    used = gene_sets[(gene_sets.database + ":" + gene_sets.set).isin(set(mem.dst))]
    sz = used.groupby(["database", "set"]).size().reset_index(name="n")
    g.node("Pathway", pd.DataFrame({"id": sz.database + ":" + sz.set, "name": sz.set, "database": sz.database,
                                    "source_species": "human (MSigDB-curated)", "size:int": sz.n}))

    # --- PTM signatures (site level, via 15-mer flanking sequence)
    ptm_sets = ms[ms.database.isin(PTM_SETS)].copy()
    ptm_sets["flank"] = ptm_sets.member.str.split(";").str[0]
    ptm_sets["sig_dir"] = ptm_sets.member.str.split(";").str[1].fillna("")
    fl = flank.assign(fid="human:prot-ph:" + flank.feature_id)
    sig_sites = set(hum_sig[hum_sig.assay == "prot-ph"].fid)
    h_sites = fl[fl.fid.isin(sig_sites)]
    # rat sites reach the same signatures through their human orthologous site
    r2h = rat_ps.assign(r="rat:PHOSPHO:" + rat_ps.ptm_id_rat_refseq, h="human:prot-ph:" + rat_ps.ptm_id_human_uniprot)
    r_sites = r2h[r2h.r.isin(set(rat_sig.fid))].merge(fl, left_on="h", right_on="fid")
    pm = pd.concat([
        h_sites.merge(ptm_sets, left_on="flanking_sequence", right_on="flank").assign(src=lambda d: d.fid, species="human", evidence="direct"),
        r_sites.merge(ptm_sets, left_on="flanking_sequence", right_on="flank").assign(src=lambda d: d.r, species="rat", evidence="via_human_site")])
    if len(pm):
        g.edge("IN_PATHWAY", pd.DataFrame({"src": pm.src, "dst": pm.database + ":" + pm.set, "species": pm.species,
                                           "evidence": pm.evidence, "signature_direction": pm.sig_dir,
                                           "weight:float": np.where(pm.evidence == "direct", 1.0, 0.8)}).drop_duplicates())
        ps_sz = ptm_sets.groupby(["database", "set"]).size().reset_index(name="n")
        ps_sz = ps_sz[(ps_sz.database + ":" + ps_sz.set).isin(set(pm.database + ":" + pm.set))]
        g.node("Pathway", pd.DataFrame({"id": ps_sz.database + ":" + ps_sz.set, "name": ps_sz.set,
                                        "database": ps_sz.database, "source_species": "human (PTM signature)",
                                        "size:int": ps_sz.n}))

    # --- human pathway enrichment -> exercise group & organ
    enr = []
    for d, meth, zc in [(cam, "CAMERA-PR", "z.std"), (ptmsea, "PTM-SEA", "NES")]:
        d = d[(d.contrast_type == "exercise_with_controls") & (d.adj_p_value < 0.05)].merge(conv, on="contrast")
        d = d[d.set.isin(set(ms.set))]
        db = d.set.map(ms.drop_duplicates("set").set_index("set").database)
        enr.append(pd.DataFrame({"pid": db + ":" + d.set, "species": "human", "tissue": d.tissue,
                                 "group": d.randomGroupCode.map(HUMAN_GROUP_NODE), "sex": "both",
                                 "timepoint": d.Timepoint, "timescale": "time relative to acute bout",
                                 "assay": d.assay, "method": meth, "score": d[zc], "adj_p_value": d.adj_p_value,
                                 "up": d.direction.str.lower().eq("up")}))
    # --- rat pathway enrichment of graphical clusters (gene membership from the intersections)
    pw = rd(f"{rat}/GRAPH_PW_ENRICH.tsv.gz")
    pw = pw[pw.adj_p_value < 0.05].copy()
    pw["pid"] = pw.term_id
    g.node("Pathway", pd.DataFrame({"id": pw.term_id, "name": pw.term_name, "database": pw.source,
                                    "source_species": "rat (g:Profiler)", "size:int": pw.term_size}))
    inter = pw.assign(g=pw.intersection.str.split(",")).explode("g")
    g.edge("IN_PATHWAY", pd.DataFrame({"src": "ENSEMBL:" + inter.g, "dst": inter.term_id, "species": "rat",
                                       "evidence": "enrichment_intersection", "weight:float": 1.0}).drop_duplicates())
    # last week in the cluster name carries the state -> per-sex direction
    last = pw.cluster.str.split(r"->|---").str[-1].str.split(":").str[-1]
    wk = last.str.extract(r"^(\dw)_F(-?\d)_M(-?\d)$")
    span = pw.cluster.str.split(":").str[1].str.findall(r"(\dw)_").str.join("-")
    for sex, col in [("female", 1), ("male", 2)]:
        st = wk[col].astype(float)
        keep = st != 0
        enr.append(pd.DataFrame({"pid": pw.term_id[keep], "species": "rat", "tissue": pw.tissue[keep],
                                 "group": "rat:endurance_training", "sex": sex, "timepoint": span[keep],
                                 "timescale": "weeks of training", "assay": pw.ome[keep], "method": "g:Profiler",
                                 "score": np.nan, "adj_p_value": pw.adj_p_value[keep], "up": (st[keep] > 0)}))
    e = pd.concat(enr, ignore_index=True)
    e["type"] = np.where(e.up, "ENRICHED_UP_IN", "ENRICHED_DOWN_IN")
    for dst_col in ["group", "tissue_node"]:
        e["tissue_node"] = e.species + ":" + e.tissue
        for et, sub in e.groupby("type"):
            g.edge(et, pd.DataFrame({"src": sub.pid, "dst": sub[dst_col], "species": sub.species,
                                     "tissue": sub.tissue, "group": sub.group, "sex": sub.sex,
                                     "timepoint": sub.timepoint, "timescale": sub.timescale, "assay": sub.assay,
                                     "method": sub.method, "score:float": sub.score,
                                     "adj_p_value:float": sub.adj_p_value,
                                     "weight:float": w_from_p(sub.adj_p_value)}))


# ------------------------------------------------------------------------------------ interactions (PPI)
PPI_SOURCES = "Hetionet v1.0 Gene-interacts-Gene (HI-II-14, Lit-BM-13, iRefIndex, Incomplete Interactome); STRING v12 physical subnetwork"


def add_interactions(g, hum, ppi_dir, rat_ps, flank, rat_sig, hum_sig, string_min=700):
    """Human protein-protein interactions among graph genes, ligand-receptor pairs, kinase-substrate sites.
    All interaction evidence is human; rat genes reach it through ORTHOLOG_OF (flagged on each edge)."""
    genes = pd.concat(g.nodes["Gene"]).drop_duplicates("id")
    hum_g = genes[genes.species == "human"].copy()
    hum_g["ez"] = hum_g.entrez_id.astype(str).str.replace(r"\.0$", "", regex=True)
    ez2id = hum_g.dropna(subset=["ez"]).drop_duplicates("ez").set_index("ez").id.to_dict()
    sym2id = hum_g.drop_duplicates("symbol").set_index("symbol").id.to_dict()
    maps = pd.concat(g.edges["MAPS_TO_GENE"])
    sig_h = set(maps[maps.species == "human"].dst)
    orth = pd.concat(g.edges["ORTHOLOG_OF"])
    rat_sig_genes = set(maps[maps.species == "rat"].dst)
    h_with_rat = set(orth[orth.src.isin(rat_sig_genes)].dst)

    # 1. physical PPI: Hetionet (curated) + STRING v12 physical subnetwork (score >= string_min), merged per pair
    h = rd(f"{ppi_dir}/hetionet_edges.sif.gz")
    h = h[h.metaedge == "GiG"]
    het = pd.DataFrame({"src": h.source.str.split("::").str[1].map(ez2id), "dst": h.target.str.split("::").str[1].map(ez2id)}).dropna()
    het["hetionet"] = True
    parts = [het]
    sl, si = f"{ppi_dir}/9606.protein.physical.links.v12.0.txt.gz", f"{ppi_dir}/9606.protein.info.v12.0.txt.gz"
    if os.path.exists(sl) and os.path.exists(si):
        info = pd.read_csv(si, sep="\t", usecols=[0, 1], names=["sid", "sym"], header=0)
        s2g = info.assign(g=info.sym.map(sym2id)).dropna(subset=["g"]).set_index("sid").g.to_dict()
        st = pd.read_csv(sl, sep=" ")
        st = st[st.combined_score >= string_min]
        st = pd.DataFrame({"src": st.protein1.map(s2g), "dst": st.protein2.map(s2g), "string_score": st.combined_score}).dropna()
        parts.append(st)
    e = pd.concat(parts, ignore_index=True)
    e = e[e.src != e.dst]
    e[["src", "dst"]] = np.sort(e[["src", "dst"]].values, axis=1)
    e = e.groupby(["src", "dst"], as_index=False).agg(hetionet=("hetionet", lambda x: x.fillna(False).astype(bool).any()) if "hetionet" in e else ("src", "size"),
                                                      string_score=("string_score", "max") if "string_score" in e else ("src", "size"))
    has_s = e.string_score.notna() if "string_score" in e else False
    src_lab = np.where(e.hetionet & has_s, "Hetionet;STRING", np.where(e.hetionet, "Hetionet", "STRING"))
    g.edge("INTERACTS_WITH", pd.DataFrame({
        "src": e.src, "dst": e.dst, "interaction": "physical", "source": src_lab, "evidence_species": "human",
        "string_score:int": e.string_score.astype("Int64") if "string_score" in e else pd.NA,
        "both_significant_human:boolean": e.src.isin(sig_h) & e.dst.isin(sig_h),
        "both_significant_rat_orthologs:boolean": e.src.isin(h_with_rat) & e.dst.isin(h_with_rat),
        # weight: STRING confidence when present (0.7-1.0), curated-only Hetionet pairs 0.8
        "weight:float": np.where(has_s, e.string_score.fillna(800) / 1000.0, 0.8)}))

    # 2. ligand -> receptor (LIANA consensus, OmniPath-derived); complexes split into subunits
    l = pd.read_csv(f"{ppi_dir}/liana_omni_resource.csv")
    l = l[l.resource == "consensus"].assign(L=lambda d: d.source_genesymbol.str.split("_"),
                                            R=lambda d: d.target_genesymbol.str.split("_")).explode("L").explode("R")
    l = l[l.L.isin(sym2id) & l.R.isin(sym2id)]
    lr = pd.DataFrame({"src": l.L.map(sym2id), "dst": l.R.map(sym2id)}).drop_duplicates()
    g.edge("LIGAND_OF", pd.DataFrame({
        "src": lr.src, "dst": lr.dst, "interaction": "ligand-receptor", "source": "LIANA consensus (OmniPath)",
        "evidence_species": "human", "ligand_significant_human:boolean": lr.src.isin(sig_h),
        "receptor_significant_human:boolean": lr.dst.isin(sig_h), "weight:float": 1.0}))

    # 3. kinase -> substrate site (PhosphoSitePlus kinase sets bundled with the human package)
    ms = molecular_signatures(hum)
    psp = ms[ms.database == "PSP"].copy()
    psp["kinase"] = psp.set.str.replace(r"^PSP_", "", regex=True)
    psp["flank"] = psp.member.str.split(";").str[0]
    fl = flank.assign(fid="human:prot-ph:" + flank.feature_id)
    sig_sites = set(hum_sig[hum_sig.assay == "prot-ph"].fid)
    hs = fl[fl.fid.isin(sig_sites)].merge(psp, left_on="flanking_sequence", right_on="flank").assign(site=lambda d: d.fid, species="human", evidence="direct")
    r2h = rat_ps.assign(r="rat:PHOSPHO:" + rat_ps.ptm_id_rat_refseq, h="human:prot-ph:" + rat_ps.ptm_id_human_uniprot)
    rs = r2h[r2h.r.isin(set(rat_sig.fid))].merge(fl, left_on="h", right_on="fid").merge(
        psp, left_on="flanking_sequence", right_on="flank").assign(site=lambda d: d.r, species="rat", evidence="via_human_orthologous_site")
    ks = pd.concat([hs, rs], ignore_index=True)
    # kinases need a Gene node even when the kinase itself does not change
    f2g = rd(f"{hum}/HUMAN_FEATURE_TO_GENE.tsv.gz").dropna(subset=["gene_symbol", "ensembl_gene"])
    s2e = f2g.drop_duplicates("gene_symbol").set_index("gene_symbol").ensembl_gene.map(strip_ver).to_dict()
    new = sorted(set(ks.kinase) - set(sym2id))
    g.node("Gene", pd.DataFrame({"id": [f"ENSEMBL:{s2e[k]}" if k in s2e else f"SYMBOL:{k}" for k in new], "symbol": new,
                                 "species": "human", "is_kinase:boolean": True}))
    for k in new:
        sym2id[k] = f"ENSEMBL:{s2e[k]}" if k in s2e else f"SYMBOL:{k}"
    g.edge("PHOSPHORYLATES", pd.DataFrame({
        "src": ks.kinase.map(sym2id), "dst": ks.site, "interaction": "kinase-substrate", "source": "PhosphoSitePlus (MSigDB bundle)",
        "site_species": ks.species, "evidence": ks.evidence,
        "weight:float": np.where(ks.evidence == "direct", 1.0, 0.8)}).drop_duplicates(["src", "dst"]))


# ------------------------------------------------------------------------------------ metabolite co-regulation
def profile_matrix(df, key, cols, val):
    return df.pivot_table(index=key, columns=cols, values=val, aggfunc="mean")


def corr_edges(M, P, r_min, top_k, min_n, meta):
    """Pearson r between rows of M (metabolites) and rows of P (partners), pairwise-complete."""
    out = []
    if M.empty or P.empty:
        return out
    common = M.columns.intersection(P.columns)
    M, P = M[common], P[common]
    for mid, mv in M.iterrows():
        ok = mv.notna().values
        if ok.sum() < min_n:
            continue
        Pm = P.loc[:, ok]
        Pm = Pm[Pm.notna().sum(axis=1) == ok.sum()]
        if Pm.empty:
            continue
        x = mv.values[ok].astype(float)
        X = Pm.values.astype(float)
        xc = x - x.mean()
        Xc = X - X.mean(axis=1, keepdims=True)
        den = np.sqrt((xc ** 2).sum() * (Xc ** 2).sum(axis=1))
        with np.errstate(invalid="ignore", divide="ignore"):
            r = (Xc @ xc) / den
        idx = np.where(np.abs(r) >= r_min)[0]
        idx = idx[np.argsort(-np.abs(r[idx]))][:top_k]
        for i in idx:
            pid = Pm.index[i]
            if pid == mid:
                continue
            out.append((mid, pid, float(r[i]), int(ok.sum()), *meta))
    return out


def add_metabolite_coregulation(g, rat_trf, rat_sig, hum_da, hum_sig, r_min, top_k):
    rows = []
    # rat: profiles = logFC over sex x week (8 points), within tissue
    rt = rat_trf.copy()
    rt["col"] = rt.sex + "_" + rt.training_group
    sig_ids = set(rat_sig.fid)
    rt = rt[rt.fid.isin(sig_ids)]
    for tissue, sub in rt.groupby("tissue"):
        met = sub[sub.assay == "METAB"]
        prot = sub[sub.assay.isin(["PROT", "IMMUNO"])]
        M = profile_matrix(met, "fid", "col", "timewise_logFC")
        meta = ("rat", tissue, "logFC over sex x week (8)")
        rows += [(*r, "metabolite") for r in corr_edges(M, M, r_min, top_k, 6, meta)]
        rows += [(*r, "protein") for r in corr_edges(M, profile_matrix(prot, "fid", "col", "timewise_logFC"),
                                                     r_min, top_k, 6, meta)]
    # human: profiles = logFC over group x timepoint (all control-adjusted contrasts), within tissue
    hd = hum_da[hum_da.fid.isin(set(hum_sig.fid)) & ~hum_da.is_clinical].copy()
    hd["col"] = hd.randomGroupCode + "_" + hd.Timepoint
    for tissue, sub in hd.groupby("tissue"):
        met = sub[sub.assay == "metab"]
        prot = sub[sub.assay.isin(["prot-pr", "prot-ol"])]
        M = profile_matrix(met, "fid", "col", "logFC")
        n = M.shape[1]
        meta = ("human", tissue, f"logFC over group x timepoint ({n})")
        rows += [(*r, "metabolite") for r in corr_edges(M, M, r_min, top_k, min(6, n), meta)]
        rows += [(*r, "protein") for r in corr_edges(M, profile_matrix(prot, "fid", "col", "logFC"),
                                                     r_min, top_k, min(6, n), meta)]
    e = pd.DataFrame(rows, columns=["src", "dst", "r", "n_points", "species", "tissue", "profile", "partner"])
    # metabolite-metabolite pairs are symmetric: keep one orientation
    mm = e.partner == "metabolite"
    a, b = e.src.where(~mm | (e.src < e.dst), e.dst), e.dst.where(~mm | (e.src < e.dst), e.src)
    e["src"], e["dst"] = a, b
    e = e.drop_duplicates(["src", "dst", "species", "tissue"])
    g.edge("CO_REGULATED_WITH", pd.DataFrame({
        "src": e.src, "dst": e.dst, "species": e.species, "tissue": e.tissue, "partner_type": e.partner,
        "pearson_r:float": e.r, "n_points:int": e.n_points, "profile": e.profile,
        "sign": np.where(e.r > 0, "positive", "negative"), "weight:float": e.r.abs()}))
    return e


# ------------------------------------------------------------------------------------ phenotypes
RAT_PHENO = {
    "calculated.variables.vo2_max_change": ("VO2max change", "%", None),
    "calculated.variables.pct_body_fat_change": ("Body fat change", "% points", None),
    "calculated.variables.pct_body_lean_change": ("Lean mass change", "% points", None),
    "calculated.variables.lactate_change_dueto_train": ("Blood lactate change (VO2max test)", "mmol/L", "BLOOD"),
    "terminal.weight.bw": ("Body weight at sacrifice", "g", None),
    "terminal.weight.lg": ("Lateral gastrocnemius mass", "mg", "SKM-GN"),
    "terminal.weight.mg": ("Medial gastrocnemius mass", "mg", "SKM-GN"),
    "terminal.weight.pl": ("Plantaris mass", "mg", None),
    "terminal.weight.sol": ("Soleus mass", "mg", None),
}
RAT_GROUP_WEEK = {"One-week program": "1w", "Two-week program": "2w", "Four-week program": "4w",
                  "Eight-week program Training Group": "8w"}


def add_phenotypes(g, rat, hum_da):
    p = rd(f"{rat}/PHENO.tsv.gz").drop_duplicates("pid")
    ctrl = p[p["key.anirandgroup"] == "Eight-week program Control Group"]
    nodes, rows = [], []
    for col, (name, unit, tissue) in RAT_PHENO.items():
        pid = "rat:pheno:" + col.split(".")[-1]
        nodes.append((pid, name, unit, "rat", "physiology"))
        if tissue:
            g.edge("MEASURED_IN", pd.DataFrame({"src": [pid], "dst": [f"rat:{tissue}"], "weight:float": [1.0]}))
        for grp, wk in RAT_GROUP_WEEK.items():
            for sex in ["female", "male"]:
                a = pd.to_numeric(p[(p["key.anirandgroup"] == grp) & (p.sex == sex)][col], errors="coerce").dropna()
                c = pd.to_numeric(ctrl[ctrl.sex == sex][col], errors="coerce").dropna()
                if len(a) < 3 or len(c) < 3:
                    continue
                t, pv = stats.ttest_ind(a, c, equal_var=False)
                rows.append((pid, sex, wk, a.mean(), c.mean(), a.mean() - c.mean(), t, pv, len(a), len(c)))
    g.node("Phenotype", pd.DataFrame(nodes, columns=["id", "name", "unit", "species", "category"]))
    r = pd.DataFrame(rows, columns=["src", "sex", "timepoint", "trained_mean", "control_mean", "diff", "t", "p", "n_trained", "n_control"])
    r["adj_p"] = stats.false_discovery_control(r.p) if len(r) else []
    r = r[r.adj_p < 0.05].copy()   # BH across all rat phenotype tests
    r["w"] = (-np.log10(r.adj_p)).rank(pct=True).round(4)
    for et, sub in r.groupby(np.where(r["diff"] > 0, "INCREASED_IN", "DECREASED_IN")):
        g.edge(et, pd.DataFrame({
            "src": sub.src, "dst": "rat:endurance_training", "species": "rat", "tissue": "whole-body", "sex": sub.sex,
            "timepoint": sub.timepoint, "timescale": "weeks of training", "assay": "physiology",
            "trained_mean:float": sub.trained_mean, "control_mean:float": sub.control_mean,
            "difference:float": sub["diff"], "t:float": sub.t, "p_value:float": sub.p, "adj_p_value:float": sub.adj_p,
            "significance_basis": "Welch t-test vs sedentary control, BH adj_p<0.05",
            "weight:float": sub.w}))
    # human clinical chemistry = phenotypes
    cl = hum_da[hum_da.platform.isin(HUMAN_CLINICAL) | hum_da.assay.isin(HUMAN_CLINICAL)].copy()
    cl["pid"] = "human:pheno:" + cl.feature_id
    nd = cl.drop_duplicates("pid")
    g.node("Phenotype", pd.DataFrame({"id": nd.pid, "name": nd.feature_id, "unit": "log2 (relative)",
                                      "species": "human", "category": "clinical chemistry (" + nd.platform + ")"}))
    g.edge("MEASURED_IN", pd.DataFrame({"src": nd.pid, "dst": "human:" + nd.tissue, "weight:float": 1.0}))
    s = cl[cl.adj_p_value < 0.05].copy()
    s["w"] = (-np.log10(s.adj_p_value)).rank(pct=True).round(4)
    for et, sub in s.groupby(np.where(s.logFC > 0, "INCREASED_IN", "DECREASED_IN")):
        g.edge(et, pd.DataFrame({
            "src": sub.pid, "dst": sub.randomGroupCode.map(HUMAN_GROUP_NODE), "species": "human", "tissue": sub.tissue, "sex": "both",
            "timepoint": sub.Timepoint, "timescale": "time relative to acute bout", "assay": sub.platform,
            "z:float": sub["z.std"], "p_value:float": sub.p_value,
            "adj_p_value:float": sub.adj_p_value, "neg_log10_adj_p:float": (-np.log10(sub.adj_p_value)).round(3),
            "significance_basis": "adj_p<0.05 (BH within contrast)",
            "weight:float": sub.w}))


# ------------------------------------------------------------------------------------ citations
HUMAN_TISSUE_REFS = {"muscle": "D3;D4;D8;M2;M1", "adipose": "D3;D5;D8;M2;M1", "blood": "D3;D6;D8;M2;M1"}
PATHWAY_DB_REFS = {"REACTOME": "R9;R8;D8", "WP": "R10;R8;D8", "KEGG_MEDICUS": "R11;R8;D8", "PID": "R12;R8;D8",
                   "BIOCARTA": "R8;D8", "MITOCARTA": "R13;R8;D8", "GOBP": "R14;R8;D8", "GOCC": "R14;R8;D8",
                   "GOMF": "R14;R8;D8", "PSP": "R6;R8;D8", "PTMSIGDB": "R7;R8;D8", "KEGG": "R11;D1;M4", "REAC": "R9;D1;M4"}


def cite_edges(g):
    """Attach `source_refs` (ids in refs/references.json) to every edge row."""
    for et, parts in g.edges.items():
        for d in parts:
            sp = d["species"] if "species" in d else pd.Series("", index=d.index)
            if et in ("UPREGULATED_IN", "DOWNREGULATED_IN"):
                r = np.where(sp == "rat", "D1;D2", d["tissue"].map(HUMAN_TISSUE_REFS).fillna("D3;D8;M1"))
            elif et in ("ENRICHED_UP_IN", "ENRICHED_DOWN_IN"):
                r = np.where(sp == "rat", "D1;D2;M4", np.where(d.get("method", "") == "PTM-SEA", "D8;R7", "D8;M3"))
            elif et == "IN_PATHWAY":
                db = d["dst"].str.split(":").str[0]
                r = db.map(PATHWAY_DB_REFS).fillna("R8;D8") + np.where(d.get("evidence", "").astype(str).str.startswith("via"), ";D2", "")
            elif et == "INTERACTS_WITH":
                r = d["source"].map({"STRING": "R1", "Hetionet": "R2;R3", "Hetionet;STRING": "R1;R2;R3"})
            elif et == "CO_REGULATED_WITH":
                r = np.where(sp == "rat", "P1;D1;D2", "P1;D3;D8")
            elif et in ("INCREASED_IN", "DECREASED_IN"):
                r = np.where(sp == "rat", "D1;D2;M5;M1;P1", "D3;D8;M1")
            elif et == "PHOSPHORYLATES":
                r = np.where(d.get("site_species", "") == "rat", "R6;D8;D2", "R6;D8")
            elif et in ("MAPS_TO_GENE", "SITE_ON"):
                r = np.where(sp == "rat", "D2", "D8")
            else:
                r = {"ORTHOLOG_OF": "D2", "ORTHOLOGOUS_SITE": "D2", "LIGAND_OF": "R4;R5", "EQUIVALENT_TISSUE": "R16",
                     "MEASURED_IN": "D2;D8"}.get(et, "P1")
            d["source_refs"] = r


def rank_edges(g):
    """Per edge type: rank of each edge among the edges leaving its source (rank_from_src) and entering its
    target (rank_at_dst), strongest weight first. Filter r.rank_from_src <= 5 AND r.rank_at_dst <= 5 for a 5-edge view."""
    for et in list(g.edges):
        d = pd.concat(g.edges[et], ignore_index=True)
        w = [c for c in d.columns if c.startswith("weight")]
        wt = d[w[0]].astype(float).fillna(0) if w else pd.Series(1.0, index=d.index)
        d = d.dropna(subset=["src", "dst"]).drop_duplicates().copy(); wt = wt.loc[d.index]
        d["rank_from_src:int"] = wt.groupby(d["src"].astype(str)).rank(ascending=False, method="first").astype("Int64")
        d["rank_at_dst:int"] = wt.groupby(d["dst"].astype(str)).rank(ascending=False, method="first").astype("Int64")
        g.edges[et] = [d]


def add_reference_nodes(g, path):
    import json
    refs = json.load(open(path))
    g.node("Reference", pd.DataFrame([{"id": x["id"], "type": x["type"], "short": x["short"], "citation": x["full"],
                                       "doi": x.get("doi", ""), "url": x.get("url") or (("https://doi.org/" + x["doi"]) if x.get("doi") else ""),
                                       "used_for": x["use"]} for x in refs]))

# ------------------------------------------------------------------------------------ checks
def check_rules(out):
    lab = {}
    import glob
    for f in glob.glob(f"{out}/nodes_*.csv.gz"):
        d = pd.read_csv(f, usecols=["id:ID", ":LABEL"])
        lab.update(dict(zip(d["id:ID"], d[":LABEL"])))
    is_met = lambda s: s.map(lambda x: "Metabolite" in str(lab.get(x, "")))
    allowed = {"Metabolite", "ProteinFeature", "AffinityProtein", "Tissue"}
    bad = 0
    for f in glob.glob(f"{out}/edges_*.csv.gz"):
        d = pd.read_csv(f, usecols=[":START_ID", ":END_ID"])
        for a, b in [(":START_ID", ":END_ID"), (":END_ID", ":START_ID")]:
            m = d[is_met(d[a])]
            other = m[b].map(lambda x: str(lab.get(x, "")).split(";")[-1])
            bad += (~other.isin(allowed)).sum()
    assert bad == 0, f"{bad} metabolite edges violate the metabolite rules"
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rat", required=True); ap.add_argument("--human", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--r-min", type=float, default=0.9); ap.add_argument("--top-k", type=int, default=10)
    ap.add_argument("--ppi-dir", default=None, help="folder with hetionet_edges.sif.gz and liana_omni_resource.csv")
    ap.add_argument("--string-min", type=int, default=700, help="minimum STRING combined score (700 = high confidence)")
    ap.add_argument("--go", choices=["enriched", "all", "none"], default="enriched",
                    help="GO terms: only those enriched in some human contrast (default), all, or none")
    a = ap.parse_args()
    g = Graph()
    add_static(g)
    rat_trf, rat_sig = add_rat(g, a.rat)
    hum_da, hum_sig, flank = add_human(g, a.human)
    ortho, rat_ps = add_orthology(g, a.rat)
    add_pathways(g, a.rat, a.human, ortho, rat_ps, flank, rat_sig, hum_sig, a.go)
    if a.ppi_dir:
        add_interactions(g, a.human, a.ppi_dir, rat_ps, flank, rat_sig, hum_sig, a.string_min)
    add_metabolite_coregulation(g, rat_trf, rat_sig, hum_da, hum_sig, a.r_min, a.top_k)
    add_phenotypes(g, a.rat, hum_da)
    BASIS = {"MAPS_TO_GENE": "annotation (feature-to-gene map)", "SITE_ON": "annotation", "ORTHOLOG_OF": "curated orthology",
             "ORTHOLOGOUS_SITE": "curated site orthology", "IN_PATHWAY": "curated membership",
             "ENRICHED_UP_IN": "adj_p<0.05 (enrichment)", "ENRICHED_DOWN_IN": "adj_p<0.05 (enrichment)",
             "CO_REGULATED_WITH": "correlation |r|>=%s, top %d per metabolite (no significance test)" % (a.r_min, a.top_k),
             "LIGAND_OF": "curated (LIANA consensus)", "PHOSPHORYLATES": "curated (PhosphoSitePlus)",
             "EQUIVALENT_TISSUE": "anatomical mapping", "MEASURED_IN": "annotation"}
    for et, parts in g.edges.items():
        for d in parts:
            if "significance_basis" not in d.columns:
                if et == "INTERACTS_WITH":
                    d["significance_basis"] = np.where(d.source.str.contains("STRING"), "STRING physical combined_score>=%d" % a.string_min, "curated (Hetionet)")
                else:
                    d["significance_basis"] = BASIS.get(et, "")
    cite_edges(g)
    rank_edges(g)
    add_reference_nodes(g, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "refs", "references.json"))
    man = g.write(a.out)
    check_rules(a.out)
    man["metabolite_rules_check"] = "passed"
    json.dump(man, open(f"{a.out}/manifest.json", "w"), indent=2)
    print(json.dumps(man, indent=2))


if __name__ == "__main__":
    main()
