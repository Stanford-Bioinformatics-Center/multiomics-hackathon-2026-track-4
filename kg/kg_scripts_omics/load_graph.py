#!/usr/bin/env python3
"""
load_graph.py — load the exported CSVs into NetworkX and run an example cross-species query.

  python load_graph.py --exports ../kg/exports

For Neo4j instead, use kg/neo4j/import.sh (one --nodes / --relationships flag per file;
a comma-separated list would treat the first file as the only header).
"""
import argparse, glob, os
import networkx as nx
import pandas as pd


def clean(col):
    return col.split(":")[0] if not col.startswith(":") else col


def load(exports):
    G = nx.MultiDiGraph()
    for f in sorted(glob.glob(f"{exports}/nodes_*.csv.gz")):
        d = pd.read_csv(f, low_memory=False)
        d.columns = [clean(c) for c in d.columns]
        for r in d.to_dict("records"):
            G.add_node(r.pop("id"), label=r.pop(":LABEL").split(";")[-1], **r)
    for f in sorted(glob.glob(f"{exports}/edges_*.csv.gz")):
        d = pd.read_csv(f, low_memory=False)
        d.columns = [clean(c) for c in d.columns]
        for r in d.to_dict("records"):
            G.add_edge(r.pop(":START_ID"), r.pop(":END_ID"), key=r.pop(":TYPE"), **r)
    return G


def conserved_responders(G, rat_tissue="rat:SKM-GN", human_tissue="human:muscle", top=20):
    """Genes whose rat product responds to training in `rat_tissue` AND whose human ortholog
    responds to an acute bout in `human_tissue`. Score = best rat weight x best human weight."""
    def best_response(gene, species, tissue):
        best = (0.0, None, None)
        for feat, _, k, e in G.in_edges(gene, keys=True, data=True):
            if k != "MAPS_TO_GENE":
                continue
            for _, con, k2, r in G.out_edges(feat, keys=True, data=True):
                if k2 != "RESPONDS_IN" or not con.startswith(species):
                    continue
                if G.nodes[con].get("tissue") and f"{species}:{G.nodes[con]['tissue']}" != tissue:
                    continue
                w = r["weight"] * e["weight"]          # response strength x feature->gene confidence
                if w > best[0]:
                    best = (w, con, r["direction"])
        return best

    rows = []
    for rg, hg, k, o in G.edges(keys=True, data=True):
        if k != "ORTHOLOG_OF":
            continue
        rw, rc, rd = best_response(rg, "rat", rat_tissue)
        if not rw:
            continue
        hw, hc, hd = best_response(hg, "human", human_tissue)
        if not hw:
            continue
        rows.append({"rat_gene": G.nodes[rg].get("symbol"), "human_gene": G.nodes[hg].get("symbol"),
                     "rat_best": rc, "rat_dir": rd, "rat_w": round(rw, 3),
                     "human_best": hc, "human_dir": hd, "human_w": round(hw, 3),
                     "same_direction": rd == hd, "score": round(rw * hw * o["weight"], 3)})
    return pd.DataFrame(rows).sort_values("score", ascending=False).head(top)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--exports", default=os.path.join(os.path.dirname(__file__), "..", "kg", "exports"))
    a = ap.parse_args()
    G = load(a.exports)
    print(f"{G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")
    print(conserved_responders(G).to_string(index=False))
