#!/usr/bin/env python3
"""Checks export CSVs against neo4j-admin bulk-import rules without needing Neo4j:
headers, typed columns parse, IDs unique across node files, every edge endpoint exists."""
import glob, sys, pandas as pd, numpy as np
D = sys.argv[1] if len(sys.argv) > 1 else "../exports_sig"
TYPES = {"int": "Int64", "float": float, "boolean": "bool"}
ok, ids, problems = True, {}, []
for f in sorted(glob.glob(f"{D}/nodes_*.csv.gz")):
    d = pd.read_csv(f, dtype=str, keep_default_na=False)
    if d.columns[0] != "id:ID" or ":LABEL" not in d.columns: problems.append(f"{f}: missing id:ID/:LABEL")
    for c in d.columns:
        if ":" in c and not c.startswith(":") and c != "id:ID":
            t = c.split(":")[1]; v = d[c][d[c] != ""]
            if t in ("int", "float"): bad = pd.to_numeric(v, errors="coerce").isna().sum()
            elif t == "boolean": bad = (~v.str.lower().isin(["true", "false"])).sum()
            else: bad = 0
            if bad: problems.append(f"{f}:{c} {bad} unparsable")
    for i in d["id:ID"]:
        if i in ids: problems.append(f"duplicate id {i} in {f} and {ids[i]}"); break
        ids[i] = f
    if d["id:ID"].eq("").any(): problems.append(f"{f}: empty ids")
n_edges = 0
for f in sorted(glob.glob(f"{D}/edges_*.csv.gz")):
    d = pd.read_csv(f, dtype=str, keep_default_na=False)
    for c in (":START_ID", ":END_ID", ":TYPE"):
        if c not in d.columns: problems.append(f"{f}: missing {c}")
    miss = (~d[":START_ID"].isin(ids)).sum() + (~d[":END_ID"].isin(ids)).sum()
    if miss: problems.append(f"{f}: {miss} dangling endpoints")
    for c in d.columns:
        if ":" in c and not c.startswith(":"):
            t = c.split(":")[1]; v = d[c][d[c] != ""]
            if t in ("int", "float") and pd.to_numeric(v, errors="coerce").isna().sum(): problems.append(f"{f}:{c} unparsable")
    n_edges += len(d)
print(f"{len(ids):,} unique node ids, {n_edges:,} edges checked in {D}")
print("\n".join(problems) if problems else "OK: files satisfy neo4j-admin bulk-import rules")
sys.exit(1 if problems else 0)
