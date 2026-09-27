#!/usr/bin/env python3
"""How the knowledge graph's gene <-> metabolite edges are built from Human-GEM.

Human-GEM (SysBioChalmers, https://github.com/SysBioChalmers/Human-GEM) is a curated genome-scale model of
human metabolism. It has no gene-metabolite table; the link runs through reactions:

    gene --(reaction gene rule)--> reaction --(stoichiometric coefficient)--> metabolite

    * Each reaction carries a Boolean gene rule over Ensembl ids, e.g. "ENSG_A or (ENSG_B and ENSG_C)":
      "or" joins isozymes, "and" joins subunits of one complex.
    * Each reaction lists its metabolites with a stoichiometric coefficient: negative = consumed
      (substrate), positive = produced (product). A reaction with lower bound < 0 is reversible.
    * Metabolite ids are compartment-specific (MAM01796c = L-lactate in the cytosol, MAM01796e = the same
      molecule outside the cell). The graph collapses compartments: MAM01796.

Edge rule. A gene G and a metabolite M are linked when at least one reaction
    (1) has G anywhere in its gene rule (isozyme or complex subunit alike), and
    (2) has M, in any compartment, among its metabolites, and
    (3) M is not a "currency" metabolite: it appears in <= 60 reactions (ATP, H2O, H+, NAD(H), CoA, ...
        would otherwise link to almost every gene).
Edge attributes, aggregated over all reactions that satisfy (1)-(3):
    n_reactions   number of such reactions
    role          substrate (M consumed), product (M produced), or either (reversible reaction, or
                  different roles in different reactions)
    is_transport  M crosses a compartment boundary in at least one of the reactions (it appears in two
                  compartments of the same reaction), or the reaction is an exchange with the
                  extracellular space and its subsystem is a transport subsystem
    reactions     the reaction ids (first 20, sorted)
    subsystems    the reactions' subsystems (first 6, sorted)

Usage
    python human_gem_edges.py download                   # fetch Human-GEM v2.0.1 (~12 MB) into ./data/human_gem
    python human_gem_edges.py explain LDHA L-lactate      # show, step by step, how one gene-metabolite edge is computed
    python human_gem_edges.py explain SLC16A1 HMDB0000190 # metabolite by name, MAM id, KEGG, HMDB or ChEBI id
    python human_gem_edges.py build                       # write every edge to ./data/human_gem/gene_metabolite_edges.csv

`explain` and `build` download the model first if it is missing. Requires Python 3.9+ and PyYAML.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

import yaml

GEM_VERSION = "v2.0.1"
RAW_URL = "https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/{version}/model/{name}"
FILES = ("Human-GEM.yml", "metabolites.tsv", "genes.tsv")
DEFAULT_DIR = Path(__file__).resolve().parent / "data" / "human_gem"

CURRENCY_DEGREE = 60      # metabolites in more reactions than this get no gene edges
MAX_LISTED_REACTIONS = 20
MAX_LISTED_SUBSYSTEMS = 6


# ---------------------------------------------------------------- step 0: download

def download(gem_dir: Path, version: str = GEM_VERSION, force: bool = False) -> None:
    """Fetch the model YAML (reactions, metabolites, gene rules) and the metabolite / gene annotation tables."""
    gem_dir.mkdir(parents=True, exist_ok=True)
    for name in FILES:
        dest = gem_dir / name
        if dest.exists() and not force:
            print(f"  have {dest}")
            continue
        url = RAW_URL.format(version=version, name=name)
        print(f"  get  {url}")
        tmp = dest.with_suffix(dest.suffix + ".part")
        with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as fh:
            fh.write(r.read())
        tmp.replace(dest)


def ensure_downloaded(gem_dir: Path) -> None:
    if not all((gem_dir / n).exists() for n in FILES):
        print(f"Human-GEM not found in {gem_dir}; downloading {GEM_VERSION}")
        download(gem_dir)


# ---------------------------------------------------------------- step 1: load the model

class _Loader(getattr(yaml, "CSafeLoader", yaml.SafeLoader)):
    """YAML 1.1 reads a bare NO as boolean False, and NO is nitric oxide's name in Human-GEM: keep such words as text."""


_Loader.yaml_implicit_resolvers = {k: [(tag, rx) for tag, rx in v if tag != "tag:yaml.org,2002:bool"]
                                   for k, v in _Loader.yaml_implicit_resolvers.items()}


def _omap(x):
    """PyYAML reads `!!omap` as a list of (key, value) pairs."""
    return dict(x) if isinstance(x, list) else (x or {})


class HumanGEM:
    def __init__(self, gem_dir: Path):
        with open(gem_dir / "Human-GEM.yml") as fh:
            model = _omap(yaml.load(fh, Loader=_Loader))
        self.version = _omap(model["metaData"]).get("version")

        # metabolites: compartment-specific id -> name, compartment
        self.met_name, self.met_comp = {}, {}
        for m in model["metabolites"]:
            m = _omap(m)
            self.met_name[m["id"]] = m.get("name", "")
            self.met_comp[m["id"]] = m.get("compartment", "")

        # reactions: id -> name, {metabolite: coefficient}, gene rule, subsystems, reversibility
        self.rxns = {}
        for r in model["reactions"]:
            r = _omap(r)
            sub = r.get("subsystem") or []
            self.rxns[r["id"]] = dict(
                name=r.get("name", ""),
                mets={k: float(v) for k, v in _omap(r.get("metabolites")).items()},
                rule=(r.get("gene_reaction_rule") or "").strip(),
                subsystems=[sub] if isinstance(sub, str) else list(sub),
                reversible=float(r.get("lower_bound", 0.0)) < 0,
            )

        # genes: Ensembl id <-> symbol
        self.symbol, self.ensembl_of = {}, {}
        with open(gem_dir / "genes.tsv", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                self.symbol[row["genes"]] = row["geneSymbols"]
                if row["geneSymbols"]:
                    self.ensembl_of.setdefault(row["geneSymbols"].upper(), row["genes"])

        # metabolite cross-references (compartment-free id), used to resolve user input
        self.xref = {}
        with open(gem_dir / "metabolites.tsv", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                base = row["metsNoComp"]
                for col in ("metKEGGID", "metHMDBID", "metChEBIID", "metPubChemID"):
                    for v in (row.get(col) or "").split(";"):
                        v = v.strip().upper()
                        if v:
                            self.xref.setdefault(v, base)
                            if v.startswith("CHEBI:"):
                                self.xref.setdefault(v[6:], base)

        # reaction degree of each compartment-free metabolite: the currency test
        self.degree = defaultdict(set)
        for rid, r in self.rxns.items():
            for m in r["mets"]:
                self.degree[base_id(m)].add(rid)
        self.degree = {k: len(v) for k, v in self.degree.items()}

        # compartment-free id -> name (first compartment's name; names do not differ across compartments)
        self.base_name = {}
        for m, n in self.met_name.items():
            self.base_name.setdefault(base_id(m), n)

    # ---- input resolution
    def resolve_gene(self, q: str) -> str | None:
        q = q.strip()
        if re.fullmatch(r"ENSG\d+", q.upper()):
            return q.upper()
        return self.ensembl_of.get(q.upper())

    def resolve_metabolite(self, q: str) -> str | None:
        s = q.strip()
        if re.fullmatch(r"MAM\d+[a-z]?", s):
            return base_id(s) if s[-1].isalpha() else s
        if s.upper() in self.xref:
            return self.xref[s.upper()]
        by_name = {n.lower(): b for b, n in self.base_name.items()}
        return by_name.get(s.lower())

    def name_suggestions(self, q: str, k: int = 8) -> list[str]:
        q = q.lower()
        return sorted({n for n in self.base_name.values() if q in n.lower()}, key=len)[:k]


def base_id(met: str) -> str:
    """MAM01796c -> MAM01796 (Human-GEM appends a one-letter compartment code)."""
    return met[:-1]


def rule_genes(rule: str) -> set[str]:
    return set(re.findall(r"ENSG\d+", rule))


# ---------------------------------------------------------------- step 2: one reaction -> (metabolite, gene) contributions

def reaction_contributions(gem: HumanGEM, rid: str):
    """Yield (base metabolite, gene, role, is_transport) for every non-currency metabolite x rule gene of one reaction."""
    r = gem.rxns[rid]
    genes = rule_genes(r["rule"])
    if not genes:
        return
    kept = {m: c for m, c in r["mets"].items() if gem.degree[base_id(m)] <= CURRENCY_DEGREE}
    rxn_comps = {m[-1] for m in kept}
    comps_of = defaultdict(set)
    for m in kept:
        comps_of[base_id(m)].add(m[-1])
    subsystem_text = ";".join(r["subsystems"])
    for m, coef in kept.items():
        b = base_id(m)
        role = "either" if r["reversible"] else ("substrate" if coef < 0 else "product")
        crosses = len(comps_of[b]) > 1
        exchange = "e" in rxn_comps and len(rxn_comps) > 1 and "ransport" in subsystem_text
        for g in genes:
            yield b, g, role, crosses or exchange


# ---------------------------------------------------------------- step 3: aggregate over reactions -> edges

def build_edges(gem: HumanGEM) -> list[dict]:
    acc = defaultdict(lambda: dict(rxns=set(), roles=set(), transport=False, subsystems=set()))
    for rid in gem.rxns:
        for b, g, role, tr in reaction_contributions(gem, rid):
            e = acc[(b, g)]
            e["rxns"].add(rid); e["roles"].add(role); e["transport"] |= tr
            e["subsystems"].update(s for s in gem.rxns[rid]["subsystems"] if s)
    edges = []
    for (b, g), e in sorted(acc.items()):
        edges.append(dict(
            met_id=b, met_name=gem.base_name.get(b, ""), ensembl=g, symbol=gem.symbol.get(g, ""),
            n_reactions=len(e["rxns"]),
            role="either" if len(e["roles"]) > 1 else next(iter(e["roles"])),
            is_transport=e["transport"],
            reactions=";".join(sorted(e["rxns"])[:MAX_LISTED_REACTIONS]),
            subsystems=";".join(sorted(e["subsystems"])[:MAX_LISTED_SUBSYSTEMS]),
        ))
    return edges


# ---------------------------------------------------------------- explain one pair

def equation(gem: HumanGEM, rid: str, highlight: str) -> str:
    r = gem.rxns[rid]
    def side(sign):
        terms = []
        for m, c in r["mets"].items():
            if (c < 0) == (sign < 0):
                n = f"{gem.met_name.get(m, m)}[{m[-1]}]"
                if base_id(m) == highlight:
                    n = f"**{n}**"
                terms.append(n if abs(c) == 1 else f"{abs(c):g} {n}")
        return " + ".join(terms) or "(nothing)"
    return f"{side(-1)} {'<=>' if r['reversible'] else '=>'} {side(1)}"


def gene_context(rule: str, gene: str) -> str:
    """'isozyme' if the gene can catalyse on its own (joined by 'or'), 'complex subunit' if it is AND-ed with others."""
    toks = re.findall(r"\(|\)|and|or|ENSG\d+", rule)
    i = toks.index(gene)
    left = toks[i - 1] if i > 0 else None
    right = toks[i + 1] if i + 1 < len(toks) else None
    return "complex subunit (AND)" if "and" in (left, right) else ("sole gene" if len(rule_genes(rule)) == 1 else "isozyme (OR)")


def explain(gem: HumanGEM, gene_q: str, met_q: str) -> int:
    g = gem.resolve_gene(gene_q)
    b = gem.resolve_metabolite(met_q)
    print(f"Human-GEM {gem.version}: {len(gem.rxns):,} reactions, {len(gem.base_name):,} metabolites (compartments collapsed), {len(gem.symbol):,} genes\n")

    print("Step 1  resolve the gene")
    if not g:
        print(f"  '{gene_q}' is not a gene symbol or Ensembl id in Human-GEM -> no edge"); return 1
    g_rxns = [rid for rid, r in gem.rxns.items() if g in rule_genes(r["rule"])]
    print(f"  {gene_q} -> {g} ({gem.symbol.get(g, '?')}); appears in the gene rule of {len(g_rxns)} reactions")
    if not g_rxns:
        print("  the gene is in no reaction rule -> no edge"); return 1

    print("\nStep 2  resolve the metabolite")
    if not b:
        sug = gem.name_suggestions(met_q)
        print(f"  '{met_q}' did not match a name, MAM id, KEGG, HMDB, ChEBI or PubChem id -> no edge")
        if sug:
            print("  names containing it: " + "; ".join(sug))
        return 1
    comps = sorted(gem.met_comp[m] for m in gem.met_name if base_id(m) == b)
    deg = gem.degree.get(b, 0)
    print(f"  {met_q} -> {b} '{gem.base_name.get(b)}'; exists in compartments {','.join(comps)}; appears in {deg} reactions")
    currency = deg > CURRENCY_DEGREE
    print(f"  currency test: {deg} {'>' if currency else '<='} {CURRENCY_DEGREE} -> " + ("CURRENCY, gets no gene edges" if currency else "kept"))

    print(f"\nStep 3  reactions that contain {b} (any compartment) AND have {g} in their gene rule")
    shared = [rid for rid in g_rxns if any(base_id(m) == b for m in gem.rxns[rid]["mets"])]
    if not shared:
        print("  none -> no edge"); return 1
    per_rxn = []
    for rid in shared:
        r = gem.rxns[rid]
        coefs = {m: c for m, c in r["mets"].items() if base_id(m) == b}
        contrib = [(role, tr) for mb, gg, role, tr in reaction_contributions(gem, rid) if mb == b and gg == g]
        print(f"  {rid}  {r['name']}  [{'; '.join(r['subsystems'])}]")
        print(f"      {equation(gem, rid, b)}")
        print(f"      coefficient(s): " + ", ".join(f"{m} {c:+g}" for m, c in coefs.items())
              + f"; {'reversible' if r['reversible'] else 'irreversible'}; gene is {gene_context(r['rule'], g)} in a rule of {len(rule_genes(r['rule']))} genes")
        if contrib:
            roles = {c[0] for c in contrib}
            role = "either" if len(roles) > 1 else roles.pop()
            tr = any(c[1] for c in contrib)
            per_rxn.append((rid, role, tr))
            print(f"      -> contributes role={role}, transport={tr}")
        else:
            print("      -> contributes nothing (metabolite is currency)")

    print("\nStep 4  aggregate over the contributing reactions -> edge")
    if not per_rxn:
        print(f"  no reaction contributes (currency metabolite) -> no edge"); return 1
    roles = {x[1] for x in per_rxn}
    subs = sorted({s for rid, _, _ in per_rxn for s in gem.rxns[rid]["subsystems"] if s})
    edge = dict(met_id=b, ensembl=g, symbol=gem.symbol.get(g, ""), n_reactions=len(per_rxn),
                role="either" if len(roles) > 1 else roles.pop(), is_transport=any(x[2] for x in per_rxn),
                reactions=";".join(sorted(x[0] for x in per_rxn)[:MAX_LISTED_REACTIONS]),
                subsystems=";".join(subs[:MAX_LISTED_SUBSYSTEMS]))
    print(f"  EDGE  {gem.symbol.get(g, g)} -- {gem.base_name.get(b)}")
    for k, v in edge.items():
        print(f"    {k:<13}{v}")
    return 0


# ---------------------------------------------------------------- CLI

def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gem-dir", type=Path, default=DEFAULT_DIR, help=f"where the Human-GEM files live (default {DEFAULT_DIR})")
    sub = p.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("download", help="fetch Human-GEM model files from GitHub")
    d.add_argument("--version", default=GEM_VERSION, help="git tag of SysBioChalmers/Human-GEM (default %(default)s)")
    d.add_argument("--force", action="store_true", help="re-download files that already exist")
    e = sub.add_parser("explain", help="trace how the edge between one gene and one metabolite is computed")
    e.add_argument("gene", help="gene symbol (LDHA) or Ensembl id (ENSG00000134333)")
    e.add_argument("metabolite", help="name (L-lactate), MAM id, KEGG (C00186), HMDB (HMDB0000190), ChEBI or PubChem id")
    b = sub.add_parser("build", help="compute all gene-metabolite edges and write them as CSV")
    b.add_argument("--out", type=Path, help="output CSV (default <gem-dir>/gene_metabolite_edges.csv)")
    a = p.parse_args(argv)

    if a.cmd == "download":
        download(a.gem_dir, a.version, a.force); return 0
    ensure_downloaded(a.gem_dir)
    gem = HumanGEM(a.gem_dir)
    if a.cmd == "explain":
        return explain(gem, a.gene, a.metabolite)
    edges = build_edges(gem)
    out = a.out or a.gem_dir / "gene_metabolite_edges.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(edges[0])); w.writeheader(); w.writerows(edges)
    n_met = len({e["met_id"] for e in edges}); n_gene = len({e["ensembl"] for e in edges})
    print(f"Human-GEM {gem.version}: {len(edges):,} gene-metabolite edges ({n_met:,} metabolites x {n_gene:,} genes) -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
