# Exercise KG — MoTrPAC multi-omics knowledge graph

**Exercise as medicine:** which molecular responses to exercise are shared between trained rats and humans after one bout, how they connect through pathways, protein interactions and signalling, and how they relate to physiology, diseases and drugs. Every link carries a citation.

|  |  |
| :---- | :---- |
| **Graph** | 77,123 nodes · 1,500,207 relationships (significant exercise responses only) |
| **Species** | Rat: 8-week endurance training, 20 organs \[D1, D2\]. Human: one endurance or resistance bout, muscle / adipose / blood \[D3–D8\] |
| **Layers** | transcripts, proteins, phospho/acetyl/ubiquityl sites, Olink/immunoassay, chromatin, methylation, metabolites; pathways; protein interactions; ligand–receptor; kinase–substrate; phenotypes; **diseases and drugs** |
| **Scores** | adjusted p-values (no fold changes); edge weight \= within-assay percentile of −log10(adj p) |
| **Citations** | 53-reference library; every edge has `source_refs`; every hypothesis sentence is cited |

## Quick start

cp .env.example .env          \# set NEO4J\_AUTH=neo4j/\<password, 8+ characters\>

bash deploy.sh                \# Neo4j on http://localhost:7474, explorer on http://localhost:8080

Details: [`DEPLOY.md`](http://DEPLOY.md).

## What's in the repo

| Path | Contents |
| :---- | :---- |
| `kg/exports_sig/` | Neo4j bulk-import CSVs (15 node files, 26 relationship files) \+ `manifest.json` |
| `kg/adjacency/` | 7 adjacency matrices in the team format (gene/metabolite × contrast, gene × tissue, metabolite × protein) — see `ADJACENCY.md` |
| `kg/neo4j/` | `docker-import.sh`, `schema.cypher` (indexes), `queries.cypher` (17 examples), `validate_import.py` |
| `kg/schema/` | Schema of the graph (`SCHEMA_v0.2_significant.md`) |
| `explorer/` | Web explorer (`index.html`, `data.json`, vendored Cytoscape.js, Dockerfile, nginx config) |
| `docs/` | `METHODS.md` (materials, methods, outcome, caveats), `edge_weights.md`, `references.md` |
| `refs/references.json` | Reference library (D \= data, R \= resources, M \= methods, L \= literature, P \= this project) |
| `scripts/` | Build pipeline (below) and `record_demo.py` |
| `deploy.sh`, `docker-compose.kg.yaml` | Docker deployment (add-on to the team's `docker-compose.yaml`) |

## Pipeline

| Step | Script | Output |
| :---- | :---- | :---- |
| 1 | `scripts/convert.R`, `scripts/convert_h.R` | MoTrPAC rat and human R packages → TSV (`processed/`) |
| 2 | `scripts/build_sig_graph.py` | significant-response graph → `kg/exports_sig/` |
| 3 | `scripts/add_disease_layer.py` | Hetionet diseases, drugs and disease tests → `kg/exports_sig/` |
| 4 | `scripts/export_explorer_data.py` | `explorer/data.json` (incl. database IDs and the disease layer) |
| 5 | `scripts/build_adjacency.py` | `kg/adjacency/` |
| 6 | `bash deploy.sh` | Neo4j \+ explorer in Docker |
| — | `scripts/load_graph.py` | load the graph into NetworkX |
| — | `scripts/record_demo.py` | re-record the captioned demo video |

Inputs (`raw/`, `processed/`, \~2 GB) are not in git; `docs/METHODS.md` §7 gives the commands and sources.

## Explorer

- **Tabs:** Pathway, Gene, Metabolite, Phenotype, **Disease**.  
- **Themes:** Diabetes, Cardiovascular, Appetite & energy balance, Brain, Reproductive, Cancer — each loads a disease, gene or pathway with the relevant organs.  
- **Filters:** species, organ presets, timepoint, sex, direction, omics layer, interaction type (incl. disease & drugs), minimum strength, organs vs exercise groups, max edges per node (5 by default), max molecules.  
- **Details:** click any node or edge for every significant response, adjusted p-values, sources, and links to NCBI Gene, UniProt, Ensembl, RGD, NCBI Protein, UCSC, KEGG, PubChem, MSigDB, Reactome, UBERON, Disease Ontology and DrugBank.  
- **Hypotheses:** rule-based cards on every tab, green \= observed in the graph, amber \= new hypothesis; each sentence, evidence item and test is cited.  
- **Claude** (when opened inside Claude): Generate insights, Expand a card, Ask Claude with 9 graph tools — cited to the same library.

## Latest changes

- **Disease & drug layer** from Hetionet v1.0 \[R2\] (Disease Ontology \[R19\], DrugBank \[R20\]): 137 diseases, 1,344 drugs; gene–disease associations and disease signatures for responding genes; drug–target and drug–disease links.  
- **New analyses** \[P1\]: per species × organ over-representation of disease genes among exercise responders (hypergeometric \[M6\], BH \[M1\]; 410 significant), and whether exercise **opposes** (23) or **mimics** (29) a disease's expression signature (binomial test, BH).  
- **Explorer:** Disease tab, theme presets, disease/drug cards and details, `disease_summary` tool for Claude, database links on every node, labels that stay readable when zoomed out, focus nodes keep more edges.  
- **Neo4j:** Disease/Compound indexes and queries 15–17.  
- **References:** R19, R20, M6, L16, L17 added (53 total, DOIs verified).  
- **Docs:** `METHODS.md` §3.12–3.13 and disease-layer caveats; schema, edge weights, deployment and Neo4j READMEs updated.

## Caveats (short)

Rat training vs one human bout are different exposures; adjusted p-values are not comparable across assays; co-regulation edges are untested hypotheses; all interaction evidence is human (rat via orthologs); disease over-representation shows that exercise touches disease genes, not that it changes the disease; one exercise bout can transiently resemble inflammatory disease signatures in blood. Full list: `docs/METHODS.md` §5.