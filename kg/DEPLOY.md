# Exercise KG: Docker deployment

Materials, methods, outcome and caveats: `kg/METHODS.md`.

Adds two services to the hackathon's Neo4j setup without changing the team's `docker-compose.yaml`:

| Service | URL | What it is |
|---|---|---|
| `neo4j` (existing) | http://localhost:7474 | Loaded with the exercise KG: 77.1k nodes, 1.50M relationships, including 137 diseases, 1,344 drugs and 53 `Reference` nodes |
| `explorer` (new) | http://localhost:8080 | Interactive explorer: Pathway / Gene / Metabolite / Phenotype / Disease tabs, one-click themes, database links on every node, readable node and edge labels, hypothesis cards with a citation on every sentence |

Everything for the KG lives in `kg/`:

| File | Role |
|---|---|
| `kg/deploy.sh` | one command: import, start, create indexes, print counts |
| `kg/docker-compose.kg.yaml` | add-on compose file: mounts `./kg` into `neo4j`, adds `kg-import` (one-off bulk loader) and `explorer` (nginx) |
| `kg/neo4j/docker-import.sh` | runs `neo4j-admin database import full` inside the container |
| `kg/neo4j/schema.cypher` | indexes, applied after the import |
| `kg/explorer/Dockerfile`, `nginx.conf`, `.dockerignore`, `vendor/` | explorer image (static site, no CDN needed) |

## Run it

1. Start Docker Desktop.
2. From the repo root:
   ```bash
   cp .env.example .env     # skip if you already have one; set NEO4J_AUTH=neo4j/<password, 8+ characters>
   bash kg/deploy.sh
   ```
   `deploy.sh` stops Neo4j, bulk-imports `kg/exports_sig` into the `neo4j` database (this **overwrites** that database), starts Neo4j and the explorer, creates the indexes, and prints the node and relationship counts (expect 77,123 and 1,500,207).

Then open http://localhost:7474 (Neo4j Browser) and http://localhost:8080 (explorer). Example queries: `kg/neo4j/queries.cypher` (17, including the disease layer).

## Day to day

```bash
DC="docker compose -f docker-compose.yaml -f kg/docker-compose.kg.yaml"
$DC up -d                      # start (graph already imported)
$DC down                       # stop; data stays in ./neo4j/data
bash kg/deploy.sh --no-import  # rebuild the explorer after editing kg/explorer/
bash kg/deploy.sh              # re-import after rebuilding kg/exports_sig
```
`docker compose up -d` on its own still starts only Neo4j, exactly as before.

## Rebuilding the graph

```bash
S=kg/kg_scripts_omics
python3 $S/build_sig_graph.py --rat processed/rat_pkg --human processed/human_pkg --ppi-dir raw/ppi --out kg/exports_sig
python3 $S/add_disease_layer.py --exports kg/exports_sig --hetionet raw/ppi --rat processed/rat_pkg --human processed/human_pkg
python3 $S/export_explorer_data.py kg/exports_sig kg/explorer/data.json
python3 $S/build_adjacency.py kg/exports_sig kg/adjacency
bash kg/deploy.sh
```
The inputs (`raw/`, `processed/`, about 2 GB) are not in git: `processed/` comes from the MoTrPAC R packages [D2, D8] via `kg/kg_scripts_omics/convert.R` and `convert_h.R`; `raw/ppi/` holds the STRING, Hetionet and LIANA files [R1, R2, R4] listed in `kg/references.md`. `export_explorer_data.py` also reads `kg/refs/references.json`.

## Notes

- **Works offline.** The explorer serves a vendored Cytoscape.js 3.30.2 (MIT, `kg/explorer/vendor/`); without internet the page falls back from IBM Plex to system fonts.
- **Claude features** (Generate insights, Expand with Claude, Ask Claude) only work when the page is opened inside Claude. In the Docker deployment those buttons are disabled; the rule-based hypothesis cards and all citations work.
- **Password.** Neo4j reads `NEO4J_AUTH` only when `./neo4j/data` is first created. If you change it later, change it in Neo4j Browser, or delete `./neo4j/data` and rerun `bash kg/deploy.sh`.
- **Port clash.** `EXPLORER_PORT=9000 bash kg/deploy.sh` if 8080 is taken.
- **Team import files.** The team's `neo4j/import/*.cypher` are still mounted at `/import` and are not touched by `deploy.sh`.
- **Citations.** Every edge has `source_refs`; `MATCH (r:Reference {id:'D1'}) RETURN r` gives the full citation. Library: `kg/refs/references.json`, `kg/references.md`.
- **Metabolite rule.** This graph has no metabolite–gene or metabolite–exercise-group edges by design.
