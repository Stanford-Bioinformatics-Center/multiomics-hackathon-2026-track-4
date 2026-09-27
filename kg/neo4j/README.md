# Neo4j for the exercise KG

**Docker (hackathon repo setup):** run `bash kg/deploy.sh` from the repo root; see `kg/DEPLOY.md`. The steps below are for Neo4j Desktop.

## One-time setup (Mac)

1. Install **Neo4j Desktop** (free): https://neo4j.com/download/
2. Create a project, then **Add → Local DBMS** (Neo4j 5.x). Set a password. **Don't start it yet.**
3. On the DBMS, open the **"…" menu → Terminal**. The terminal opens inside the DBMS folder.
4. Run:

   ```bash
   NEO4J_HOME=$PWD <repo>/kg/neo4j/import.sh
   ```

   This imports `kg/exports_sig` (77.1k nodes, 1.50M relationships) in about a minute.
   To import the full v0.1 graph instead, add `exports` at the end of the command.
5. **Start** the DBMS, click **Open** (Neo4j Browser), and paste the lines from `schema.cypher` to create the indexes.
6. Try the queries in `queries.cypher`. Query 2 draws oxidative phosphorylation in muscle as a graph.

## Tips

- **Styling:** in Neo4j Browser, click a label chip (e.g. `Gene`) under the result to set its colour, size and caption (use `symbol` for genes, `name` for pathways).
- **Point-and-click exploration:** open **Neo4j Bloom** from Desktop. Search "Gene PPARGC1A" and expand its neighbours.
- **Performance:** keep `LIMIT` on visual queries (300–500 nodes draws comfortably). Tables can be as big as you like.
- **Check before importing:** `python validate_import.py ../exports_sig` checks the files against the bulk-import rules without Neo4j.
- **Re-importing:** `import.sh` overwrites the database, so stop the DBMS first.

## Model cheat-sheet

| Label | Key properties |
|---|---|
| `MolecularFeature` + one of `TranscriptFeature`, `ProteinFeature`, `PTMSite`, `AffinityProtein`, `ChromatinRegion`, `MethylationRegion`, `Metabolite` | `id`, `species`, `omics_layer`, `assay`, `source_id` |
| `Gene` | `id` (ENSEMBL:…), `symbol`, `species`, `is_TF` |
| `Pathway` | `id`, `name`, `database`, `source_species` |
| `Tissue`, `ExerciseGroup`, `Phenotype` | `id`, `name`, `species` |

Regulation relationships (`UPREGULATED_IN`, `DOWNREGULATED_IN`, `ENRICHED_UP_IN`, `ENRICHED_DOWN_IN`, `INCREASED_IN`, `DECREASED_IN`) carry `species`, `tissue`, `group`, `sex`, `timepoint`, `timescale`, `assay`, `logFC`, `z`, `p_value`, `adj_p_value` and `weight`.


## Disease & drug layer

`Disease` and `Compound` nodes come from Hetionet [R2]. Queries 15–17 in `queries.cypher` show where disease genes are over-represented among exercise responders, where exercise opposes or mimics a disease signature, and which drugs for a disease target genes that exercise moves. Method: `docs/METHODS.md` §3.13.
