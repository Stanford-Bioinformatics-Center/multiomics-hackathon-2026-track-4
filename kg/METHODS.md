# Exercise KG: materials, methods, outcome and caveats

*Stanford Multi-omics Hackathon 2026, Track 4 (MoTrPAC Knowledge Graph). Documentation date 2026-09-26.*

Bracketed IDs (for example \[D1\], \[R1\], \[M1\], \[L1\], \[P1\]) refer to the reference library in `refs/references.json` and `docs/references.md`. The same IDs appear on every graph edge (`source_refs`), as `Reference` nodes in Neo4j, and on every sentence of the explorer's hypothesis cards. \[P1\] marks anything computed in this project that has not been peer reviewed.

---

## 1\. Summary

**Question.** Which molecular responses to exercise are shared between trained rats and humans after a single bout, how are they wired together (pathways, protein interactions, signalling), and which of them track physiological change? The frame is "exercise as medicine".

**What was built.**

- A **significant-response knowledge graph**: 77,123 nodes and 1,500,207 relationships. It contains only molecules that change significantly with exercise. They are anchored on pathways and linked to organs, exercise groups, phenotypes, protein interactions, diseases, drugs and literature sources.  
- An **interactive web explorer** with Pathway, Gene, Metabolite and Phenotype tabs:  
  - Every node keeps at most 5 edges by default.  
  - Rule-based hypothesis cards on every tab carry citations at the level of each sentence.  
  - Claude-generated insights and a chat that uses graph tools are available when the page is opened inside Claude.  
- A **Neo4j deployment** in Docker: bulk import, indexes and 14 example queries, plus the explorer served by nginx. It is started with one command (`bash deploy.sh`).

**Data used.**

- **Rat:** the MoTrPAC endurance-training study (6-month-old rats, 1/2/4/8 weeks of progressive treadmill training vs sedentary controls, 20 organs) \[D1, D2\].  
- **Human:** the MoTrPAC pre-suspension acute-exercise study in sedentary adults (one endurance or resistance bout vs non-exercise controls, in muscle, adipose and blood) \[D3–D8\].

---

## 2\. Materials

### 2.1 Datasets

| ID | Dataset | Version used | What was taken |
| :---- | :---- | :---- | :---- |
| D1, D2 | MoTrPAC rat endurance training; R package `MotrpacRatTraining6moData` | GitHub commit f831a4f (2025-08-13), original release | Training-regulated features (training q), per-week graphical state calls, per-tissue differential analysis tables (timewise adj p), graphical-cluster pathway enrichment, feature→gene maps, rat→human gene and phosphosite orthologs, phenotypes |
| D3–D8 | MoTrPAC human pre-suspension acute exercise; R package `MotrpacHumanPreSuspensionAnalysis` | GitHub commit b7695e5 (2026-09-25), data release c2.0 | Differential analysis (`*_DA`, control-adjusted contrasts), CAMERA-PR and PTM-SEA enrichment, feature→gene maps, bundled molecular signatures (MSigDB, PhosphoSitePlus, PTMsigDB), clinical chemistry |
| D9 | MoTrPAC consortium design \[D9\] | — | Study-design context only |

The two packages together give 210 rat and 55 human tables. They were converted from `.rda` to `.tsv.gz` with R 4.3.3 (`scripts/convert.R`, `scripts/convert_h.R`).

### 2.2 Knowledge resources

| ID | Resource | Used for |
| :---- | :---- | :---- |
| R1 | STRING v12, human physical subnetwork (`9606.protein.physical.links.v12.0`, `9606.protein.info.v12.0`) | Physical protein–protein interactions, combined score ≥ 700 |
| R2, R3 | Hetionet v1.0 Gene–interacts–Gene (built on HI-II-14, Lit-BM-13 and other curated sets) | Curated physical interactions |
| R4, R5 | LIANA consensus ligand–receptor resource (OmniPath-derived) | Ligand → receptor edges |
| R6 | PhosphoSitePlus kinase–substrate sets (bundled in \[D8\]) | Kinase → phosphosite edges |
| R7 | PTMsigDB and PTM-SEA | Site-level signatures; human PTM enrichment |
| R8–R14 | MSigDB-curated collections (bundled in \[D8\]): Reactome, WikiPathways, KEGG / KEGG MEDICUS, BioCarta, PID, MitoCarta 3.0, GO | Pathway nodes and membership |
| R15 | RefMet | Metabolite names shared across species |
| R16 | UBERON | Matching tissues across species |
| R17, R18 | Rhea, HMDB | Suggested look-ups in hypothesis tests only (not loaded) |
| R2, R19, R20 | Hetionet v1.0 disease and compound nodes and edges (`hetionet-v1.0-nodes.tsv` and the edge file); Disease Ontology; DrugBank | Disease & drug layer (section 3.13) |

### 2.3 Software

- **Data processing:** Python 3.11 (pandas 3.0, NumPy 2.4, SciPy 1.17, NetworkX 3.6) and R 4.3.3.  
- **Explorer:**  
  - Cytoscape.js 3.30.2 (vendored, MIT licence).  
  - The Claude artifact runtime (`sample` capability) powers the Claude features.  
- **Database:** Neo4j (`neo4j:latest` Docker image, `neo4j-admin database import full`).  
- **Serving:** nginx 1.27-alpine, Docker Compose.  
- **Scripts:**

| Script | Role |
| :---- | :---- |
| `scripts/build_sig_graph.py` | Builds the graph (all methods below) |
| `scripts/export_explorer_data.py` | Writes the compact explorer data (`explorer/data.json`) |
| `scripts/load_graph.py` | NetworkX loader with an example query |
| `kg/neo4j/validate_import.py` | Checks headers and IDs before import |
| `deploy.sh` | Docker deployment |

---

## 3\. Methods

### 3.1 What counts as a significant response

Only significant changes enter the graph. The significance rules come from the source studies, not from this project.

- **Rat** \[D1, D2\]:  
  - A molecule must be training-regulated: MoTrPAC's training q \< 0.05, one test per feature across all sex × week groups.  
  - Its direction at each sex and week comes from MoTrPAC's graphical state call (F/M \= −1, 0, \+1). A regulation edge is created only where the call is non-zero.  
  - The 3% of training-regulated features with no state call, mostly ovary and testes, fall back to a timewise p \< 0.05 at that week \[P1\].  
  - Each edge records which rule applied in `significance_basis`:

| Rule | Edges |
| :---- | :---- |
| State call | 125,126 |
| Timewise p fallback | 2,414 |

- **Human** \[D3, D8\]:  
  - Only control-adjusted contrasts (`exercise_with_controls`) are used: the change in the exercise group minus the change in the control group, each against pre-exercise.  
  - They come from the MoTrPAC dream mixed models \[M2\], with BH adjusted p \< 0.05 within each contrast \[M1\]. This gives 38,405 edges.  
  - Contrasts that are not control-adjusted (exercise without controls, baseline, controls only) are excluded.

### 3.2 Scores and edge weights

- **Headline score \= adjusted p-value.** Fold changes are not used as scores (team decision). Direction is carried by the edge type (`UPREGULATED_IN` / `DOWNREGULATED_IN`).  
  - **Human:** BH adj p within the contrast \[M1\].  
  - **Rat:** MoTrPAC training q \[D1\]. The per-week BH adj p is kept as `timewise_adj_p` \[D2\]. It exceeds 0.05 for 95% of significant rat changes (median 1.0) because each group has about 5 animals \[P1\], which is why it is not the headline score.  
  - z, p and −log10(adj p) are also stored.  
- **Weight (0–1)** \= the percentile of −log10(adj p) within species × assay × organ, with ties broken by |z| \[P1\].  
  - Adjusted p is not comparable across assays: −log10(adj p) reaches 98 for human metabolites but only 7 for Olink proteins. The percentile puts every assay and species on one scale.  
  - The median weight is 0.50 in both species.  
- **Other layers.** Weights are defined in `docs/edge_weights.md`:  
  - enrichment: min(−log10 adj p, 10)/10  
  - STRING: score/1000; Hetionet-only interactions: 0.8  
  - orthologs: 1/number of orthologs  
  - pathway membership via an ortholog: 0.8  
  - co-regulation: |r|  
- **Per-node ranks.** `rank_from_src` and `rank_at_dst` rank each edge by weight within its edge type at each end. They support the 5-edges-per-node view in Neo4j: `WHERE r.rank_from_src <= 5 AND r.rank_at_dst <= 5`.

### 3.3 Nodes and identifiers

- **Molecular features** are species-specific:  
  - Transcript, Protein, PTM site (phospho, acetyl, ubiquityl), Affinity protein (Olink and Luminex), Chromatin region (ATAC) and Methylation region.  
  - Their IDs are `<species>:<assay>:<source feature ID>`.  
  - They map to genes with MoTrPAC's feature-to-gene tables \[D2, D8\]. ATAC and methylation regions carry their relationship to the gene (`distance_to_gene`).  
- **Metabolites** are single nodes keyed by RefMet name (`REFMET:<name>`), shared by rat and human \[R15\]. Metabolites without a RefMet name keep their platform ID.  
- **Genes:** 16,978 human and 11,454 rat. Rat and human genes are linked by `ORTHOLOG_OF` from MoTrPAC's ortholog table \[D2\]. Rat phosphosites map to human sites through `ORTHOLOGOUS_SITE` \[D2\].  
- **Organs, groups and timepoints:**  
  - 20 rat organs and 3 human tissues (muscle, adipose, blood), with `EQUIVALENT_TISSUE` links \[R16\].  
  - 3 exercise groups (rat endurance training; human acute endurance; human acute resistance).  
  - Timepoints are stored with their order and timescale:  
    - rat: weeks of training  
    - human: time relative to the bout (during 20/40 min; post 10 min, 15–45 min, 3.5–4 h, 24 h)

### 3.4 Regulation edges

- Every significant response becomes an `UPREGULATED_IN` or `DOWNREGULATED_IN` edge to the organ and, for non-metabolites, a second edge to the exercise group.  
- Each edge carries:  
  - species, tissue, group, sex (rat F/M; human "both"), timepoint, timepoint order, timescale, assay and omics layer  
  - adj p, z, p, `significance_basis`, weight and `source_refs`

**Metabolite rules**, applied as a design choice and checked automatically at the end of every build:

- no metabolite → exercise group edges  
- no metabolite → gene edges  
- no metabolite → pathway edges  
- metabolites connect only to organs, other metabolites and proteins

### 3.5 Pathways

- **Human curated sets** \[R8–R13, D8\]:  
  - All KEGG MEDICUS, Reactome, WikiPathways, BioCarta, PID and MitoCarta sets are included.  
  - GO terms \[R14\] are included only if they are enriched in at least one human contrast.  
  - Rat genes join these sets through their human ortholog (weight 0.8).  
- **Phosphosite signatures** (PhosphoSitePlus, PTMsigDB) \[R6, R7\] are matched to sites by their 15-amino-acid flanking sequence. Rat sites are matched through the orthologous human site.  
- **Rat KEGG and Reactome terms** come from MoTrPAC's g:Profiler enrichment of graphical clusters \[D1, M4\].  
- **Enrichment edges** (`ENRICHED_UP_IN` / `ENRICHED_DOWN_IN`, adj p \< 0.05) come from:  
  - human CAMERA-PR \[M3, D8\] and PTM-SEA \[R7, D8\]  
  - rat cluster enrichment, with a per-sex direction taken from the cluster state \[D1\]

### 3.6 Interaction layers

All interaction evidence is human. Rat genes reach it through orthologs.

- **`INTERACTS_WITH`** (165,639 edges) merges the two sources for each pair and records `source` and `string_score`:

| Source | Edges |
| :---- | :---- |
| Hetionet only \[R2, R3\] | 101,294 |
| STRING only \[R1\] (physical, combined score ≥ 700\) | 47,162 |
| Both | 17,183 |

STRING protein IDs were mapped to gene symbols with STRING's own protein-info file.

- **`LIGAND_OF`** (3,077): from the LIANA consensus resource \[R4, R5\]. Complexes are split into subunits.  
- **`PHOSPHORYLATES`** (745): a kinase linked to a significant phosphosite, from PhosphoSitePlus \[R6\] matched by flanking sequence. Kinase genes get a node even if the kinase itself does not change, so kinase activity can be inferred from substrates \[L6, L12\].

### 3.7 Metabolite co-regulation

- **Method** \[P1\]: Pearson correlation of response profiles within one species and one tissue.  
  - Rat profile: logFC across 2 sexes × 4 weeks.  
  - Human profile: logFC across group × timepoint (6 points in muscle and adipose, 10 in blood).  
- **Edges are kept** for metabolite–metabolite and metabolite–protein pairs with |r| ≥ 0.9, top 10 partners per metabolite. Weight \= |r|, and the sign is stored.  
- **Count:** 31,676 edges (rat 24,588; human 7,088).  
- **No significance test is applied.** These edges are hypotheses (see Caveats).

### 3.8 Phenotypes

- **Rat** \[D1, D2\]:  
  - Nine measures: VO2max change, body-fat change, lean-mass change, blood-lactate change, body weight, and lateral gastrocnemius, medial gastrocnemius, plantaris and soleus mass.  
  - Each trained group is compared with 8-week sedentary controls within sex, using Welch's t-test \[M5\] with at least 3 animals per group.  
  - The BH correction \[M1\] runs across all rat phenotype tests, and edges are kept at adj p \< 0.05 \[P1\].  
- **Human** \[D3, D8\]:  
  - Nine clinical chemistry markers: lactate, glycerol, ketones, NEFA, cortisol, glucose (metabolite panel); glucagon, creatine kinase, insulin (protein panel).  
  - Taken from the control-adjusted contrasts at adj p \< 0.05, linked to blood.  
  - They are modelled as phenotypes, not metabolites.  
- **Result:** 65 phenotype edges (49 increased, 16 decreased).

### 3.9 Citations and provenance

- **Edges:** every edge has `source_refs`, assigned by rule from the edge type, species, tissue and resource (`cite_edges`).  
- **Reference library:** 53 references with verified DOIs, exported as `Reference` nodes. They are grouped as:

| Group | IDs |
| :---- | :---- |
| Datasets | D1–D9 |
| Resources | R1–R20 |
| Statistical methods | M1–M6 |
| Literature | L1–L17 |
| This project | P1 |

- **Explorer:** each hypothesis card cites every sentence, evidence item and suggested test separately:  
  - Data statements cite the datasets.  
  - Interpretations cite the literature behind them.  
  - Heuristics of this project cite \[P1\] and say so.  
- **Claude:** Claude is instructed to cite only library IDs and to mark anything else `[unsourced]`. Citations outside the library are flagged "unverified" in the page.

### 3.10 Five edges per node

- **In the explorer**, each node keeps at most 5 edges by default; a slider sets 1 to no limit.  
- **Edge priority:**  
  1. The core chain is kept first, pathway → gene → molecule → organ or group, plus enrichment, co-regulation and phenotypes.  
  2. Protein interactions, ligand–receptor, kinase and ortholog edges then fill the remaining slots, strongest first \[P1\].

### 3.11 Hypothesis generation

**Rule-based cards** are computed from what is on screen under the current filters \[P1\]. They are labelled "From the graph" (observations) or "New hypothesis":

| Tab | Cards |
| :---- | :---- |
| Pathway | Genes conserved in direction across species \[L1–L3\]; acute vs trained divergence \[L1, L2\]; interaction hubs \[L7\]; kinase activity from substrates \[L6, L12\]; enrichment in both species |
| Gene | Rat vs human comparison; cross-omics agreement or discordance \[L10, L5\]; responding interaction module \[L7\]; ligand–receptor route \[L4\]; upstream kinase \[L6, L12\] |
| Metabolite | Organs affected; species comparison \[L1, L11\]; systemic signal \[L4, L8\]; tracking proteins \[P1\] |
| Phenotype | Rat physiology; muscle pathways at the same weeks \[L13, L1\]; human clinical chemistry; co-moving blood metabolites \[L8, L9, L11\] |
| Disease | Organs where the disease's genes are over-represented among responders \[M6, M1\]; exercise opposes the disease signature \[L16\]; exercise transiently mimics it \[L17, L11\]; drugs for the disease act on genes exercise moves \[R20\] |

Each hypothesis card proposes a test, and the test is cited \[L14, L15, R17, R18, etc.\].

**Claude** (inside Claude only) works from the current view's evidence:

- **Generate insights** separates existing insights from new hypotheses and gives a novelty call and confidence for each.  
- **Expand** elaborates a single card.  
- **Ask Claude** answers questions with 9 graph tools (including `disease_summary`). Every tool result returns its sources with it.

### 3.12 Explorer themes and labels

- **Themes** load a disease, gene or pathway with the organs that matter for it:

| Theme | Loads | Organs |
| :---- | :---- | :---- |
| Diabetes | Disease: type 2 diabetes mellitus | rat gastrocnemius, vastus lateralis, subcutaneous fat, liver, adrenal; human muscle, adipose, blood |
| Cardiovascular | Disease: coronary artery disease | rat heart, vena cava, blood, adrenal, brown fat; human muscle, adipose, blood |
| Appetite & energy balance | Gene: rat Lep (plasma leptin, with human LEP) | all |
| Brain | Disease: schizophrenia | rat cortex, hippocampus, hypothalamus |
| Reproductive | The estrogen/steroid-hormone pathway with the most responding genes in ovary and testes | rat ovary, testes; human blood |
| Cancer | Disease: colon cancer | rat colon, small intestine; human blood, muscle |

- **Labels** for pathways, organs, phenotypes, diseases and drugs grow when the view is zoomed out, so they stay readable at any zoom.  
- **Focus nodes** (the selected pathway, gene with its orthologs, metabolite or disease) may keep up to 3 × the per-node edge limit; disease → organ findings are kept first.

### 3.13 Disease and drug layer

Built by `scripts/add_disease_layer.py` from Hetionet v1.0 \[R2\], which integrates the Disease Ontology \[R19\] and DrugBank \[R20\]. Rat genes enter through their human ortholog \[D2\].

- **Nodes:** 137 diseases and 1,344 compounds. A compound is kept if it binds a gene that responds to exercise, or treats or palliates a disease.  
- **Curated edges** (limited to genes that respond to exercise in either species):

| Edge | n | Hetionet source |
| :---- | :---- | :---- |
| `ASSOCIATED_WITH` Gene → Disease | 10,677 | DaG: GWAS Catalog, DISEASES, DisGeNET, DOAF |
| `UP_IN_DISEASE` / `DOWN_IN_DISEASE` Gene → Disease | 6,585 / 6,477 | DuG / DdG: disease-vs-control expression (STARGEO) |
| `BINDS` Compound → Gene | 7,429 | CbG: DrugBank, ChEMBL, BindingDB, DrugCentral |
| `TREATS` / `PALLIATES` Compound → Disease | 755 / 390 | CtD / CpD: curated indications |

- **Over-representation** (`DISEASE_GENES_ENRICHED_IN`, Disease → Tissue) \[P1\]:  
  - For each disease and each species × organ, a hypergeometric test \[M6\] asks whether its associated genes are over-represented among the genes that respond to exercise there.  
  - The universe is every gene measured in that organ, taken from all tested features in the MoTrPAC differential-analysis tables \[D2, D8\]. ATAC, methylation and metabolites are left out because the rat package has no full tested-feature list for them.  
  - Tests need at least 5 disease genes measured and 3 responding. BH correction runs across all 1,557 tests \[M1\]; 410 edges pass adj p \< 0.05 with fold \> 1 (373 rat, 37 human).  
- **Direction** (`EXERCISE_OPPOSES` / `EXERCISE_MIMICS`, Disease → Tissue) \[P1\]:  
  - For genes that are dysregulated in the disease (DuG/DdG) and also respond to exercise in an organ, an exact binomial test asks whether exercise moves them the opposite way (or the same way) more often than 50:50.  
  - Tests need at least 10 genes; BH correction runs across 585 tests \[M1\]. 23 organ × disease pairs oppose (15 rat, 8 human) and 29 mimic (11 rat, 18 human).

### 3.14 Validation checks

- The metabolite-rule check passes at the end of every build.  
- Dangling and duplicate edges are removed and counted in `manifest.json`.  
- `validate_import.py` checks the Neo4j headers and ID integrity.  
- An automated browser test covered 51 cards across the Pathway, Gene, Metabolite and Phenotype tabs, and every theme preset on the Disease tab: 0 uncited sentences, evidence items or tests, and 0 page errors.  
- The Docker explorer was tested with no network access.

---

## 4\. Outcome

### 4.1 Graph size

| Nodes | n | Relationships | n |
| :---- | :---- | :---- | :---- |
| TranscriptFeature | 21,901 | IN\_PATHWAY | 845,278 |
| PTMSite | 9,051 | UPREGULATED\_IN | 176,138 |
| ProteinFeature | 2,946 | DOWNREGULATED\_IN | 141,692 |
| Metabolite | 1,833 | INTERACTS\_WITH | 165,639 |
| ChromatinRegion | 1,330 | ENRICHED\_UP\_IN / DOWN\_IN | 39,680 / 11,430 |
| MethylationRegion | 1,052 | MAPS\_TO\_GENE | 37,895 |
| AffinityProtein | 269 | CO\_REGULATED\_WITH | 31,676 |
| Gene | 28,432 | ORTHOLOG\_OF | 11,607 |
| Pathway | 8,731 | LIGAND\_OF | 3,077 |
| Phenotype | 18 | SITE\_ON | 2,323 |
| Tissue | 23 | PHOSPHORYLATES | 745 |
| ExerciseGroup | 3 | ORTHOLOGOUS\_SITE | 170 |
| Reference | 53 | INCREASED\_IN / DECREASED\_IN | 49 / 16 |
| Disease | 137 | MEASURED\_IN / EQUIVALENT\_TISSUE | 12 / 5 |
| Compound | 1,344 | ASSOCIATED\_WITH | 10,677 |
|  |  | UP\_IN\_DISEASE / DOWN\_IN\_DISEASE | 6,585 / 6,477 |
|  |  | BINDS | 7,429 |
|  |  | TREATS / PALLIATES | 755 / 390 |
|  |  | DISEASE\_GENES\_ENRICHED\_IN | 410 |
|  |  | EXERCISE\_OPPOSES / EXERCISE\_MIMICS | 23 / 29 |
| **Total** | **77,123** | **Total** | **1,500,207** |

### 4.2 Significant responses (organ edges)

| Species | Layer | Features | Response edges | Organs |
| :---- | :---- | :---- | :---- | :---- |
| Rat | Transcript | 9,521 | 66,161 | 19 |
| Rat | Protein | 2,609 | 15,764 | 7 |
| Rat | Phosphosite | 2,210 | 11,766 | 7 |
| Rat | Acetylsite | 2,188 | 9,793 | 2 |
| Rat | Ubiquitylsite | 134 | 461 | 2 |
| Rat | Metabolite | 1,471 | 12,528 | 19 |
| Rat | Immunoassay | 46 | 387 | 16 |
| Rat | ATAC region | 1,330 | 6,115 | 8 |
| Rat | Methylation region | 1,052 | 4,565 | 8 |
| Human | Transcript | 12,380 | 29,719 | 3 |
| Human | Phosphosite | 4,519 | 6,470 | 2 |
| Human | Protein | 337 | 344 | 2 |
| Human | Olink protein | 223 | 340 | 1 |
| Human | Metabolite | 529 | 1,532 | 3 |

- **Totals:** 20,561 rat and 17,988 human features respond significantly (127,540 rat and 38,405 human organ edges).  
- **Pathways:** 8,731 pathway nodes:

| Source | Pathway nodes |
| :---- | :---- |
| GO, enriched only | 3,381 |
| Reactome | 1,691 |
| WikiPathways | 791 |
| PTMsigDB | 622 |
| KEGG MEDICUS | 619 |
| Rat Reactome | 566 |
| BioCarta | 292 |
| Rat KEGG | 286 |
| PID | 196 |
| MitoCarta | 148 |
| PhosphoSitePlus | 139 |

### 4.3 Deliverables

| Deliverable | Where |
| :---- | :---- |
| Neo4j bulk-import files | `kg/exports_sig/` |
| Schema, weights, references | `kg/schema/SCHEMA_v0.2_significant.md`, `docs/edge_weights.md`, `docs/references.md` |
| Neo4j indexes and 17 example queries | `kg/neo4j/schema.cypher`, `kg/neo4j/queries.cypher` |
| Web explorer | `explorer/` (Docker: [http://localhost:8080](http://localhost:8080); also published as a private Claude artifact) |
| Docker deployment | `docker-compose.kg.yaml`, `deploy.sh`, `DEPLOY.md` |
| Adjacency matrices (gene/metabolite × contrast, gene × tissue, metabolite × protein), team format | `kg/adjacency/` (see `kg/adjacency/ADJACENCY.md`) |
| Build pipeline | `scripts/` |

---

## 5\. Caveats and limitations

**Study design**

- **Rat and human differ in exposure.** Rat is chronic training; human is one acute bout. Opposite directions can be real biology, since acute signals are transient and add up over training \[L1, L2, L3\]. The divergence cards are hypotheses, not contradictions.  
- **Different organs and timescales.** Rat covers 20 organs; human covers 3 tissues. Weeks of training and hours after a bout are not comparable time axes.  
- **Human sex.** Human edges are "both sexes" (group-level contrasts); rat edges are sex-specific.

**Disease layer**

- **Disease-named is not disease-proven.** Over-representation shows that exercise touches many genes linked to a disease; it does not show that exercise changes the disease. Hetionet's disease genes are dominated by immune and signalling hubs, so immune-rich organs (rat lung, brown fat, spleen) and blood light up for many unrelated diseases \[R2, P1\].  
- **Disease signatures come from other tissues.** Hetionet's up/down disease genes (STARGEO) were measured in disease-relevant tissues of patients, not in the organ being compared. "Opposes" and "mimics" compare directions only \[R2, P1\].  
- **One bout can look like disease in blood.** Most human "mimics" are in blood (malaria, Crohn's disease, COPD), consistent with the transient immune and inflammatory response to acute exercise \[L17, L11\]. It is not evidence of harm.  
- **Rat enters through human orthologs**, and the rat universe excludes ATAC and methylation regions.

**Statistics**

- **Two significance systems.** Rat uses MoTrPAC's training q plus state calls. Human uses BH per contrast. A rat week call is not a per-week test \[D1\]. The 3% of rat features without a state call use a nominal timewise p \< 0.05, which is weaker; filter on `significance_basis` to exclude them.  
- **Weights are ranks within an assay.** A weight of 0.9 means "strong for this assay and organ", not "stronger than 0.8 in another assay". Adjusted p-values are not comparable across assays or species \[P1\].  
- **Co-regulation edges are not tested.** They rest on 6–10 points per profile with no multiple-testing control. Treat them as hypotheses \[P1\].  
- **Phenotype edges are group-level.** Rat phenotypes compare trained groups with controls; they are not per-animal correlations with molecules. Linking phenotypes to molecules is therefore by timing, not by individual.

**Mapping and annotation**

- **All interactions are human.** Protein–protein, ligand–receptor and kinase edges between rat genes are inferred through orthologs and flagged in the explorer. STRING physical links include predicted and transferred evidence \[R1\]; Hetionet has literature bias \[R2, R3\].  
- **Kinase activity is inferred** from the direction of annotated substrate sites. PhosphoSitePlus coverage is uneven across kinases \[R6, L12\]. Rat sites depend on flanking-sequence orthology (170 orthologous sites map).  
- **Rat and human pathway IDs are not unified.** Rat KEGG and Reactome terms and human MSigDB sets are separate nodes. Rat genes still reach human sets through orthologs.  
- **Metabolite identity.** Only RefMet-named metabolites are shared between species. Unnamed features stay species-specific. Human clinical chemistry is modelled as phenotype, so, for example, plasma lactate appears both as a clinical phenotype and, where measured untargeted, as a metabolite.  
- **Genome builds.** Rat features follow the original package release (believed to be Rn6); the Data Hub's newer Rn7 bundles were not used. ATAC and methylation regions are only as current as that release.

**Scope and tools**

- **Epigenomics.** Rat ATAC and methylation include only the regions MoTrPAC reports as changed at 5% FDR. Human ATAC-seq and MethylCap-seq differential results are **not included** (see section 6).  
- **Claude-generated content** is labelled as such, and anything outside the reference library is flagged. The novelty calls come from the model's own knowledge and need checking against the literature before use.  
- **Docker deployment.** It was validated for configuration and serving, but the full image build and Neo4j import had not been run end to end at the time of writing. Run `bash deploy.sh` once to confirm.

---

## 6\. Datasets available but not included

| Dataset | Why not included |
| :---- | :---- |
| Human epigenomics differential results: ATAC (blood cells, muscle), MethylCap (blood, muscle, adipose) \[D8\] | Not in the R package (only group means are); needs the Data Hub human Epigenomics bundle (\~1.4 GB), which could not be downloaded here |
| Full rat epigenomics (all regions) | Data Hub rat Epigenomics bundle (14.6 GB); only significant regions are needed by design |
| Human splicing differential results \[D8\] | Available; not yet added |
| Human cluster enrichment by CAMERA (`FCM_CAMERA`) \[D8\] | Available; only the ORA version is used |
| Human group means and SDs (`*_SUM_STATS`) \[D8\] | Available; effect sizes not used by team choice |
| Rat per-animal data (normalised data, raw counts, and \~500 unused phenotype columns) \[D2\] | Available; would allow per-animal phenotype–molecule correlation |
| Phosphosite enrichment input (`PTMSEA_INPUT`) \[D8\] | Skipped at conversion (needs cmapR); the enrichment results themselves are included |
| Non-significant features; human contrasts without controls | Excluded by design |
| Human individual-level and full clinical data | Controlled access (dbGaP) |
| STRING rat (10116) and STRING functional links; Hetionet disease and drug edges; MetaMEx \[L3\] | Not added; possible extensions (a disease/drug layer would support the "exercise as medicine" framing) |

---

## 7\. Reproducing

\# 1\. inputs (not in git; \~2 GB)

git clone https://github.com/MoTrPAC/MotrpacRatTraining6moData raw/rat/MotrpacRatTraining6moData            \# commit f831a4f

git clone https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis raw/human/MotrpacHumanPreSuspensionAnalysis   \# commit b7695e5

Rscript scripts/convert.R && Rscript scripts/convert\_h.R                                           \# \-\> processed/rat\_pkg, processed/human\_pkg

\# raw/ppi/: 9606.protein.physical.links.v12.0.txt.gz, 9606.protein.info.v12.0.txt.gz (STRING), hetionet\_edges.sif.gz, hetionet\_nodes.tsv

\#          (github.com/hetio/hetionet, hetnet/tsv/hetionet-v1.0-nodes.tsv), liana\_omni\_resource.csv

\# 2\. build

python3 scripts/build\_sig\_graph.py \--rat processed/rat\_pkg \--human processed/human\_pkg \--ppi-dir raw/ppi \\

        \--out kg/exports\_sig \--r-min 0.9 \--top-k 10 \--go enriched \--string-min 700

python3 scripts/add\_disease\_layer.py \--exports kg/exports\_sig \--hetionet raw/ppi \--rat processed/rat\_pkg \--human processed/human\_pkg   \# needs raw/ppi/hetionet\_nodes.tsv

python3 scripts/export\_explorer\_data.py kg/exports\_sig explorer/data.json

python3 scripts/build\_adjacency.py kg/exports\_sig kg/adjacency

\# 3\. deploy (Docker Desktop running, .env with NEO4J\_AUTH=neo4j/\<password\>)

bash deploy.sh

`convert.R` and `convert_h.R` read from `raw/…/data` and write to `processed/…` by default; pass other paths as arguments (`Rscript scripts/convert.R <src> <out>`).

---

## 8\. References

The full citations with DOIs are in `docs/references.md`. The IDs used here are:

- **Data:** D1 MoTrPAC Nature 2024; D2 rat R package; D3–D7 MoTrPAC human preprints (Katz, Keshishian, Ahn, Robbins, Brandt; bioRxiv 2026); D8 human R package; D9 Sanford et al., Cell 2020\.  
- **Resources:** R1 STRING v12; R2 Hetionet; R3 HI-II-14 / Lit-BM-13; R4 LIANA; R5 OmniPath; R6 PhosphoSitePlus; R7 PTMsigDB / PTM-SEA; R8 MSigDB; R9 Reactome; R10 WikiPathways; R11 KEGG; R12 PID; R13 MitoCarta 3.0; R14 GO; R15 RefMet; R16 UBERON; R17 Rhea; R18 HMDB; R19 Disease Ontology; R20 DrugBank.  
- **Statistical methods:** M1 Benjamini–Hochberg; M2 dream; M3 CAMERA; M4 g:Profiler; M5 Welch; M6 Fisher (exact / hypergeometric test).  
- **Literature:** L1 Egan & Zierath 2013; L2 Perry et al. 2010; L3 Pillon et al. 2020; L4 Chow et al. 2022; L5 Hoffman et al. 2015; L6 Casado et al. 2013; L7 Barabási et al. 2011; L8 Brooks 2018; L9 Horowitz & Klein 2000; L10 Liu, Beyer & Aebersold 2016; L11 Contrepois et al. 2020; L12 Hernandez-Armenta et al. 2017; L13 Holloszy & Coyle 1984; L14 Nikolić et al. 2017; L15 Narkar et al. 2008; L16 Pedersen & Saltin 2015; L17 Nieman & Wentz 2019\.  
- **This project:** P1 (`scripts/build_sig_graph.py`).