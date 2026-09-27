# Scores, significance and edge weights — exercise-kg (v0.5)

Bracketed IDs refer to `docs/references.md`. Every edge also stores its sources in `source_refs`.

Every edge has four kinds of information:

- **`significance_basis`**: the rule that let it into the graph.
- **Statistics**: adjusted p, z, p.
- **`weight`** (0–1): for ranking.
- **`rank_from_src` / `rank_at_dst`**: the edge's rank at each end, used for the 5-edges-per-node view.

## Molecular responses (`UPREGULATED_IN` / `DOWNREGULATED_IN`)

The team reports **adjusted p-values, not fold changes**. Direction is carried by the edge type.

| | Human (one acute bout) | Rat (8-wk training) |
|---|---|---|
| Data | Control-adjusted contrasts, muscle / adipose / blood [D3, D4, D5, D6, D8] | Training-regulated features [D1, D2] |
| Model | Mixed model (dream) [M2] | MoTrPAC timewise models and training test [D1] |
| Gate | BH adj p < 0.05 within each contrast [M1, D8] | training q < 0.05 plus MoTrPAC per-week state call [D1]; 3% of features without a state call use timewise p < 0.05 [D2, P1] |
| `adj_p_value` | BH adj p [M1] | MoTrPAC training q [D1] |
| Also stored | z, p, `neg_log10_adj_p` | z, timewise p, `timewise_adj_p` (per-week BH adj p [D2, M1]; none released for ATAC or methylation [D2]) |
| `weight` | percentile of −log10(adj p) within species × assay × organ [P1] | same [P1] |

**Why rat uses training q.** MoTrPAC tests each rat molecule once across all sex × week groups, then assigns a direction to each week [D1]. Per-week adjusted p-values [D2] exceed 0.05 for 95% of significant rat changes (median 1.0), as computed in this project [P1]. They are kept as `timewise_adj_p`.

**Why the weight is a percentile.** In this graph, −log10(adj p) reaches 98 for human metabolites but only 7 for Olink proteins [D8, P1]. The multiple-testing correction [M1] and sample size both shift adjusted p between assays. A within-assay percentile puts every assay and species on one scale; the median weight is 0.50 in both species [P1].

## Other layers

| Edge | Gate / basis | Weight | Sources |
|---|---|---|---|
| `ENRICHED_UP_IN` / `ENRICHED_DOWN_IN` | adj p < 0.05 | min(−log10 adj p, 10)/10 [P1] | human CAMERA-PR [D8, M3], PTM-SEA [D8, R7]; rat g:Profiler on graphical clusters [D1, M4] |
| `IN_PATHWAY` | curated membership | 1 direct, 0.8 via ortholog [P1] | MSigDB bundle [R8, D8]; Reactome [R9], WikiPathways [R10], KEGG [R11], PID [R12], MitoCarta [R13], GO [R14], PhosphoSitePlus/PTMsigDB [R6, R7]; rat orthologs [D2] |
| `INTERACTS_WITH` | STRING physical combined score ≥ 700 and/or curated | STRING score/1000; Hetionet-only 0.8 [P1] | STRING v12 [R1]; Hetionet [R2] from HI-II-14/Lit-BM-13 and others [R3] |
| `LIGAND_OF` | curated consensus | 1 | LIANA [R4] built on OmniPath [R5] |
| `PHOSPHORYLATES` | curated | 1 direct, 0.8 via orthologous site [P1] | PhosphoSitePlus [R6] via the MoTrPAC bundle [D8]; rat site orthology [D2] |
| `CO_REGULATED_WITH` | Pearson \|r\| ≥ 0.9, top 10 per metabolite, **no significance test** | \|r\| | computed here [P1] from response profiles [D1, D2, D3, D8] |
| `ORTHOLOG_OF`, `ORTHOLOGOUS_SITE` | curated orthology | 1/number of orthologs | [D2] |
| `MAPS_TO_GENE`, `SITE_ON` | annotation | 1; ATAC/methylation by distance [P1] | [D2, D8] |
| `EQUIVALENT_TISSUE` | anatomy | 1 exact, 0.5 proxy [P1] | UBERON [R16] |
| Rat phenotypes (`INCREASED_IN` / `DECREASED_IN`) | Welch t-test vs sedentary controls, BH adj p < 0.05 | percentile of −log10(adj p) [P1] | [D1, D2, M5, M1] |
| Human clinical chemistry | adj p < 0.05 | percentile of −log10(adj p) [P1] | [D3, D8, M1] |
| `ASSOCIATED_WITH`, `UP_IN_DISEASE`, `DOWN_IN_DISEASE` | curated (responding genes only) | 1 | Hetionet [R2], Disease Ontology [R19] |
| `BINDS`, `TREATS`, `PALLIATES` | curated | 1 | Hetionet [R2], DrugBank [R20] |
| `DISEASE_GENES_ENRICHED_IN` | hypergeometric, BH adj p < 0.05 across 1,557 disease × organ tests | percentile of −log10(adj p) [P1] | [R2, R19, M6, M1, P1] + organ dataset refs |
| `EXERCISE_OPPOSES` / `EXERCISE_MIMICS` | exact binomial on direction agreement, BH adj p < 0.05 across 585 tests | fraction opposite (or same) [P1] | [R2, R19, M1, P1] + organ dataset refs |

## Five edges per node

**In the explorer**, each node keeps at most 5 edges by default; a slider sets the limit anywhere from 1 to no limit [P1]. Edges are picked strongest first, with the core chain given priority: pathway → gene → molecule → organ or group, then enrichment, co-regulation and phenotypes. Protein interactions, ligand–receptor links, kinase links and orthologs fill whatever capacity is left [P1].

**In Neo4j**, use `WHERE r.rank_from_src <= 5 AND r.rank_at_dst <= 5`. The ranks are computed per edge type by `weight`.

## Caveats

- **Co-regulation edges are hypotheses** [P1].
- **Rat–human comparisons contrast chronic training with one acute bout.** Acute signals are transient and add up over training, so opposite directions can be real biology [L1, L2, L3].
