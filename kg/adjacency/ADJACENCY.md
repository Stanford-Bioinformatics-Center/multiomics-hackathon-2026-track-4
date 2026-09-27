# Exercise KG adjacency matrices: contents, rules and provenance

These seven matrices were built on 2026-09-26 from the significant-response graph (`kg/exports_sig/`) by `scripts/build_adjacency.py`. They use the same layout as the other adjacency files in `neo4j/`: rows × columns, the first column holds the row label, and every matrix comes with a `_rows.csv`, a `_cols.csv` and an `_edges.csv` (one row per non-zero cell). `_index.csv` lists all seven. Bracketed IDs refer to `docs/references.md`.

| File | Shape | Cell value | Non-zero cells |
|---|---|---|---|
| `human_gene_contrast_adjacency_12815x21.csv` | 12,815 genes × 21 contrasts | signed −log10(adj p) | 32,735 (12.2%): 19,338 up, 13,397 down |
| `rat_gene_contrast_adjacency_11283x150.csv` | 11,283 genes × 150 contrasts | signed −log10(training q) | 98,476 (5.8%): 52,836 up, 45,640 down |
| `human_metab_contrast_adjacency_529x21.csv` | 529 metabolites × 21 contrasts | signed −log10(adj p) | 1,532 (13.8%) |
| `rat_metab_contrast_adjacency_1471x140.csv` | 1,471 metabolites × 140 contrasts | signed −log10(training q) | 12,446 (6.0%) |
| `gene_tissue_adjacency_15894x23.csv` | 15,894 human-ortholog genes × 23 species:tissues | signed −log10(adj p), strongest response at any timepoint | 39,834 (10.9%) |
| `human_metab_protein_coreg_adjacency_431x426.csv` | 431 metabolites × 426 proteins | Pearson r | 2,875 (1.6%) |
| `rat_metab_protein_coreg_adjacency_1153x1980.csv` | 1,153 metabolites × 1,980 proteins | Pearson r | 12,264 (0.5%) |

**Cell value.** In the response matrices, the value is sign(direction) × −log10(adjusted p). Positive means up, negative means down, and 0 means no significant change. No fold changes are used [P1].

- **Human:** the adjusted p is BH within a control-adjusted contrast [D3, D8, M1, M2].
- **Rat:** it is MoTrPAC's training q [D1]. Training q is one test per molecule across all sexes and weeks, so a rat row has the same magnitude in every week where it changes. The sign comes from MoTrPAC's per-week, per-sex call, which means the rat matrices show **where and when** a molecule changes, not how strongly at each week. For per-week strength, the `_edges.csv` files carry `z`.

Only significant responses enter (section 1 of `docs/METHODS.md`).

**Row order.** Rows are sorted from strongest to weakest (by the maximum |cell|). Columns are ordered by tissue, then contrast, then timepoint.

## 1. Gene × contrast (`human_gene_…`, `rat_gene_…`)

- **Rows:** every gene with at least one significant molecule. A molecule can be a transcript, protein, phospho-, acetyl- or ubiquityl-site, Olink/immunoassay protein, or an ATAC or methylation region.
  - Molecules are mapped to one gene each, the best-weighted gene in MoTrPAC's feature-to-gene table [D2, D8].
  - Row labels are gene symbols. A duplicated symbol gets its Ensembl ID appended.
- **Columns:**
  - **Human:** `tissue|contrast|timepoint`.
    - Contrast `EE-CON` = acute endurance minus control; `RE-CON` = acute resistance minus control.
    - Muscle and adipose have only post-exercise timepoints. Blood also has `during_20_min` and `during_40_min` (endurance only).
  - **Rat:** `tissue|sex|week`, where `F`/`M` × `1w`, `2w`, `4w`, `8w` compare training with sedentary controls.
- **Aggregation:** when several molecules of one gene fall in the same cell, the cell keeps the molecule with the smallest adjusted p. The `_edges.csv` file then records:
  - that molecule (`best_molecule`) and its `layer`
  - `n_molecules`, `n_up` and `n_down`, which flag cells where the gene's molecules disagree in direction
  - `z`, `significance_basis` and `source_refs`
- **Row metadata** (`_rows.csv`): gene ID, symbol, Entrez ID, human ortholog (rat only), number of non-zero cells, the omics layers that respond, and the best adjusted p.

## 2. Metabolite × contrast (`human_metab_…`, `rat_metab_…`)

- **Rows:** metabolites with at least one significant response, labelled by RefMet name [R15]. A metabolite without a RefMet name keeps its platform ID.
  - RefMet metabolites are one node shared by rat and human, so the same name can appear in both matrices. `_rows.csv` marks these with `shared_with_other_species`.
  - When one RefMet compound is measured on several platforms in the same tissue, the cell keeps the smallest adjusted p.
- **Columns:** as for the gene matrices.
- **Clinical chemistry is not included.** Lactate, glucose, NEFA, glycerol, ketones, cortisol, insulin, glucagon and CK from the clinical panels are Phenotype nodes in this graph. Untargeted "Lactic acid" is a metabolite and is included.
- This is the adjusted-p counterpart of the team's `exercise_metab_246x10_weighted_log2FC_fdr05.csv`. It uses signed −log10(adj p) instead of log2FC and covers all three tissues and all significant metabolites.

## 3. Gene × tissue (`gene_tissue_…`)

- **Rows:** human genes. Rat genes are placed on the row of their first human ortholog [D2]. Rat genes without a human ortholog are **left out** of this matrix but stay in the rat gene × contrast matrix.
  - 6,989 genes respond in both species, 5,826 in human only and 3,079 in rat only (`responds_in_human` and `responds_in_rat` in `_rows.csv`).
- **Columns:**
  - human `muscle`, `adipose` and `blood`
  - 20 rat organs [D1], with rat vastus lateralis and gastrocnemius alongside human muscle, rat subcutaneous white fat alongside human adipose, and rat blood and plasma alongside human blood [R16]
- **Cell:** the strongest significant response in that species and tissue at any timepoint (and, for rat, either sex), with its sign.
  - The `_edges.csv` file gives `n_up` and `n_down`. A cell can hide a response that goes in opposite directions at different times, such as rat female down and male up.

## 4. Metabolite × protein co-regulation (`human_metab_protein_…`, `rat_metab_protein_…`)

- **Rows:** metabolites. **Columns:** protein features (global proteome, Olink and immunoassay), labelled `symbol|assay`; a duplicated label gets the protein accession appended.
- **Cell:** the Pearson r of the two response profiles within one tissue [P1].
  - Rat profiles: logFC over 2 sexes × 4 weeks.
  - Human profiles: logFC over contrast × timepoint.
  - Only pairs with |r| ≥ 0.9 are kept, with each metabolite's top 10 partners. For a pair that appears in several tissues, the strongest one is kept, and `_edges.csv` names it.
- **No significance test is applied.** These are hypotheses about producers, transporters or consumers, not established links.
- They are the only metabolite–protein links in this graph. By design there is no metabolite–gene matrix: the graph has no metabolite–gene edges. The Human-GEM-based `metab_gene_adjacency_*` files in `neo4j/` are a separate, curated source.

## Caveats

- **Rat and human are not comparable per cell.** Rat is 8 weeks of training with a feature-level q; human is one bout with a per-contrast BH adjusted p. Compare the pattern (sign, tissue, timing), not the magnitudes across species [D1, D8, L1].
- **Cell magnitudes depend on the assay and on sample size** [M1]. A large value in the metabolomics columns is not directly comparable with one in Olink. For a scale that works across assays, use the within-assay percentile `weight` on the graph edges (`docs/edge_weights.md`) [P1].
- **One gene per molecule.** ATAC and methylation regions are assigned to their nearest annotated gene [D2, D8]. That assignment is positional, not functional.
- **Rat features without a MoTrPAC state call** (3%) enter at timewise p < 0.05 (`significance_basis` in `_edges.csv`) [D1, P1].

## Rebuild

```bash
python3 scripts/build_adjacency.py kg/exports_sig kg/adjacency
```
