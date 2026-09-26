# Metabolite x gene adjacency matrix: subset rules and provenance

Built 2026-09-26 from the MoTrPAC 2026 human acute-exercise release (aggregate tables), Human-GEM and GTEx v8.
Transfer bundle: `results/transfer/`. Working copies: `results/metab/adjacency/tierA_metab_x_gtex_gene__*.csv`.

| File | Shape | Content |
|---|---|---|
| `metab_gene_adjacency_246x647.csv` | 246 rows x 647 cols (+ header, + row-label column) | binary matrix, `1` = gene product acts on the metabolite in Human-GEM |
| `metab_gene_adjacency_rows.csv` | 246 rows | row metadata: RefMet class, KEGG/RefMet ids, platform, exercise response stats, link count |
| `metab_gene_adjacency_cols.csv` | 647 rows | column metadata: Ensembl id, GTEx median TPM (blood/muscle/adipose), tissue specificity tau, eGene flags, link count |
| `metab_gene_adjacency_edges.csv` | 1,562 rows | one row per `1` cell: role (substrate/product/either/curated), transport flag, reaction count, subsystem, `gem_degree` |

Row order = metabolites sorted by smallest exercise-vs-control FDR. Column order = order of the source gene list.
Cell density 0.98 %. 83 rows and 79 columns are all zero; dropping them leaves a 163 x 568 connected core.

## 1. Rows: 246 "Tier A" blood metabolites

Source: `BLOOD_METAB_DA` (1,143 named metabolites, 10 platforms) and `HUMAN_FEATURE_TO_GENE` (assay = `metab`).
Rule (all three must hold), implemented in the session script that wrote `results/metab_tierA.csv`:

1. **Has a KEGG compound id** in `HUMAN_FEATURE_TO_GENE.kegg_id`. This is the requirement for a pathway / reaction join. 583 of the 1,143 have any external id; 424 have KEGG.
2. **Not a complex lipid species.** Excluded if the name matches lipid shorthand `^(TG|DG|MG|PC|PE|PI|PS|PG|LPC|LPE|LPI|SM|Cer|HexCer|Hex2Cer|Hex3Cer|GlcCer|CE|CAR|FA|ST)\b` or the RefMet class is one of TG, DG, MG, PC, PE, PI, PS, PG, LPC, LPE, LPI, SM, Cer, HexCer, Hex2Cer, CE, O-PC, O-PE, O-LPC, O-LPE, DHCer, Chol. esters, Acyl carnitines, Unsaturated FA, Saturated FA, Oxidized FA. Reason: sum-composition names such as `TG 54:6` or `PC P-36:3 or PC O-36:4` are unresolved isomer pools with no gene mapping; they are only interpretable at class level. Free fatty acids and acylcarnitines are excluded by this rule too (they belong to the lipid tiers, see below).
3. **Not a contaminant.** Names containing `phthalate` (2 plasticizers) are dropped.

Result: 246 metabolites, 154 of them exercise-responsive (FDR < 0.05 in at least one `exercise_with_controls` contrast).
The same rule applied to muscle and adipose gives 235 and 171; those are in `results/metab_tierA.csv` but are not rows of this matrix because the Human-GEM mapping was built on the blood dictionary.

Tiers that were considered and *not* kept (for reference):
- Tier B, 126: named, RefMet class but no KEGG id (modified amino acids, dipeptides, extra bile acids, hydroxy-FA).
- Tier C1, 112: acylcarnitines and free fatty acids, individually meaningful but lipid-coded.
- Tier C2, 656: complex lipid species, recommended to be collapsed to class nodes.
- Tier D, 2: phthalates.

Note that the "baseline" contrasts (`contrast_type == baseline`, group differences before the bout) were **never** used for the response statistics; only `exercise_with_controls` (EE-CON, RE-CON) contrasts count.

## 2. Columns: 647 GTEx-annotated Human-GEM genes

Source list: `results/metab/gene_list_gtex_gem_tierA.csv`, produced by the metabolism pipeline (`motrpac/metab/genes.py`, `scripts/gtex_metab.py`). Its construction:

1. **Human-GEM enzyme / transporter genes** linked to any metabolite of the harmonised human + rat metabolite dictionary (`results/metab/metabolite_dictionary.parquet`), where the metabolite was mapped to a Human-GEM species through KEGG, ChEBI, HMDB, PubChem or exact name. 1,126 genes at this stage.
2. **Present in the GTEx v8 median-TPM subset** prepared by `scripts/gtex_metab.py`: the union of MSigDB Reactome metabolism sets (names matching METABOLISM, GLYCOLYSIS, GLUCONEOGENESIS, CITRIC_ACID/TCA, FATTY_ACID, BETA_OXIDATION, KETONE, UREA_CYCLE, RESPIRATORY_ELECTRON; 2,738 genes, 97 sets, taken from the MoTrPAC `MOLECULAR_SIGNATURES` object) plus extra GEM genes passed with `--genes`; 3,121 genes have a TPM row. 453 of the 647 are in the Reactome universe, the rest entered through the GEM extras.
3. **eGene in at least one MoTrPAC tissue**: every one of the 647 is a GTEx v8 eGene (q <= 0.05) in whole blood, skeletal muscle or subcutaneous adipose. This holds for all 647 as observed; the writer of the list is not in the repository, so treat it as a property of the list rather than a documented filter.

Attached per gene (from `gtex_features()`): median TPM in blood / muscle / adipose / liver / kidney / heart, dominant tissue, tissue-specificity tau, eGene flag and lead variant per tissue, number of significant eQTL pairs per tissue.

No expression threshold was applied. A stricter variant (median TPM >= 10 in blood, muscle or adipose) would keep 561 columns.

## 3. Cells: Human-GEM reaction membership

A cell is `1` when at least one Human-GEM reaction contains the metabolite and carries the gene in its gene rule
(`motrpac/metab/gem.py`, reaction table parsed from the Human-GEM YAML). Details:

- **Metabolite matching**: Tier A name -> Human-GEM species via KEGG id (176 of the 182 mappable), then ChEBI (3), HMDB (2), exact name (1). 64 Tier A metabolites have no Human-GEM species and are all-zero rows.
- **Currency filter removed.** The pipeline default (`gem.CURRENCY_DEGREE = 60`) drops metabolites in more than 60 reactions. For this matrix that filter was **lifted** at the user's request, which added back 31 metabolites (ATP, ADP, AMP, CMP, FAD, glutathione, carnitine and the proteinogenic amino acids) and raised the non-zero count from 954 to 1,562. `gem_degree` in the edge file records each metabolite's total Human-GEM reaction count so hubs can be down-weighted later.
- **Role**: `substrate` (coef < 0, irreversible), `product` (coef > 0, irreversible), `either` (reversible or mixed across reactions). `is_transport` is true when the metabolite appears in more than one compartment of the reaction or the reaction is an extracellular transport subsystem.
- **Curated additions**: for clinical analytes that are also Tier A names (Glucose) the pipeline's canonical genes (SLC2A4, HK2, G6PC1, PCK1) are added with role `curated` when they are among the 647 columns.
- Compartments are ignored (`met_nocomp`), so a cytosolic and a mitochondrial reaction on the same metabolite count once.

## 4. What the matrix means and does not mean

- An edge is **structural prior knowledge**: the gene encodes an enzyme or transporter acting on the metabolite. It is not derived from MoTrPAC data and carries no exercise information by itself.
- Direction for the graph: Gene -> Metabolite (enzyme acts on metabolite). Sign is not available; use `role` if a producing / consuming distinction is needed.
- Exercise evidence lives on the rows (`responsive`, `min_adj_p`, `peak`) and, for the genes, in the MoTrPAC DA tables via `motrpac.metab.genes.motrpac_gene_da`.
- 45 of the 154 responsive metabolites remain unlinked because they have no Human-GEM species (methylated nucleosides, N-acetyl amino acids, isomers). Extending coverage needs a second reaction source (KEGG RPAIR / Recon3D by KEGG id), not a threshold change.

## 5. Reproduce

```bash
conda activate motrpac
python scripts/metab_pipeline.py            # dictionary, GEM links, GTEx gene list (results/metab/)
# Tier A rows: see the rule in section 1 (session script; output results/metab_tierA.csv)
# Matrix: session script, Human-GEM links rebuilt without the currency filter, restricted to rows x columns
```
