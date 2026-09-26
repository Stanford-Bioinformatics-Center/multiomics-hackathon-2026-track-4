# Rat metabolite x gene adjacency: same pipeline as the human matrix

Built 2026-09-26 from the MoTrPAC 2024 rat endurance-training release (`MotrpacRatTraining6moData`: `METAB_PLASMA_DA`,
`METAB_PLASMA_DA_METAREG`, `TRAINING_REGULATED_FEATURES`), the harmonised metabolite dictionary from the metabolism
pipeline, Human-GEM and the same 647-gene GTEx list used for the human matrix. Read `METAB_GENE_ADJACENCY.md` first;
this file lists only what differs. Transfer bundle: `results/transfer/rat_metab_gene_adjacency_*.csv`.

| File | Shape | Content |
|---|---|---|
| `rat_metab_gene_adjacency_219x647.csv` | 219 x 647 | binary matrix, `1` = gene product acts on the metabolite in Human-GEM |
| `rat_metab_gene_adjacency_rows.csv` | 219 | RefMet name / id, KEGG, main and sub class, platform, training response, link count |
| `rat_metab_gene_adjacency_cols.csv` | 647 | as the human column file, plus `rat_symbol` (MoTrPAC `RAT_TO_HUMAN_GENE` orthology) |
| `rat_metab_gene_adjacency_edges.csv` | 1,186 | one row per `1` cell: role, transport flag, reaction count, subsystem, `gem_degree`, rat gene symbol |
| `metab_human_rat_crosswalk_by_kegg.csv` | 313 | outer join of human and rat Tier A rows on KEGG id: 171 matched pairs over 153 shared KEGG ids (a few ids carry two names on one side), 92 human-only, 50 rat-only; 48 pairs responsive in both species |

Density 0.84 %. 166 rows and 444 columns have at least one link; the connected core is 166 x 444.

## 1. Rows: 219 rat plasma "Tier A" metabolites

Universe: 1,084 named rat plasma metabolites (all sedentary-vs-trained platforms, plus the 112 `meta-reg` features that
the consortium meta-regressed across platforms). Rule, identical to the human one:

1. **KEGG id** present in the harmonised dictionary (262 of 1,084; the dictionary took ids from RefMet, the MoTrPAC feature
   map and Metabolomics Workbench).
2. **Not a complex lipid**: name or RefMet name matching the lipid shorthand regex, or RefMet **sub_class** in the same
   exclusion set as the human matrix (TG, DG, MG, PC, PE, PI, PS, PG, LPC, LPE, LPI, SM, Cer, HexCer, Hex2Cer, CE, O-PC,
   O-PE, O-LPC, O-LPE, DHCer, Chol. esters, Acyl carnitines, Unsaturated FA, Saturated FA, Oxidized FA). 39 KEGG-mapped
   lipids removed. For the human matrix the class came from the MoTrPAC REFMET signature sets; for rat it comes from the
   dictionary's RefMet sub_class. The two vocabularies are the same RefMet sub-class names.
3. **No phthalates** (none present in rat).

Row order: responsive first, then ascending pooled p-value.

Differences from the human rows to keep in mind:
- **Names are not harmonised across species.** Rat feature ids keep the source spelling (`pyruvate`, `vitamin A`,
  `1-methyladenosine`) whereas human ids are RefMet-style (`Pyruvic acid`, `Retinol`, `1-Methyladenosine`). Only 71 row
  names match exactly; **153 match by KEGG id**. Always join the two matrices through `kegg_id`, or use the crosswalk file.
- **Fewer rows** (219 vs 246) because the rat panel has fewer KEGG-annotated small molecules, not because of a stricter rule.

## 2. Training response on the rows (rat, not exercise bout)

The rat study is 1, 2, 4 and 8 weeks of treadmill training versus sedentary controls, sampled 48 h after the last bout,
both sexes, with the consortium's per-sex limma fits. Row statistics use the **pooled** summary already built in
`results/metab/adjacency/R_ET__vs_sedentary_{pooled__log2FC, pooled__p, _edge}.csv`:

- `n_sig_weeks`: weeks (of 1w, 2w, 4w, 8w) flagged as an edge in the pooled table.
- `min_p_pooled`, `max_abs_log2FC`, `peak_week`: over the four weeks.
- `training_regulated`: feature is in the consortium's `TRAINING_REGULATED_FEATURES` (plasma metabolomics, training q < 0.05 after
  IHW). This is the consortium's own call and is the more conservative flag.
- `responsive` = `n_sig_weeks > 0` OR `training_regulated`. 80 of 219 rows; here the two criteria coincide.

These are chronic-training effects, not acute-bout effects. A rat row marked responsive and a human row marked responsive
are answering different questions (adaptation over weeks vs response within hours); the crosswalk marks 153 shared KEGG
compounds so the comparison can be made explicitly.

## 3. Columns and cells

- Columns are the **same 647 human genes**, in the same order, so the human and rat matrices can be stacked or compared
  column-wise. Human-GEM is a human model; rat enzymes are represented through their human orthologs. `rat_symbol` in the
  column file gives the rat gene from MoTrPAC's `RAT_TO_HUMAN_GENE` map; 13 of the 647 have no rat ortholog listed and
  their columns should be treated as human-only.
- Cell rule, mapping routes and the lifted currency filter are exactly as in the human matrix (KEGG 179, ChEBI 2, HMDB 1,
  name 1 of the 183 mappable rows; hubs such as AMP, glucose, glutamate and succinate are kept; `gem_degree` records the
  hub size).
- Unlinked rows: 36 rows have no Human-GEM species (methylated and acetylated amino acids, homoarginine, pipecolic acid,
  5-methylcytosine); a further 17 map to a species whose reactions have no gene in the 647 list. 20 of the 80 responsive
  rows are unlinked for these reasons.

## 4. Summary numbers

| | Human | Rat |
|---|---|---|
| Universe | 1,143 blood | 1,084 plasma |
| Tier A rows | 246 | 219 |
| Responsive rows | 154 (acute bout, FDR < 0.05 vs control) | 80 (training, pooled sexes) |
| Mapped to Human-GEM | 182 | 183 |
| Non-zero cells | 1,562 | 1,186 |
| Rows / columns linked | 163 / 568 | 166 / 444 |
| Top metabolite hubs | ATP 117, ADP 101, GTP 60 | AMP 44, glucose 25, vitamin A 24 |
| Top gene hubs | CTSA 24, TPP1 23, SLC7A6 22 | TPP1 30, CTSA 29, SLC7A6 28 |
