# Exercise x metabolite adjacency: contrasts, values and rules

Built 2026-09-26 from `BLOOD_METAB_DA` (MoTrPAC 2026 human acute-exercise release). Rows are the same 246 Tier A blood
metabolites, in the same order, as the metabolite x gene matrix (see `METAB_GENE_ADJACENCY.md`), so the two matrices chain:
Exercise -> Metabolite -> Gene. Transfer bundle: `results/transfer/exercise_metab_*.csv`.

| File | Shape | Content |
|---|---|---|
| `exercise_metab_246x10_weighted_log2FC_fdr05.csv` | 246 x 10 | **primary edge matrix**: signed weight = log2 fold change where FDR < 0.05, else `0` |
| `exercise_metab_246x10_weighted_z_fdr05.csv` | 246 x 10 | alternative weight: z statistic where FDR < 0.05, else `0` (significance-scaled) |
| `exercise_metab_246x10_p_nominal.csv` | 246 x 10 | unadjusted p, in case a nominal threshold is preferred |
| `exercise_metab_246x10_log2FC.csv` | 246 x 10 | effect size, log2 fold change of the delta-delta contrast |
| `exercise_metab_246x10_z.csv` | 246 x 10 | standardized z statistic (`z.std`), comparable across features |
| `exercise_metab_246x10_fdr.csv` | 246 x 10 | Benjamini-Hochberg adjusted p (`adj_p_value`) |
| `exercise_metab_246x33_z_allcontrasts.csv` | 246 x 33 | z for **every** contrast in the release, for co-response correlation only |
| `exercise_metab_cols.csv` | 33 rows | column metadata: contrast type, category, group, timepoint, full contrast string, n per arm, `causal`, `in_primary_10` |

Row-label header is `metabolite`; column labels are `<contrast_category>|<Timepoint>`, e.g. `EE-CON|post_10_min`.

## 1. Columns of the primary matrix: 10 exercise-vs-control contrasts

Only `contrast_type == exercise_with_controls` is used. Each column is a **difference-in-change** estimate from the
consortium's dream mixed model (`~ 0 + group_timepoint + Sex + calculatedAge + BMI + site + (1|pid)`):

    (exercise arm at timepoint  -  exercise arm at pre)  -  (control arm at timepoint  -  control arm at pre)

so the time-of-day, fasting and repeated-sampling effects that the control arm also experiences are removed.
Because arms were randomized, the column is a **causal, population-level** effect of one bout of that exercise type.

| Column | Exercise arm n | Control arm n |
|---|---|---|
| `EE-CON\|during_20_min` | 63 | 36 |
| `EE-CON\|during_40_min` | 58 | 37 |
| `EE-CON\|post_10_min` | 64 | 34 |
| `EE-CON\|post_15_30_45_min` | 64 | 37 |
| `EE-CON\|post_3.5_4_hr` | 28 | 16 |
| `EE-CON\|post_24_hr` | 33 | 19 |
| `RE-CON\|post_10_min` | 68 | 34 |
| `RE-CON\|post_15_30_45_min` | 72 | 37 |
| `RE-CON\|post_3.5_4_hr` | 35 | 16 |
| `RE-CON\|post_24_hr` | 38 | 19 |

EE = endurance bout (5 min warm-up, 25-30 min cycle, 25-30 min treadmill); RE = resistance bout (8 exercises, 3 x 10RM).
Resistance has no during-bout draws, hence 6 EE and 4 RE columns. Pre-exercise is the reference and is not a column.
The 3.5-4 h and 24 h cells rest on roughly half the cohort; their few edges reflect lower power as well as recovery.

## 2. Cell rule

- **Weighted, not binary.** `weight = log2FC` if `adj_p_value < 0.05`, else `0`. Sign is the direction of change, magnitude
  is the effect size in log2 units (delta-delta vs control). FDR is the consortium's, computed per tissue x ome x contrast.
- The z-weighted file uses the same mask with `z.std` as the value; use it when edges should be scaled by confidence
  rather than by effect size (z folds in the standard error, so a small but precise effect scores higher).
- Threshold is the **adjusted** p. A nominal p < 0.05 would keep 723 cells instead of 501; `exercise_metab_246x10_p_nominal.csv`
  is provided so that choice can be made downstream.
- No fold-change floor is applied. 142 of the 501 kept cells have |log2FC| < 0.2 (significant but small); the older
  1,143-row tables in `results/metab/adjacency/` additionally required |log2FC| >= log2(1.2).
- Values come straight from the consortium tables; nothing was refit.

## 3. What is in the result

- 501 non-zero weights in 246 x 10 (20 % dense): 315 positive, 186 negative. 154 metabolites have at least one edge, 92 have none.
- |weight| ranges 0.05 to 4.85 log2 units, median 0.32. The largest are inosine and hypoxanthine after resistance exercise
  (4.1 to 4.85 at 10 and 15-45 min), xanthosine (2.3) and 11-deoxycortisol after endurance (2.1).
- Edges concentrate in the first 45 minutes: 63 / 72 / 71 / 55 for the four early EE columns, 92 / 92 for early RE, then 35 and 17 for RE at 3.5 h and 24 h, and
  only 1-3 per column at 3.5 h and 24 h for EE.
- 26 metabolites respond to endurance only, 40 to resistance only, 88 to both. 19 change direction across the time course
  (e.g. ketone-related and fatty-acid intermediates suppressed during the bout, elevated at 3.5 h).
- Strongest and most persistent rows: purine breakdown (hypoxanthine, xanthine, xanthosine, inosine), branched-chain
  keto-acids (ketoleucine, 3-methyl-2-oxovaleric acid), valine, serine, GABA, pregnenolone sulfate, lactate, TCA acids.

## 4. The 33-column z matrix: use only for co-response

`exercise_metab_246x33_z_allcontrasts.csv` adds the other four contrast types:

| contrast_type | Columns | Meaning | Causal? |
|---|---|---|---|
| `exercise_with_controls` | 10 | delta-delta vs control (above) | yes |
| `exercise_no_controls` | 10 | within-arm change from pre (EE-EE, RE-RE) | no, confounded with day effects |
| `control_only` | 6 | control arm change from pre (CON-CON) | isolates the day effect |
| `Endur_vs_Resist` | 4 | EE change minus RE change | modality contrast, no control |
| `baseline` | 3 | arm differences before the bout | should be null; a randomization / batch check |

These columns are linear combinations of each other (EE-CON = EE-EE minus CON-CON), so the 33 are not independent
observations. They exist because 10 points are too few for metabolite-metabolite profile correlation (a |r| > 0.9 pair is
barely above the noise floor with 10 columns, and 3 to 4 orders of magnitude rarer than noise with all 33). Do **not**
read the non-causal columns as exercise effects, and build permutation nulls by shuffling the actual columns.

## 5. Interpretation in the graph

- Node `Exercise:EE` or `Exercise:RE` -> node `Metabolite`, one edge per non-zero cell, weight = log2FC (or z), attributes:
  timepoint, FDR, n per arm. Direction is from the intervention to the molecule; the weight's sign is the direction of change.
- Combined with the metabolite x gene matrix this yields the path Exercise -> Metabolite <- Gene, where the second hop is
  prior knowledge (Human-GEM), not an exercise result. A gene's own exercise response has to be added from the transcript
  or protein DA tables separately.
- The nine clinical analytes (lactate, glucose, insulin, ...) are not rows here; they live in `BLOOD_*_CLINICAL_DA`. Lactate
  and glucose appear as their mass-spec counterparts `Lactic acid` and `Glucose` (log2FC agreement with the clinical
  assay r = 0.998 across contrasts, see `results/metab/summary.md`).

## 6. Reproduce

Session script; inputs `data/cache/human/BLOOD_METAB_DA.parquet`, `BLOOD_METAB_SUM_STATS.parquet` (for n) and the row
order of `results/metab/adjacency/tierA_metab_x_gtex_gene__binary.csv`.
