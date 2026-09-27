# References — exercise-kg

Every edge in the graph carries `source_refs` pointing to these IDs. The explorer, the Hypotheses cards and Claude's answers cite the same IDs. Claude is instructed to cite only from this list and to mark anything else `[unsourced]`.

## Datasets

- **D1** — MoTrPAC Study Group. Temporal dynamics of the multi-omic response to endurance exercise training. Nature 629, 174–183 (2024). [doi:10.1038/s41586-023-06877-w](https://doi.org/10.1038/s41586-023-06877-w)  
  *Used for:* Rat 8-week endurance training multi-omics; training-regulated features (training q), per-week state calls, graphical clusters and their pathway enrichment
- **D2** — MoTrPAC Bioinformatics Center. MotrpacRatTraining6moData: data for the MoTrPAC endurance-training study in 6-month-old rats. GitHub, commit f831a4f. [link](https://github.com/MoTrPAC/MotrpacRatTraining6moData)  
  *Used for:* Rat tables used here: TRAINING_REGULATED_FEATURES, *_DA, GRAPH_STATES, GRAPH_PW_ENRICH, FEATURE_TO_GENE, RAT_TO_HUMAN_GENE/PHOSPHO, PHENO
- **D3** — Katz DH, Jin CA, Many GM, et al. Multi-Omic, Multi-Tissue Responses to Acute Exercise in Sedentary Adults: Findings from the Molecular Transducers of Physical Activity Consortium. bioRxiv (2026). [doi:10.64898/2026.02.27.702183](https://doi.org/10.64898/2026.02.27.702183)  
  *Used for:* Human acute endurance/resistance vs control; blood, muscle, adipose multi-omics (pre-suspension cohort)
- **D4** — Keshishian H, Many GM, Smith G, et al. Integrative Multi-omics Analysis of the Human Skeletal Muscle Response to Endurance or Resistance Exercise: Findings from MoTrPAC. bioRxiv (2026). [doi:10.64898/2026.03.04.705181](https://doi.org/10.64898/2026.03.04.705181)  
  *Used for:* Human skeletal muscle acute responses (transcriptome, proteome, phosphoproteome, metabolome)
- **D5** — Ahn C, Jin CA, Whytock KL, et al.; MoTrPAC Study Group. Multi-omic responses to acute exercise in abdominal subcutaneous adipose tissue of sedentary adults: findings from MoTrPAC. bioRxiv (2026). [doi:10.64898/2026.03.05.702363](https://doi.org/10.64898/2026.03.05.702363)  
  *Used for:* Human subcutaneous adipose acute responses
- **D6** — Robbins JM, Katz DH, Many GM, et al. Blood Biochemical Responses to Acute Exercise: Findings from MoTrPAC. bioRxiv (2026). [doi:10.64898/2026.03.02.704798](https://doi.org/10.64898/2026.03.02.704798)  
  *Used for:* Human blood acute responses (Olink plasma proteomics, plasma metabolomics, whole-blood RNA-seq)
- **D7** — Brandt AR, Fleg J, Goodpaster BH, et al. Molecular Transducers of Physical Activity Consortium (MoTrPAC): Initial Insights into the Dynamic Human Responses to Exercise. bioRxiv (2026). [doi:10.64898/2026.03.02.705347](https://doi.org/10.64898/2026.03.02.705347)  
  *Used for:* Human study design, randomisation and sampling timeline
- **D8** — MoTrPAC Bioinformatics Center. MotrpacHumanPreSuspensionAnalysis: public summary statistics, differential analysis, enrichment and clustering for the human pre-suspension acute-exercise cohort (collection c2.0). GitHub. [link](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis)  
  *Used for:* Human tables used here: *_DA (control-adjusted contrasts, BH adj p), CAMERA_RESULTS, PTMSEA_RESULTS, FCM_CLUSTERS, HUMAN_FEATURE_TO_GENE, MOLECULAR_SIGNATURES, clinical chemistry
- **D9** — Sanford JA, Nogiec CD, Lindholm ME, et al. Molecular Transducers of Physical Activity Consortium (MoTrPAC): Mapping the Dynamic Responses to Exercise. Cell 181, 1464–1474 (2020). [doi:10.1016/j.cell.2020.06.004](https://doi.org/10.1016/j.cell.2020.06.004)  
  *Used for:* MoTrPAC consortium design

## Knowledge resources

- **R1** — Szklarczyk D, Kirsch R, Koutrouli M, et al. The STRING database in 2023: protein–protein association networks and functional enrichment analyses for any sequenced genome of interest. Nucleic Acids Res 51, D638–D646 (2023). [doi:10.1093/nar/gkac1000](https://doi.org/10.1093/nar/gkac1000)  
  *Used for:* Physical PPI subnetwork, combined score >= 700
- **R2** — Himmelstein DS, Lizee A, Hessler C, et al. Systematic integration of biomedical knowledge prioritizes drugs for repurposing. eLife 6, e26726 (2017). [doi:10.7554/eLife.26726](https://doi.org/10.7554/eLife.26726)  
  *Used for:* Curated human PPI (Gene-interacts-Gene)
- **R3** — Rolland T, Taşan M, Charloteaux B, et al. A proteome-scale map of the human interactome network. Cell 159, 1212–1226 (2014). [doi:10.1016/j.cell.2014.10.050](https://doi.org/10.1016/j.cell.2014.10.050)  
  *Used for:* Binary interactome underlying Hetionet PPI
- **R4** — Dimitrov D, Türei D, Garrido-Rodriguez M, et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022). [doi:10.1038/s41467-022-30755-0](https://doi.org/10.1038/s41467-022-30755-0)  
  *Used for:* Ligand-receptor consensus resource
- **R5** — Türei D, Valdeolivas A, Gul L, et al. Integrated intra- and intercellular signaling knowledge for multicellular omics analysis. Mol Syst Biol 17, e9923 (2021). [doi:10.15252/msb.20209923](https://doi.org/10.15252/msb.20209923)  
  *Used for:* Source of the LIANA consensus ligand-receptor pairs
- **R6** — Hornbeck PV, Zhang B, Murray B, et al. PhosphoSitePlus, 2014: mutations, PTMs and recalibrations. Nucleic Acids Res 43, D512–D520 (2015). [doi:10.1093/nar/gku1267](https://doi.org/10.1093/nar/gku1267)  
  *Used for:* Kinase-substrate relationships
- **R7** — Krug K, Mertins P, Zhang B, et al. A curated resource for phosphosite-specific signature analysis. Mol Cell Proteomics 18, 576–593 (2019). [doi:10.1074/mcp.TIR118.000943](https://doi.org/10.1074/mcp.TIR118.000943)  
  *Used for:* PTM signatures and PTM-SEA enrichment
- **R8** — Liberzon A, Subramanian A, Pinchback R, et al. Molecular signatures database (MSigDB) 3.0. Bioinformatics 27, 1739–1740 (2011). [doi:10.1093/bioinformatics/btr260](https://doi.org/10.1093/bioinformatics/btr260)  
  *Used for:* Curated gene-set collections bundled with the human package
- **R9** — Milacic M, Beavers D, Conley P, et al. The Reactome Pathway Knowledgebase 2024. Nucleic Acids Res 52, D672–D678 (2024). [doi:10.1093/nar/gkad1025](https://doi.org/10.1093/nar/gkad1025)  
  *Used for:* Reactome pathways
- **R10** — Agrawal A, Balcı H, Hanspers K, et al. WikiPathways 2024: next generation pathway database. Nucleic Acids Res 52, D679–D689 (2024). [doi:10.1093/nar/gkad960](https://doi.org/10.1093/nar/gkad960)  
  *Used for:* WikiPathways
- **R11** — Kanehisa M, Furumichi M, Sato Y, et al. KEGG for taxonomy-based analysis of pathways and genomes. Nucleic Acids Res 51, D587–D592 (2023). [doi:10.1093/nar/gkac963](https://doi.org/10.1093/nar/gkac963)  
  *Used for:* KEGG / KEGG MEDICUS pathways
- **R12** — Schaefer CF, Anthony K, Krupa S, et al. PID: the Pathway Interaction Database. Nucleic Acids Res 37, D674–D679 (2009). [doi:10.1093/nar/gkn653](https://doi.org/10.1093/nar/gkn653)  
  *Used for:* PID pathways
- **R13** — Rath S, Sharma R, Gupta R, et al. MitoCarta3.0: an updated mitochondrial proteome now with sub-organelle localization and pathway annotations. Nucleic Acids Res 49, D1541–D1547 (2021). [doi:10.1093/nar/gkaa1011](https://doi.org/10.1093/nar/gkaa1011)  
  *Used for:* MitoCarta pathways
- **R14** — Gene Ontology Consortium; Aleksander SA, Balhoff J, Carbon S, et al. The Gene Ontology knowledgebase in 2023. Genetics 224, iyad031 (2023). [doi:10.1093/genetics/iyad031](https://doi.org/10.1093/genetics/iyad031)  
  *Used for:* GO terms
- **R15** — Fahy E, Subramaniam S. RefMet: a reference nomenclature for metabolomics. Nat Methods 17, 1173–1174 (2020). [doi:10.1038/s41592-020-01009-y](https://doi.org/10.1038/s41592-020-01009-y)  
  *Used for:* Metabolite names shared across species
- **R16** — Mungall CJ, Torniai C, Gkoutos GV, et al. Uberon, an integrative multi-species anatomy ontology. Genome Biol 13, R5 (2012). [doi:10.1186/gb-2012-13-1-r5](https://doi.org/10.1186/gb-2012-13-1-r5)  
  *Used for:* Cross-species tissue mapping
- **R17** — Bansal P, et al. Rhea, the reaction knowledgebase in 2022. Nucleic Acids Res 50(D1), D693–D700 (2022). [doi:10.1093/nar/gkab1016](https://doi.org/10.1093/nar/gkab1016)  
  *Used for:* Enzyme reactions linking proteins to metabolites (test step for co-regulation hypotheses)
- **R18** — Wishart DS, et al. HMDB 5.0: the Human Metabolome Database for 2022. Nucleic Acids Res 50(D1), D622–D631 (2022). [doi:10.1093/nar/gkab1062](https://doi.org/10.1093/nar/gkab1062)  
  *Used for:* Metabolite–enzyme and transporter annotations (test step for co-regulation hypotheses)
- **R19** — Schriml LM, et al. Human Disease Ontology 2018 update: classification, content and workflow expansion. Nucleic Acids Res 47(D1), D955–D962 (2019). [doi:10.1093/nar/gky1032](https://doi.org/10.1093/nar/gky1032)  
  *Used for:* Disease nodes (DOID) in the Hetionet disease layer
- **R20** — Wishart DS, et al. DrugBank 5.0: a major update to the DrugBank database for 2018. Nucleic Acids Res 46(D1), D1074–D1082 (2018). [doi:10.1093/nar/gkx1037](https://doi.org/10.1093/nar/gkx1037)  
  *Used for:* Compound nodes and drug–target / drug–disease links (via Hetionet)

## Statistical methods

- **M1** — Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc B 57, 289–300 (1995). [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)  
  *Used for:* BH adjusted p-values
- **M2** — Hoffman GE, Roussos P. Dream: powerful differential expression analysis for repeated measures designs. Bioinformatics 37, 192–201 (2021). [doi:10.1093/bioinformatics/btaa687](https://doi.org/10.1093/bioinformatics/btaa687)  
  *Used for:* Human mixed-model differential analysis (dream-acute)
- **M3** — Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene correlation. Nucleic Acids Res 40, e133 (2012). [doi:10.1093/nar/gks461](https://doi.org/10.1093/nar/gks461)  
  *Used for:* Human pathway enrichment
- **M4** — Kolberg L, Raudvere U, Kuzmin I, et al. g:Profiler—interoperable web service for functional enrichment analysis and gene identifier mapping (2023 update). Nucleic Acids Res 51, W207–W212 (2023). [doi:10.1093/nar/gkad347](https://doi.org/10.1093/nar/gkad347)  
  *Used for:* Rat cluster pathway enrichment
- **M5** — Welch BL. The generalization of 'Student's' problem when several different population variances are involved. Biometrika 34, 28–35 (1947). [doi:10.1093/biomet/34.1-2.28](https://doi.org/10.1093/biomet/34.1-2.28)  
  *Used for:* Rat phenotype trained-vs-control test (computed in this project)
- **M6** — Fisher RA. On the interpretation of χ2 from contingency tables, and the calculation of P. J R Stat Soc 85(1), 87–94 (1922). [doi:10.1111/j.2397-2335.1922.tb00768.x](https://doi.org/10.1111/j.2397-2335.1922.tb00768.x)  
  *Used for:* Hypergeometric over-representation of disease genes among exercise-responsive genes (computed in this project)

## Literature cited for interpretations and hypotheses

- **L1** — Egan B, Zierath JR. Exercise metabolism and the molecular regulation of skeletal muscle adaptation. Cell Metab 17, 162–184 (2013). [doi:10.1016/j.cmet.2012.12.012](https://doi.org/10.1016/j.cmet.2012.12.012)  
  *Used for:* Repeated acute exercise signals accumulate into training adaptation
- **L2** — Perry CGR, Lally J, Holloway GP, et al. Repeated transient mRNA bursts precede increases in transcriptional and mitochondrial proteins during training in human skeletal muscle. J Physiol 588, 4795–4810 (2010). [doi:10.1113/jphysiol.2010.199448](https://doi.org/10.1113/jphysiol.2010.199448)  
  *Used for:* Transient acute transcript bursts vs sustained protein changes with training
- **L3** — Pillon NJ, Gabriel BM, Dollet L, et al. Transcriptomic profiling of skeletal muscle adaptations to exercise and inactivity. Nat Commun 11, 470 (2020). [doi:10.1038/s41467-019-13869-w](https://doi.org/10.1038/s41467-019-13869-w)  
  *Used for:* Acute vs training muscle transcriptome responses differ
- **L4** — Chow LS, Gerszten RE, Taylor JM, et al. Exerkines in health, resilience and disease. Nat Rev Endocrinol 18, 273–289 (2022). [doi:10.1038/s41574-022-00641-2](https://doi.org/10.1038/s41574-022-00641-2)  
  *Used for:* Exercise-induced secreted signals mediating inter-organ communication
- **L5** — Hoffman NJ, Parker BL, Chaudhuri R, et al. Global phosphoproteomic analysis of human skeletal muscle reveals a network of exercise-regulated kinases and AMPK substrates. Cell Metab 22, 922–935 (2015). [doi:10.1016/j.cmet.2015.09.001](https://doi.org/10.1016/j.cmet.2015.09.001)  
  *Used for:* Acute exercise signals largely through phosphorylation in human muscle
- **L6** — Casado P, Rodriguez-Prados JC, Cosulich SC, et al. Kinase-substrate enrichment analysis provides insights into the heterogeneity of signaling pathway activation in leukemia cells. Sci Signal 6, rs6 (2013). [doi:10.1126/scisignal.2003573](https://doi.org/10.1126/scisignal.2003573)  
  *Used for:* Inferring kinase activity from substrate phosphosite changes
- **L7** — Barabási AL, Gulbahce N, Loscalzo J. Network medicine: a network-based approach to human disease. Nat Rev Genet 12, 56–68 (2011). [doi:10.1038/nrg2918](https://doi.org/10.1038/nrg2918)  
  *Used for:* Interacting, co-regulated proteins form functional modules; hubs
- **L8** — Brooks GA. The science and translation of lactate shuttle theory. Cell Metab 27, 757–785 (2018). [doi:10.1016/j.cmet.2018.03.008](https://doi.org/10.1016/j.cmet.2018.03.008)  
  *Used for:* Lactate production and exchange between organs during exercise
- **L9** — Horowitz JF, Klein S. Lipid metabolism during endurance exercise. Am J Clin Nutr 72(2 Suppl), 558S–563S (2000). [link](https://pubmed.ncbi.nlm.nih.gov/10919960/)  
  *Used for:* Adipose lipolysis raises plasma glycerol and fatty acids during endurance exercise
- **L10** — Liu Y, Beyer A, Aebersold R. On the dependency of cellular protein levels on mRNA abundance. Cell 165(3), 535–550 (2016). [doi:10.1016/j.cell.2016.03.014](https://doi.org/10.1016/j.cell.2016.03.014)  
  *Used for:* mRNA and protein levels are often decoupled; basis for cross-omics agreement/disagreement cards
- **L11** — Contrepois K, et al. Molecular choreography of acute exercise. Cell 181(5), 1112–1130 (2020). [doi:10.1016/j.cell.2020.04.043](https://doi.org/10.1016/j.cell.2020.04.043)  
  *Used for:* Acute exercise drives large, transient multi-omic shifts in human blood
- **L12** — Hernandez-Armenta C, Ochoa D, Gonçalves E, Saez-Rodriguez J, Beltrao P. Benchmarking substrate-based kinase activity inference using phosphoproteomic data. Bioinformatics 33(12), 1845–1851 (2017). [doi:10.1093/bioinformatics/btx082](https://doi.org/10.1093/bioinformatics/btx082)  
  *Used for:* Kinase activity can be inferred from changes in known substrate sites
- **L13** — Holloszy JO, Coyle EF. Adaptations of skeletal muscle to endurance exercise and their metabolic consequences. J Appl Physiol 56(4), 831–838 (1984). [doi:10.1152/jappl.1984.56.4.831](https://doi.org/10.1152/jappl.1984.56.4.831)  
  *Used for:* Muscle mitochondrial and metabolic adaptation underlies endurance-training gains
- **L14** — Nikolić N, et al. Electrical pulse stimulation of cultured skeletal muscle cells as a model for in vitro exercise – possibilities and limitations. Acta Physiol 220(3), 310–331 (2017). [doi:10.1111/apha.12830](https://doi.org/10.1111/apha.12830)  
  *Used for:* In vitro exercise model (EPS) proposed in hub tests
- **L15** — Narkar VA, et al. AMPK and PPARδ agonists are exercise mimetics. Cell 134(3), 405–415 (2008). [doi:10.1016/j.cell.2008.06.051](https://doi.org/10.1016/j.cell.2008.06.051)  
  *Used for:* AICAR as an exercise mimetic proposed in hub tests
- **L16** — Pedersen BK, Saltin B. Exercise as medicine – evidence for prescribing exercise as therapy in 26 different chronic diseases. Scand J Med Sci Sports 25(Suppl 3), 1–72 (2015). [doi:10.1111/sms.12581](https://doi.org/10.1111/sms.12581)  
  *Used for:* Exercise as therapy for chronic diseases; basis for 'exercise opposes a disease signature' hypotheses
- **L17** — Nieman DC, Wentz LM. The compelling link between physical activity and the body's defense system. J Sport Health Sci 8(3), 201–217 (2019). [doi:10.1016/j.jshs.2018.09.009](https://doi.org/10.1016/j.jshs.2018.09.009)  
  *Used for:* Acute exercise transiently mobilises immune cells and inflammatory signals; basis for 'one bout mimics a disease signature in blood'

## This project

- **P1** — Exercise KG build, this project: percentile weights, metabolite co-regulation (Pearson |r| >= 0.9 over response profiles), rule-based hypothesis generation.   
  *Used for:* Derived scores and correlations computed here (not peer reviewed)

