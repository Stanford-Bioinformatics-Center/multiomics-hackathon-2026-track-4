// Run once after import. Indexes make look-ups by id / symbol / name instant.
CREATE INDEX gene_id        IF NOT EXISTS FOR (n:Gene)             ON (n.id);
CREATE INDEX gene_symbol    IF NOT EXISTS FOR (n:Gene)             ON (n.symbol);
CREATE INDEX feature_id     IF NOT EXISTS FOR (n:MolecularFeature) ON (n.id);
CREATE INDEX feature_src    IF NOT EXISTS FOR (n:MolecularFeature) ON (n.source_id);
CREATE INDEX metab_name     IF NOT EXISTS FOR (n:Metabolite)       ON (n.refmet_name);
CREATE INDEX pathway_id     IF NOT EXISTS FOR (n:Pathway)          ON (n.id);
CREATE INDEX pathway_name   IF NOT EXISTS FOR (n:Pathway)          ON (n.name);
CREATE INDEX tissue_id      IF NOT EXISTS FOR (n:Tissue)           ON (n.id);
CREATE INDEX group_id       IF NOT EXISTS FOR (n:ExerciseGroup)    ON (n.id);
CREATE INDEX pheno_id       IF NOT EXISTS FOR (n:Phenotype)        ON (n.id);
CREATE FULLTEXT INDEX pathway_text IF NOT EXISTS FOR (n:Pathway) ON EACH [n.name];
CREATE INDEX disease_id     IF NOT EXISTS FOR (n:Disease)          ON (n.id);
CREATE INDEX disease_name   IF NOT EXISTS FOR (n:Disease)          ON (n.name);
CREATE INDEX compound_id    IF NOT EXISTS FOR (n:Compound)         ON (n.id);
CREATE INDEX compound_name  IF NOT EXISTS FOR (n:Compound)         ON (n.name);
