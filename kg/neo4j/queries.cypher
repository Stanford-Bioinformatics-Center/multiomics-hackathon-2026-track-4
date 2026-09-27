// Example queries for the significant-response graph (kg/exports_sig). Paste into Neo4j Browser.

// 1. Overview: node and relationship counts
MATCH (n) RETURN labels(n)[-1] AS label, count(*) AS n ORDER BY n DESC;
MATCH ()-[r]->() RETURN type(r) AS type, count(*) AS n ORDER BY n DESC;

// 2. Oxidative phosphorylation in muscle: genes, their molecules and responses (draws as a graph)
MATCH (p:Pathway {id:"WP:WP_OXIDATIVE_PHOSPHORYLATION"})<-[:IN_PATHWAY]-(g:Gene)<-[:MAPS_TO_GENE]-(f)-[r:UPREGULATED_IN|DOWNREGULATED_IN]->(t:Tissue)
WHERE t.id IN ["rat:SKM-GN","rat:SKM-VL","human:muscle"] AND abs(r.z) >= 3
RETURN p, g, f, r, t LIMIT 300;

// 3. Conserved responders: rat 8-week training AND human acute exercise, same direction, in muscle
MATCH (rf)-[rr]->(:Tissue {id:"rat:SKM-GN"}), (rf)-[:MAPS_TO_GENE]->(rg:Gene)-[:ORTHOLOG_OF]->(hg:Gene)<-[:MAPS_TO_GENE]-(hf)-[hr]->(:Tissue {id:"human:muscle"})
WHERE rr.timepoint = "8w" AND type(rr) = type(hr) AND type(rr) IN ["UPREGULATED_IN","DOWNREGULATED_IN"]
RETURN rg.symbol AS rat_gene, hg.symbol AS human_gene, type(rr) AS direction,
       collect(DISTINCT rf.assay)[..3] AS rat_assays, collect(DISTINCT hr.group + " " + hr.timepoint)[..4] AS human_contrasts,
       max(abs(rr.z)) * max(abs(hr.z)) AS score
ORDER BY score DESC LIMIT 25;

// 4. Opposite-direction pairs (acute up, chronic down or vice versa)
MATCH (rf)-[rr]->(:Tissue {id:"rat:SKM-GN"}), (rf)-[:MAPS_TO_GENE]->(rg:Gene)-[:ORTHOLOG_OF]->(hg:Gene)<-[:MAPS_TO_GENE]-(hf)-[hr]->(:Tissue {id:"human:muscle"})
WHERE type(rr) <> type(hr)
RETURN rg.symbol, hg.symbol, type(rr) AS rat, type(hr) AS human, count(*) AS evidence ORDER BY evidence DESC LIMIT 25;

// 5. A metabolite's world: organs where it changes + co-regulated metabolites/proteins
MATCH (m:Metabolite) WHERE toLower(m.refmet_name) = "lactic acid"
OPTIONAL MATCH (m)-[r:UPREGULATED_IN|DOWNREGULATED_IN]->(t:Tissue)
OPTIONAL MATCH (m)-[c:CO_REGULATED_WITH]-(x)
RETURN m, r, t, c, x LIMIT 200;

// 6. Pathways enriched up in both species (rat any week, human any acute contrast)
MATCH (p:Pathway)-[a:ENRICHED_UP_IN]->(:ExerciseGroup {id:"rat:endurance_training"}),
      (p)-[b:ENRICHED_UP_IN]->(h:ExerciseGroup) WHERE h.species = "human"
RETURN p.name, collect(DISTINCT a.tissue) AS rat_tissues, collect(DISTINCT b.tissue + " " + b.timepoint)[..5] AS human, count(*) AS n
ORDER BY n DESC LIMIT 25;

// 7. Phenotypes and when they change
MATCH (ph:Phenotype)-[r:INCREASED_IN|DECREASED_IN]->(g:ExerciseGroup)
RETURN ph.name, g.id, type(r), r.sex, r.timepoint, r.adj_p_value, r.difference AS rat_change ORDER BY ph.name, r.timepoint;

// 8. Kinase signatures hit by significant phosphosites (human direct + rat via orthologous site)
MATCH (s:PTMSite)-[m:IN_PATHWAY]->(p:Pathway) WHERE p.database IN ["PTMSIGDB","PSP"]
RETURN p.name, count(DISTINCT s) AS sites, collect(DISTINCT m.species) AS species ORDER BY sites DESC LIMIT 20;

// ---------------- Interaction layers (INTERACTS_WITH, LIGAND_OF, PHOSPHORYLATES) ----------------

// 9. Interaction hubs among conserved muscle responders (human PPI; rat via orthologs)
MATCH (rg:Gene)-[:ORTHOLOG_OF]->(h:Gene)<-[:MAPS_TO_GENE]-(hf)-[:UPREGULATED_IN|DOWNREGULATED_IN]->(:Tissue {id:"human:muscle"}),
      (rg)<-[:MAPS_TO_GENE]-(rf)-[:UPREGULATED_IN|DOWNREGULATED_IN]->(:Tissue {id:"rat:SKM-GN"})
WITH DISTINCT h
MATCH (h)-[:INTERACTS_WITH]-(p:Gene) WHERE p.species = "human"
WITH h, collect(DISTINCT p.symbol) AS partners
RETURN h.symbol, size(partners) AS degree, partners[..10] ORDER BY degree DESC LIMIT 20;

// 10. Kinases whose substrate sites change (activity proxy), by species
MATCH (k:Gene)-[p:PHOSPHORYLATES]->(s:PTMSite)-[r:UPREGULATED_IN|DOWNREGULATED_IN]->(t:Tissue)
RETURN k.symbol AS kinase, s.species AS site_species, count(DISTINCT s) AS sites,
       sum(CASE type(r) WHEN "UPREGULATED_IN" THEN 1 ELSE 0 END) AS up_events,
       sum(CASE type(r) WHEN "DOWNREGULATED_IN" THEN 1 ELSE 0 END) AS down_events
ORDER BY sites DESC LIMIT 25;

// 11. Candidate inter-organ axes: ligand changes in one organ, its receptor in another
MATCH (lf)-[:MAPS_TO_GENE]->(l:Gene)-[:LIGAND_OF]->(rc:Gene)<-[:MAPS_TO_GENE]-(rf),
      (lf)-[a:UPREGULATED_IN|DOWNREGULATED_IN]->(t1:Tissue), (rf)-[b:UPREGULATED_IN|DOWNREGULATED_IN]->(t2:Tissue)
WHERE t1 <> t2 AND a.species = b.species
RETURN l.symbol AS ligand, t1.name AS ligand_organ, type(a) AS ligand_change, rc.symbol AS receptor, t2.name AS receptor_organ, type(b) AS receptor_change, a.species
LIMIT 50;

// 12. A gene's interaction neighbourhood with responses (draws as a graph)
MATCH (g:Gene {symbol:"PPARGC1A"})-[i:INTERACTS_WITH|LIGAND_OF]-(p:Gene)<-[:MAPS_TO_GENE]-(f)-[r:UPREGULATED_IN|DOWNREGULATED_IN]->(t:Tissue)
RETURN g, i, p, f, r, t LIMIT 200;

// ---------------- Five edges per node, and sources ----------------

// 13. The 5-strongest-edges view around a pathway (ranks are stored on every edge)
MATCH (p:Pathway {id:"WP:WP_OXIDATIVE_PHOSPHORYLATION"})<-[m:IN_PATHWAY]-(g:Gene)<-[f:MAPS_TO_GENE]-(x)-[r:UPREGULATED_IN|DOWNREGULATED_IN]->(t:Tissue)
WHERE m.rank_at_dst <= 5 AND f.rank_at_dst <= 5 AND r.rank_from_src <= 5 AND r.rank_at_dst <= 5
RETURN p, m, g, f, x, r, t;

// 14. Full citations for any edge: split source_refs and join to Reference nodes
MATCH (x)-[r:UPREGULATED_IN]->(t:Tissue {id:"human:muscle"}) WITH r LIMIT 1
UNWIND split(r.source_refs, ";") AS rid
MATCH (ref:Reference {id: rid})
RETURN rid, ref.citation, ref.url;

// ---------------- Disease & drug layer (Hetionet) ----------------

// 15. Organs where a disease's genes are over-represented among exercise-responsive genes
MATCH (d:Disease {name:'type 2 diabetes mellitus'})-[r:DISEASE_GENES_ENRICHED_IN]->(t:Tissue)
RETURN t.id AS organ, r.overlap AS genes, r.fold_enrichment AS fold, r.adj_p_value AS adj_p, r.n_up AS up, r.n_down AS down
ORDER BY adj_p;

// 16. Where exercise opposes (or mimics) a disease's expression signature
MATCH (d:Disease)-[r:EXERCISE_OPPOSES|EXERCISE_MIMICS]->(t:Tissue)
RETURN d.name AS disease, type(r) AS effect, t.id AS organ, r.fraction_opposite AS frac_opposite, r.genes_compared AS n, r.adj_p_value AS adj_p
ORDER BY adj_p LIMIT 30;

// 17. Drugs for a disease whose targets respond to exercise, with the response
MATCH (c:Compound)-[:TREATS]->(d:Disease {name:'coronary artery disease'}), (c)-[:BINDS]->(g:Gene)<-[:MAPS_TO_GENE]-(f)-[r:UPREGULATED_IN|DOWNREGULATED_IN]->(t:Tissue)
RETURN c.name AS drug, g.symbol AS target, type(r) AS response, t.id AS organ, min(r.adj_p_value) AS best_adj_p
ORDER BY best_adj_p LIMIT 25;
