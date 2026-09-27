//node constraints
//CREATE CONSTRAINT constraint_name FOR (n:Label) REQUIRE n.property IS UNIQUE
//Tissue
CREATE CONSTRAINT tissue_name
FOR (n:Tissue)
REQUIRE (n.name, n.species) IS UNIQUE;

//Gene
CREATE CONSTRAINT gene_ncbi
FOR (n:Gene)
REQUIRE n.ncbi IS UNIQUE;
//CREATE CONSTRAINT gene_ncbi_id FOR (n:Gene) REQUIRE n.ncbi IS :: INTEGER;
CREATE CONSTRAINT gene_ensembl
FOR (n:Gene)
REQUIRE n.ensembl IS UNIQUE;
CREATE INDEX gene_symbol
FOR (a:Gene)
ON (a.symbol);

//Protein
CREATE CONSTRAINT protein_uniprot
FOR (n:Protein)
REQUIRE n.uniprot IS UNIQUE;
CREATE INDEX protein_symbol
FOR (a:Protein)
ON (a.symbol);

//Metabolite
CREATE CONSTRAINT metabolite_name
FOR (n:Metabolite)
REQUIRE n.name IS UNIQUE;
CREATE CONSTRAINT metabolite_refmet
FOR (n:Metabolite)
REQUIRE n.refmet IS UNIQUE;

//edge constraints
//CREATE CONSTRAINT constraint_name FOR ()-[r:REL_TYPE]-() REQUIRE r.property IS UNIQUE;