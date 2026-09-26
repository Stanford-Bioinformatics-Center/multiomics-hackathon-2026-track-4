# Stanford Multi-omics Hackathon 2026 Track 4

## MoTrPAC Knowledge Graph

*Can MoTrPAC data be made easier to connect, query, and explore?*

### Challenge

Transform a focused subset of MoTrPAC results and metadata into a documented knowledge graph, load it into [Neo4j](https://neo4j.com/), and customize the open-source [Knowledge Graph UI](https://pubmed.ncbi.nlm.nih.gov/40956607/) for interactive exploration.

### Data

Selected MoTrPAC molecules, tissues, omics layers, contrasts, pathways, and provenance, with approved identifier and ontology resources.

### Potential Outputs

- A Neo4j database and graph schema
- Reproducible ingestion workflow
- Searchable web interface

### Deployment

Deploy the Neo4j database container with Docker:

```docker-compose up -d neo4j```

You can then login to access the local Neo4j deplyment at <http://localhost:7474/>
