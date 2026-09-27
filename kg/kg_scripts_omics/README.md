Overview of scripts listed in scripts_omics

1. convert.R: Converts the rat MOTrPAC R package (.rda) into TSV tables.

2. convert_h.R: Does the same for the human package.

3. build_sig_graph.py: Builds the knowledge graph (kg/exports_sig/): significance rules, weights, pathways, interactions, co-regulation, phenotypes, citations.

4. export_explorer_data.py: Writes explorer/data.json for the web explorer.

5. build_adjacency.py: Writes the adjacency matrices in kg/adjacency/

6. build_kg.py: Older version 0.1 build that keeps all measured molecules (build_sig_graph.py reuses its helper functions).

7. load_graph.py: Loads the graph into Python (NetworkX), with an example query.
