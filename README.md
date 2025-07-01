# Plasmodium falciparum egress disrupts endothelial junctions and activates JAK-STAT signaling in a microvascular 3D blood-brain barrier model

Analysis described in https://www.biorxiv.org/content/10.1101/2024.10.15.618439v2.full


## Scripts
Main analysis scripts

- iRBC-egress analysis: MB10x01_analysis_final
- iRBC analysis: MB10x02_analysis_final, MB10x02_analysis_onlyBBB_final
- additional analyses: additional_analyses_final

The raw data are available in the ArrayExpress database under accession code E-MTAB-14463.
The pre-processed singlecellexperiment objects can be found in [Pre-processed Objects](data/pre-processed_objects/) and can be used to re-run downstream analysis and re-produce figures. Pre-processed objects are loaded in *Section 2: Downstream analysis* of every script.

## Data

- pre-processed objects: pre-processed objects to be used to re-run downstream analysis and re-produce figures
- demultiplex_results: results of demultiplexing MULTI-seq barcodes using the deMULTIplex package
- MAST_DE_results: results of MAST differential expression analysis run on the cluster ([cluster](scripts/cluster/))