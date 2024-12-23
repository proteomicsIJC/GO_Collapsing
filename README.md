# Functional Annotation and Functional Annotation Collapse

This script is designed to summarize Gene Ontology (GO) enrichment results. After performing a Functional Annotation Enrichment Analysis or Over-Representation Analysis (ORA), the results can be collapsed to generate a summarized list of over-represented pathways in a given set of genes.

Once the primary ORA results are obtained using `clusterProfiler`, users can apply the `CollapseGO` function to streamline the pathways. For each pathway, `CollapseGO` checks whether, using the genes of a pathway (and the specified user-defined background), any other pathway achieves a nominal p-value lower than the reference value set by the user at the ORA level. If so, the pathways with a lower p-value will be collapsed into the main pathway.

The function outputs a two-element list:
- **mainPaths**: A list of pathways that cannot be collapsed into others or have other pathways collapsed into them.
- **parent_paths**: The original list of pathways with an indication of the pathways they have been collapsed into.

Additionally, this repository includes the `minestrone` and `treemapping` functions, which can be used to generate a treemap plot for visualization.
