################ 
## CollapseReactome ##
###############

## functional_annot = BP@results table
## pathways = BP@GeneSets
## genes = gene symbols (entrez or whatever) entered to extract the paths to collapse 
## mingsize = minimum genesize to take i count (putting the same as in the first versio shall be good)
## ontology to work with = BP,CC,MF

CollapseReactome <- function(functional_annot, pathways, genes, mingsize, organism = "human"){
  ## filter enrichment result
  functional_annot <- functional_annot %>% 
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::arrange(pvalue) 
  
  ## set the universe
  universe <- genes
  
  ## do some lapplies to the paths
  pathways <- pathways[functional_annot$ID]
  pathways <- lapply(pathways, intersect, universe)
  
  ## set a parent paths
  parentPaths <- setNames(rep(NA, length(pathways)), names(pathways))
  
  ## start looping
  for (i in 1:length(pathways)) {
    ## set the path to check if is parent
    p <- names(pathways)[i]
    if (!is.na(parentPaths[p])) {
      next
    }
    
    ## set the paths to check if are childs 
    paths_to_check <- setdiff(names(which(is.na(parentPaths))), p)
    
    ## Initialize the minPval vector
    minPval <- setNames(rep(1, length(paths_to_check)), paths_to_check)
    
    ## Our universe (u2)
    u2 <- pathways[[p]]

    gobp_u2 <- signatureSearch::enrichReactome(gene = u2, universe = genes,
                        organism = organism,
                        pvalueCutoff = 0.05, qvalueCutoff = 0.01, minGSSize = mingsize)
    
    gobp_u2 <- gobp_u2@result
    gobp_u2 <- gobp_u2 %>% 
      dplyr::filter(p.adjust < 0.05)
    
    # subset the pvalues
    gobp_u2_pval <- gobp_u2$pvalue
    names(gobp_u2_pval) <- gobp_u2$ID
    gobp_u2_pval <- gobp_u2_pval[names(gobp_u2_pval) %in% paths_to_check]
    
    ## minPval new
    # common names
    common_minPval <- intersect(names(gobp_u2_pval), names(minPval))
    # minimum of common names
    min_of_common <- pmin(gobp_u2_pval[common_minPval], minPval[common_minPval])
    # find uncommon names
    uncommon_names <- setdiff(names(minPval), names(gobp_u2_pval))
    # get the uncommon values from minPval
    uncommon_minPval <- minPval[uncommon_names]
    minPval <- c(uncommon_minPval,min_of_common)
    parentPaths[names(which(minPval < 0.05))] <- p
    
  }
  return(list(mainPaths = names(which(is.na(parentPaths))),
              parent_paths = parentPaths))
}
