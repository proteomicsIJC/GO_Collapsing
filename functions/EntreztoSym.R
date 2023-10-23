################ 
## EntreztoSym ##
###############

## data_table = A ErnichR@result data table type data with itemID as column of accession in Entrez

EntreztoSym <- function(data_table){
  for (i in 1:length(rownames(data_table))){
    entrezIDs <- data_table$itemID[i]
    entrezIDs <- unlist(strsplit(entrezIDs, "\\/"), use.names = F)
    geneSymbols <- AnnotationDbi::select(org.Hs.eg.db, keys = entrezIDs, columns = "SYMBOL")
    geneSymbols <- geneSymbols$SYMBOL
    geneSymbols <- paste0(geneSymbols, collapse = "/")
    data_table$Symbol[i] <- geneSymbols
  }
  return(data_table)
}
