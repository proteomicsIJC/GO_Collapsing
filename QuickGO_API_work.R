#########################################
### Get annotations from QuickGO API ###
########################################

### Load libraries and data
## Libraries
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
library(httr)
library(jsonlite)
library(xml2)

## Functions
source("./functions/get_go_genes.R")
source("./functions/get_child_terms.R")

## Data
# Data could also be simply a string vector of GO terms to extract info from 
GOs <- openxlsx::read.xlsx("./skeletal_muscle.xlsx")


skeletal_proteins <- get_go_genes(go_id = GOs$Term, 
                               go_usage = "descendants", proteome = "gcrpCan", assigned_by = c("Uniprot","InterPro","GO_Central","BHF-UCL"), 
                               geneProductSubset = "Swiss-Prot", 
                               taxon = "9606", 
                               taxon_usage = "exact")



