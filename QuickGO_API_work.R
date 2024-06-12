#########################################
### Get annotations from QuickGO API ###
########################################

### Load libraries and data
## Libraries
library(rstudioapi)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggsankey)
library(english)
setwd(dirname(getActiveDocumentContext()$path))
## For the petitions to the API
library(httr)
library(jsonlite)
library(xml2)

## Functions
source("./functions/get_childs.R")
source("./functions/get_childs_recursive.R")
source("./functions/get_go_genes.R")

#### Get the GOs
### Get the childs of GO:0005576 - (extracellular region)
GOs_childs <- get_childs(go_id = c("GO:0005576")) 
### Get the childs recursively
GOs_childs_recursive <- get_childs_recursive(childs_data = GOs_childs)
## Remove redundant requests where parents are already childs of another term in the search to get the 
# best and cleaned search for genes !
GOs_childs_recursive <- GOs_childs_recursive %>%
  # Group by all cols except recursive pass
  group_by(across(c(-recursive_pass))) %>%
  # If we have more than a path that is repeated, get the one w the biggest recursive pass
  dplyr::slice(which.max(recursive_pass)) %>%
  ungroup()

# filter the childs
GOs_childs_recursive <- GOs_childs_recursive %>% 
  # filter only cellular component
  filter(aspect_parent == "cellular_component") %>% 
  filter(aspect == "cellular_component")
### NOW ALL RELATIONS ARE = part_of and is_a !
unique(GOs_childs_recursive$relation)

#### Get the proteins
### Proteins from GO:0005576 only (All but Reactome)
## Proteins with IEA
GO_proteins_for_extracell_NO_reactome_IEA <- get_go_genes(go_id = c("GO:0005576"), go_usage = "exact",
                                          proteome = "gcrpCan", geneProductSubset = "Swiss-Prot",
                                          assigned_by = c("Uniprot","InterPro","BHF-UCL","Ensembl","GO_Central","GOC","MGI","ARUK-UCL","ComplexPortal",
                                                          "ParkinsonsUK-UCL","HGNC-UCL"),
                                          taxon = "9606", taxon_usage = "exact", 
                                          evidences_to_remove = "",evidence_usage = "descendants",
                                          to_curate = F)
# filter for unique proteins
GO_proteins_for_extracell_NO_reactome_IEA <- GO_proteins_for_extracell_NO_reactome_IEA %>% 
  distinct()

## Proteins with IEA
GO_proteins_for_extracell_NO_reactome <- get_go_genes(go_id = c("GO:0005576"), go_usage = "exact",
                                                      proteome = "gcrpCan", geneProductSubset = "Swiss-Prot",
                                                      assigned_by = c("Uniprot","InterPro","BHF-UCL","Ensembl","GO_Central","GOC","MGI","ARUK-UCL","ComplexPortal",
                                                                      "ParkinsonsUK-UCL","HGNC-UCL"),
                                                      taxon = "9606", taxon_usage = "exact", 
                                                      evidence_usage = "descendants",
                                                      to_curate = F)
# filter for unique proteins
GO_proteins_for_extracell_NO_reactome <- GO_proteins_for_extracell_NO_reactome %>% 
  distinct()

