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

## Data
# Data could also be simply a string vector of GO terms to extract info from 
GOs_childs <- get_childs(go_id = c("GO:0015727","GO:0010714","GO:0015129")) 
## GOs_childs <- get_childs(go_id = c("GO:0015727","GO:0010714")) 
### Remove redundant requests where parents are
### already cholds of another term in the search to get the best and cleaned search for genes !
GOs_childs_recursive <- get_childs_recursive(childs_data = GOs_childs) 
GOs_childs_recursive <- GOs_childs_recursive %>% 
  group_by(across(c(-recursive_pass))) %>% 
  slice(which.max(recursive_pass)) %>% 
  ungroup()

GO_proteins_for_GO_childs <- get_go_genes(go_id = c(unique(GOs_childs$parent_id,
                                                           GOs_childs$id)), 
                               go_usage = "exact", proteome = "gcrpCan", assigned_by = c("Uniprot","InterPro","GO_Central","BHF-UCL"), 
                               geneProductSubset = "Swiss-Prot", 
                               taxon = "9606", 
                               taxon_usage = "exact")


GO_proteins_for_GO_childs_recursive <- get_go_genes(go_id = c(unique(GOs_childs_recursive$parent_id,
                                                              GOs_childs_recursive$id)), 
                                          go_usage = "exact", proteome = "gcrpCan", assigned_by = c("Uniprot","InterPro","GO_Central","BHF-UCL"), 
                                          geneProductSubset = "Swiss-Prot", 
                                          taxon = "9606", 
                                          taxon_usage = "exact")

### Sankey plot for the passes
## https://stackoverflow.com/questions/64146971/sankey-alluvial-diagram-with-percentage-and-partial-fill-in-r 
## CODE FOR THE SANKEY PLOT
## Remove some data to work better, this might be eliminated
GOs_childs_recursive2 <- GOs_childs_recursive %>% 
  subset(select = c(parent_name,
                    name, recursive_pass))

## Get the data of the trajectories that we will get to explain the extraction of paths
GOs_childs_recursive_trajectories <- GOs_childs_recursive2 %>% 
  group_by(name, parent_name) %>% 
  slice(which.max(recursive_pass)) %>% 
  ungroup() %>% 
  filter(!name %in% parent_name) %>%
  arrange((recursive_pass))

## Initialize the trajectory list
trajectory_list <- list()
for (i in 1:length(rownames(GOs_childs_recursive_trajectories))){
  trajectory_list[[i]] <- as.data.frame(GOs_childs_recursive_trajectories[i,])
}





