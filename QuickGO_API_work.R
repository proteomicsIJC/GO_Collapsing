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
GOs_childs <- get_childs(go_id = "GO:0015727") 
GOs_childs_recursive <- get_childs_recursive(childs_data = GOs_childs)

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
last_pass <- max(GOs_childs_recursive$recursive_pass)
GOs_childs_recursive <-  GOs_childs_recursive %>% 
  mutate(next_pass = ifelse(last_pass == recursive_pass,
                            NA,recursive_pass))

## Passes work
# Transform it to words
GOs_childs_recursive$recursive_pass <- as.character(ordinal(GOs_childs_recursive$recursive_pass))
GOs_childs_recursive$next_pass <- as.character(ordinal(GOs_childs_recursive$next_pass))

# Passes as factor
GOs_childs_recursive$recursive_pass <- as.factor(GOs_childs_recursive$recursive_pass)
GOs_childs_recursive$next_pass <- as.factor(GOs_childs_recursive$next_pass)

# Reorder of the levels
GOs_childs_recursive$recursive_pass <- factor(GOs_childs_recursive$recursive_pass, levels = c("first","second","third"))
GOs_childs_recursive$next_pass <- factor(GOs_childs_recursive$next_pass, levels = c("first","second","third"))

## Plot the thing !
pathways_connections <- ggplot(data = GOs_childs_recursive,
                               aes(x = recursive_pass, 
                                   next_x = next_pass, 
                                   node = parent_id, 
                                   next_node = parent_id,
                                   fill = factor(parent_id)))+ 
  geom_sankey(flow_alpha = 1,
              node.color = "black",           
              show.legend = TRUE)
pathways_connections




