#########################################
### Reactome annotation and collapse ###
########################################

### Load libraries and data
## Libraries
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
library(ggplot2)
library(ggfortify)
library(limma)
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(stringr)
library(janitor)
library(dplyr)
library(tidyr)
library(mice)
library(reshape2)
library(openxlsx)
library(NbClust)
library(factoextra)
library(readr)
library("org.Hs.eg.db")
library(clusterProfiler)
library(scales)
library(grid)
library(gridExtra)
library(cowplot)
library(httr)
library(jsonlite)
library(xml2)

## Functions
source("./functions/collapseGO.R")
source("./functions/minestrone.R")
source("./functions/treemaping.R")
source("./functions/CollapseReactome.R")
source("./functions/get_go_genes.R")

## Data
GOs <- openxlsx::read.xlsx("./skeletal_muscle.xlsx")

skeletal_proteins <- get_go_genes(go_id = GOs$Term, 
                               go_usage = "descendants", proteome = "gcrpCan", assigned_by = c("Uniprot","InterPro","GO_Central","BHF-UCL"), 
                               geneProductSubset = "Swiss-Prot", 
                               taxon = "9606", 
                               taxon_usage = "exact")



