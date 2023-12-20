#####################################
### GO:BP annotation and collapse ###
#####################################

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

## Functions
source("../functions/general/fi")
source("./functions/collapseGO.R")
source("./functions/minestrone.R")
source("./functions/treemaping.R")

## Data
raw_data <- read.table("./raw_data/protein_clusters_k6.tsv", sep = "\t", dec = ".", header = T)

## Clean genes
g_acs_1 <- strsplit(raw_data$protein, "\\;")
good_names_1 <- first_accession(g_acs_1)
raw_data$Accession_1 <- good_names_1
raw_data <- raw_data %>%
  relocate(Accession_1, .after = protein)

## Clusters
# here, define a vector of UNIPROT IDs
cluster1 <- raw_data$Accession_1[raw_data$cluster == "group1"]

### GO:BP
## Universe
organism = "org.Hs.eg.db"

backgGenes <- keys(org.Hs.eg.db) # EntrezID  
backgGenes <- AnnotationDbi::select(org.Hs.eg.db, keys = backgGenes, columns = "UNIPROT") # Converts entrezID to GeneSymbol
backgGenes <- backgGenes$UNIPROT
backgGenes_kk <- AnnotationDbi::select(org.Hs.eg.db, keys = backgGenes, columns = "SYMBOL", keytype = "UNIPROT") # Converts entrezID to GeneSymbol
backgGenes_kk <- backgGenes_kk$SYMBOL

## Annotation of the genes
c1.go <- AnnotationDbi::select(org.Hs.eg.db, keys = cluster1, columns = "SYMBOL", keytype = "UNIPROT")

# Do the analysis
BP_go1 <- enrichGO(gene = c1.go$SYMBOL,
                   universe = backgGenes_kk,
                   OrgDb = "org.Hs.eg.db",
                   ont = "BP",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1, minGSSize = 5, keyType = "SYMBOL")


## Filter by count
BP_go1 <- gsfilter(BP_go1, by = "Count", min = 5)
+
## Write the data
openxlsx::write.xlsx(x = BP_go1@result, file = "./results/go_bp/go_bp1.xlsx")

## Plot results
ggplot(top_n(BP_go1@result, n = 5, wt = -p.adjust), aes(x=Count, y=reorder(Description,Count), fill=p.adjust))+
  geom_bar(stat = "identity", width = 0.5)+
  ggtitle(paste0("Cluster 1, GO:BP"))+
  ylab("")+xlab("")+theme_bw()+ labs(fill='FDR') +scale_fill_continuous(labels = scientific_format())+
  theme(axis.text = element_text(size = 15))


## Figure-wise
go1 <- ggplot(top_n(BP_go1@result, n = 5, wt = -p.adjust), aes(x=Count, y=reorder(Description,Count), fill=p.adjust))+
  geom_bar(stat = "identity", width = 0.5)+
  ggtitle(paste0("Cluster 1"))+
  ylab("")+xlab("")+theme_bw()+ labs(fill='FDR') +scale_fill_continuous(labels = scientific_format())+
  theme(legend.text = element_text(size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", size = 1, fill = NA),
        axis.text = element_text(size = 23, face = "bold"),
        legend.key.height = unit(1.8,"cm"),
        legend.key.width = unit(1.8,"cm"),
        legend.title = element_blank())

pdf("./plots/go_bp/go1.pdf", width = 16, height = 4.5, bg = NULL )
go1
dev.off()

### GO:BP collapse
## Get the info
# Tables
gobp1 <- BP_go1@result

# filter the data
gobp1 <- gobp1 %>% 
  filter(p.adjust < 0.05)

# Get the genes
genes1 <- c1.go$SYMBOL

# Get the paths
pathways1 <- BP_go1@geneSets

## Collapsing
collapsed1 <- collapseGO(functional_annot = gobp1, pathways = pathways1,genes = genes1,  mingsize = 5, ontology_to_look = "BP")

## Work the collapsed data
main_functions1 <- gobp1[gobp1$ID %in% collapsed1$mainPaths,]

## Save the collapsed data
openxlsx::write.xlsx(x = main_functions1, file = "./results/go_bp_collasped/main_paths1.xlsx")

## Do the plots again with the collapsed data
ggplot(top_n(main_functions1, n = 5, wt = -p.adjust), aes(x=Count, y=reorder(Description,Count), fill=p.adjust))+
  geom_bar(stat = "identity", width = 0.5)+
  ggtitle(paste0("Cluster 1, GO:BP"))+
  ylab("")+xlab("")+theme_bw()+ labs(fill='FDR') +scale_fill_continuous(labels = scientific_format())+
  theme(axis.text = element_text(size = 15))

## Figure-wise for the collapsed data
go1_collapsed <- ggplot(top_n(main_functions1, n = nrow(main_functions1), wt = -p.adjust), aes(x=Count, y=reorder(Description,Count), fill=p.adjust))+
  geom_bar(stat = "identity", width = 0.5)+
  ggtitle(paste0("Cluster 1"))+
  ylab("")+xlab("")+theme_bw()+ labs(fill='FDR') +scale_fill_continuous(labels = scientific_format())+
  theme(legend.text = element_text(size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", size = 1, fill = NA),
        axis.text = element_text(size = 15, face = "bold"),
        legend.key.height = unit(1.8,"cm"),
        legend.key.width = unit(1.8,"cm"),
        legend.title = element_blank())

pdf("./plots/go_bp_collapsed/go1_collapsed.pdf", width = 16, height = 4.5, bg = NULL )
go1_collapsed
dev.off()

### GO:BP collapsing process
# Collapse to data_frame
soup1 <- minestrone(collapsed_list = collapsed1, original_data = gobp1)

## Do the treemap
# Do the trick and add a size column to the soup dataframe
soup1$size <- 1

# You can also make size proportional to the Count of each functional term
x <- soup1$child
soup1$size <- gobp1$Count[match(x, gobp1$Description)]

# Plot
pdf("./plots/tree1.pdf", width = 10, height = 10, bg = NULL )
treemaping(soup = soup1, graph_title = "GO:BP collapse cluster 1")
dev.off()

### Save results table
## Save the treemap data
openxlsx::write.xlsx("./results/go_bp_collasped/tree1.xlsx", x = soup1)

## Save the gobp data
openxlsx::write.xlsx("./results/go_bp/go_bp1.xlsx", x = gobp1)













