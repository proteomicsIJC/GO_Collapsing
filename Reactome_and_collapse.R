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

## Functions
source("./functions/first_accession.R")
source("./functions/collapseGO.R")
source("./functions/minestrone.R")
source("./functions/treemaping.R")
source("./functions/CollapseReactome.R")

## Data
raw_data <- read.table("./raw_data/protein_clusters_k6.tsv", sep = "\t", dec = ".", header = T)

## Clean genes
g_acs_1 <- strsplit(raw_data$protein, "\\;")
good_names_1 <- first_accession(g_acs_1)
raw_data$Accession_1 <- good_names_1
raw_data <- raw_data %>%
  relocate(Accession_1, .after = protein)

## Define a list of proteins
cluster1 <- raw_data$Accession_1[raw_data$cluster == "group1"]

### GO:BP
## Universe
organism = "org.Hs.eg.db"

backgGenes <- keys(org.Hs.eg.db) # EntrezID  
backgGenes <- AnnotationDbi::select(org.Hs.eg.db, keys = backgGenes, columns = "UNIPROT") # Converts entrezID to GeneSymbol
backgGenes <- backgGenes$UNIPROT
backgGenes_kk <- AnnotationDbi::select(org.Hs.eg.db, keys = backgGenes, columns = "ENTREZID", keytype = "UNIPROT") # Converts entrezID to GeneSymbol
backgGenes_kk <- backgGenes_kk$ENTREZID

## Annotation of the genes
c1.react <- AnnotationDbi::select(org.Hs.eg.db, keys = cluster1, columns = "ENTREZID", keytype = "UNIPROT")

# Do the analysis
reactome_1 <- signatureSearch::enrichReactome(gene = c1.react$ENTREZID,
                                              universe = backgGenes_kk,
                                              organism = "human",
                                              pvalueCutoff = 0.05,
                                              qvalueCutoff = 0.1, minGSSize = 5)



## Filter by count
reactome_1 <- reactome_1@result %>%
  filter(Count >= 5)

## Make gene name readable
reactome_1$Symbol <- "Not done"
reactome_1$UNIPROT <- "Not done"

## Entrez to symbol
reactome_1 <- EntreztoSym(data_table = reactome_1)

## Entrez to Uniprot
reactome_1 <- EntrezToUniprot(data_table = reactome_1)

## Write the data
## openxlsx::write.xlsx(x = reactome_1, file = "./results/reactome/reactome_c1.xlsx")

## Plot results
ggplot(top_n(reactome_1, n = 5, wt = -p.adjust), aes(x=Count, y=reorder(Description,Count), fill=p.adjust))+
  geom_bar(stat = "identity", width = 0.5)+
  ggtitle(paste0("Cluster 1, Reactome"))+
  ylab("")+xlab("")+theme_bw()+ labs(fill='FDR') +scale_fill_continuous(labels = scientific_format())+
  theme(axis.text = element_text(size = 15))

## Figure-wise
go1 <- ggplot(top_n(reactome_1, n = 5, wt = -p.adjust), aes(x=Count, y=reorder(Description,Count), fill=p.adjust))+
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

pdf("./plots/reactome/reactome_c1.pdf", width = 16, height = 4.5, bg = NULL )
go1
dev.off()

### Pathways as a list 
## get the paths and work with them
reactome_paths <- read.table("./raw_data/NCBI_Reactome_all.txt", sep = "\t", header = F)
colnames(reactome_paths) <- c("EntrezGene","Path","Link")
reactome_paths <- reactome_paths %>% 
  subset(select = c(EntrezGene,Path))

## subset human paths
reactome_paths <- reactome_paths[grepl("HSA", reactome_paths$Path),]

## collapse the dataset
reactome_paths_collapse <- reactome_paths %>% 
  dplyr::group_by(Path) %>%
  summarize_all(paste, collapse=",")

## do the list
path_list <- list()
for (i in 1:length(rownames(reactome_paths_collapse))){
  path_list[[i]] <- unlist(strsplit(reactome_paths_collapse$EntrezGene[i], "\\,"), use.names = F)
  names(path_list)[i] <- reactome_paths_collapse$Path[i]
}

### Reactome collapse
## Get the info
# Tables, they have already been created

# Get the genes
genes1 <- c1.react$ENTREZID

# Get the paths
pathways1 <- path_list[names(path_list) %in% reactome_1$ID]

## Collapsing
collapsed1 <- CollapseReactome(functional_annot = reactome_1, pathways = pathways1, genes = genes1, mingsize = 5, organism = "human")

## Work the collapsed data
main_functions1 <- reactome_1[reactome_1$ID %in% collapsed1$mainPaths,]

## Save the collapsed data
### openxlsx::write.xlsx(x = main_functions1, file = "./results/reactome_collapsed/main_paths1.xlsx")

## Do the plots again with the collapsed data
ggplot(top_n(main_functions1, n = 5, wt = -p.adjust), aes(x=Count, y=reorder(Description,Count), fill=p.adjust))+
  geom_bar(stat = "identity", width = 0.5)+
  ggtitle(paste0("Cluster 1, Reactome"))+
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

pdf("./plots/reactome_collapsed/reactome1_collapsed.pdf", width = 25, height = 4.5, bg = NULL )
go1_collapsed
dev.off()

### GO:BP collapsing process
# Collapse to data_frame
soup1 <- minestrone(collapsed_list = collapsed1, original_data = reactome_1)

## Do the treemap
# Do the trick and add a size column to the soup dataframe
soup1$size <- 1

# If wanted, put the size as the Count
x <- soup1$child
soup1$size <- gobp1$Count[match(x, gobp1$Description)]

# Plot
pdf("./plots/reactome_collapsed/tree1.pdf", width = 10, height = 10, bg = NULL )
treemaping(soup = soup1, graph_title = "Reactome: collapse cluster 1")
dev.off()

### Save results table
## Save the treemap data
### openxlsx::write.xlsx("./results/reactome_collapsed/tree1.xlsx", x = soup1)

