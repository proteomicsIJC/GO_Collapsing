#####################################
### GO:BP annotation and collapse ###
#####################################

#### Load libraries and data
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
source("./functions/collapseGO.R")
source("./functions/minestrone.R")
source("./functions/treemaping.R")

## Data
raw_data <- openxlsx::read.xlsx("./raw_data/example_data.xlsx", sheet = 1)

#### Data manipulation
### Create new columns in case it is required
# The names of the columns we mutate/merge or do not merge may change !
# the final idea of this manipulation is just to generate a vector of Gene.names (SYMBOLS)
# where each entry is realy a Gene.name and not a Gene.name followed by another like for
# example dmrt1;cyp19a1
raw_data <- raw_data %>%
  group_by(Accession) %>% 
  mutate(Gene.names_1 = unlist(strsplit(Gene.names, split = "\\;"), use.names = F)[1]) %>% 
  ungroup() %>% 
  relocate(Gene.names_1, .after = Gene.names) %>% 
  filter(!duplicated(Gene.names_1)) %>% 
  filter(!is.na(Gene.names_1))

#### Perform the enrichment analysis
### Initialize the analysis 
## Get the Background
# This is an example where the whole proteome is used as background, this may not
# be the best parctice but it gives us a base for the analysis
backgGenes <- keys(org.Hs.eg.db) # EntrezID  
backgGenes <- AnnotationDbi::select(org.Hs.eg.db, keys = backgGenes, columns = "UNIPROT") # Converts entrezID to GeneSymbol
backgGenes <- backgGenes$UNIPROT
backgGenes_kk <- AnnotationDbi::select(org.Hs.eg.db, keys = backgGenes, columns = "SYMBOL", keytype = "UNIPROT") # Converts entrezID to GeneSymbol
backgGenes_kk <- backgGenes_kk$SYMBOL

### Do the analysis
BP_go <- enrichGO(gene = raw_data$Gene.names_1, universe = backgGenes_kk,
                  OrgDb = "org.Hs.eg.db", ont = "BP",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.1, minGSSize = 5,
                  keyType = "SYMBOL")

## Filter by count
BP_go <- gsfilter(BP_go, by = "Count", min = 5)

## Write the data
openxlsx::write.xlsx(x = BP_go@result, file = "./results/GO_BP_non_collapsed.xlsx")

## Plot results
go_bp_plot <- ggplot(top_n(BP_go@result, n = 5, wt = -p.adjust), 
              aes(x=Count, y=reorder(Description,Count), fill=p.adjust))+
  geom_bar(stat = "identity", width = 0.5)+
  ggtitle(paste0("GO:Biological Process"))+
  ylab("")+
  xlab("")+
  theme_bw()+
  scale_fill_continuous(labels = scientific_format())+
  labs(fill='FDR') + 
  theme(legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        axis.ticks = element_line(colour = "black", linewidth =  1),
        panel.border = element_rect(colour = "black", linewidth = 1, fill = NA),
        axis.text = element_text(size = 10, face = "bold"),
        legend.key.height = unit(1.8,"cm"),
        legend.key.width = unit(0.6,"cm"))

pdf("./plots/GO_BP_non_collapsed.pdf", width = 16, height = 4.5, bg = NULL )
go_bp_plot
dev.off()

#### Gene Ontology collapse
ontology_table <- BP_go@result
ontology_table <- ontology_table %>% 
  filter(p.adjust < 0.05) %>%
  arrange(pvalue) 

collapsed_ontology <- collapseGO(functional_annot = BP_go@result,
                                 pathways = BP_go@geneSets,
                                 genes = raw_data$Gene.names_1,
                                 mingsize = 5, ## The same as in the filtering performed by gsfilter
                                 ontology_to_look = "BP")

## Work the collapsed data
main_functions <- BP_go@result[BP_go@result$ID %in% collapsed_ontology$mainPaths,]

## Save the collapsed data
openxlsx::write.xlsx(x = main_functions, file = "./results/GO_BP_collapsed.xlsx")

## Figure-wise for the collapsed data
go1_collapsed <- ggplot(top_n(main_functions, n = nrow(main_functions), wt = -p.adjust), 
                        aes(x=Count, y=reorder(Description,Count), fill=p.adjust))+
  geom_bar(stat = "identity", width = 0.5)+
  ggtitle(paste0("GO:Biological Process Collapsed"))+
  ylab("")+
  xlab("")+
  theme_bw()+
  scale_fill_continuous(labels = scientific_format())+
  labs(fill='FDR') + 
  theme(legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        axis.ticks = element_line(colour = "black", linewidth =  1),
        panel.border = element_rect(colour = "black", linewidth = 1, fill = NA),
        axis.text = element_text(size = 10, face = "bold"),
        legend.key.height = unit(1.8,"cm"),
        legend.key.width = unit(0.6,"cm"))

pdf("./plots/GO_BP_collapsed.pdf", width = 16, height = 4.5, bg = NULL )
go1_collapsed
dev.off()

### GO:BP collapsing process
# Collapse named list to data_frame
soup <- minestrone(collapsed_list = collapsed_ontology, original_data = ontology_table)

## Do a treemap for the collapsing process
# Do the trick and add a size column to the soup dataframe
soup$size <- 1

# You can also make size proportional to the Count of each functional term
soup <- soup %>% 
  group_by(child) %>% 
  mutate(size = ontology_table$Count[match(child,
                                           ontology_table$Description)])
## Plot
pdf("./plots/Treemap_GO_BP.pdf", width = 10, height = 10, bg = NULL )
treemaping(soup = soup, graph_title = "GO:BP collapse")
dev.off()

## Save the treemap data
openxlsx::write.xlsx("./results/GO_BP_collapsing_process.xlsx", x = soup)

## Final Save of the gobp data
ontology_table <- ontology_table %>% 
  mutate(collapsed_to = soup$parent[match(Description,
                                          soup$child)])
openxlsx::write.xlsx(x = ontology_table, file = "./results/GO_BP_all_data.xlsx")

######## The idea of having so much code that writes that much of excel files is that simply you create
######## save the excels you want but at first save all that things that you may need and then, if you do
######## do not need it or find it redudant (as it may be) remove it.









