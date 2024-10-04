################## 
## get_go_terms ##
##################

## uniprot_id = character vector of uniprot_ids to get genes from.
## proteome = character vector of proteome reference to use: none, 
# gcrpCan (Gene Centric Reference Proteome Canonical),
# gcrpIso (Gene Centric Reference Proteome IsoForm).
## assigned_by = character vector of labs that have assigned the gene/protein to the path if NULL all labs will be into the analysis.
## geneProductSubset = character vector of the Uniprot version to use reviewed (Swiss-Prot) and/or unreviewed (TrEMBL).
## taxon = character vector of taxon ids to retrive information from.
## taxon_usage = character vector for the taxon usage, it can be exact or descendants
## evidence_to_remove = character vector of those evidences which must NOT be retrieved from the API.
## evidence usage = exact or descendant
## to_curate = T or F. If T look only for "ECO:0000245","ECO:0000307","ECO:0000501" and you'll be able to remove the "bad" genes from annotation

get_go_terms <- function(uniprot_id, proteome, 
                         assigned_by = NULL,
                         geneProductSubset = "Swiss-Prot",
                         taxon, taxon_usage = "exact",
                         evidences_to_remove =  c("ECO:0000245","ECO:0000307","ECO:0000501"), evidence_usage = "exact",
                         to_curate = F){
  ## GO_list
  final_list <- list()
  
  ### Asigned by
  centers <- c("ARUK-UCL","ASAP","AgBase","Alzheimers_University_of_Toronto","AspGD",
               "BHF-UCL","CACAO","CAFA","CGD","CollecTF","ComplexPortal","DisProt","EcoCyc","EcoliWiki","Ensembl","EnsemblFungi","EnsemblMetazoa","EnsemblPlants",
               "EnsemblProtists","GOC","GO_Central","GR","GeneDB","HGNC","HGNC-UCL","HPA","IntAct","InterPro","JaponicusDB","LIFEdb","MGI","NTNU_SB","PHI-base",
               "ParkinsonsUK-UCL","PomBase","PseudoCAP","RGD","RHEA","RNAcentral","Reactome","SGD","SYSCILIA_CCNET","SynGO",
               "SynGO-UCL","TAIR","UniProt","WB","YuBioLab","ZFIN","dictyBase","iBB")
  
  if (is.null(assigned_by)){
    assigned_by <- centers
  }
  
  ### Evidences
  evidences <- c("Inferred from Experiment (EXP)	 ECO:0000269",
                 "Inferred from Direct Assay (IDA)	 ECO:0000314",
                 "Inferred from Physical Interaction (IPI)	 ECO:0000353",
                 "Inferred from Mutant Phenotype (IMP)	 ECO:0000315",
                 "Inferred from Genetic Interaction (IGI)	 ECO:0000316",
                 "Inferred from Expression Pattern (IEP)	 ECO:0000270",
                 "Inferred from High Throughput Experiment (HTP)	 ECO:0006056",
                 "Inferred from High Throughput Direct Assay (HDA)	 ECO:0007005",
                 "Inferred from High Throughput Mutant Phenotype (HMP)	 ECO:0007001",
                 "Inferred from High Throughput Genetic Interaction (HGI)	 ECO:0007003",
                 "Inferred from High Throughput Expression Pattern (HEP)	 ECO:0007007",
                 "Inferred from Sequence or structural Similarity (ISS)	 ECO:0000250",
                 "Inferred from Sequence Orthology (ISO)	 ECO:0000266",
                 "Inferred from Sequence Alignment (ISA)	 ECO:0000247",
                 "Inferred from Sequence Model (ISM)	 ECO:0000255",
                 "Inferred from Genomic Context (IGC)	 ECO:0000317",
                 "Inferred from Biological aspect of Ancestor (IBA)	 ECO:0000318",
                 "Inferred from Biological aspect of Descendant (IBD)	 ECO:0000319",
                 "Inferred from Key Residues (IKR)	 ECO:0000320",
                 "Inferred from Rapid Divergence (IRD)	 ECO:0000321",
                 "Inferred from Reviewed Computational Analysis (RCA)	 ECO:0000245",
                 "Traceable Author Statement (TAS)	 ECO:0000304",
                 "Non-traceable Author Statement (NAS)	 ECO:0000303",
                 "Inferred by Curator (IC)	 ECO:0000305",
                 "No biological Data available (ND)	 ECO:0000307",
                 "Inferred from Electronic Annotation (IEA) 	ECO:0000501")
  
  evidences <- gsub(pattern = " ", replacement = "", x = evidences)
  evidences_data <- data.frame(data.frame(do.call(rbind, strsplit(evidences, "\t")), stringsAsFactors = FALSE))
  names(evidences_data) <- c("Description", "Code")
  
  evidences_data <- evidences_data |> 
    dplyr::filter(!Code %in% evidences_to_remove, .preserve = F)
  
  evidences <- evidences_data$Code
  
  if (to_curate){
    evidences <- c("ECO:0000245","ECO:0000307","ECO:0000501")
  }
  
  ##### GET THE TERMS !
  base_url <- "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?"
  result_list <- list()
  for (i in 1:length(uniprot_id)){
    ### Proteome
    unique_proteomes <- paste(proteome, collapse = "%2C")
    unique_proteomes <- paste0("proteome=",unique_proteomes)
    ###Include fields
    unique_include_fields <- "includeFields=goName"
    ### Assigned by
    unique_assigned_by <- paste(assigned_by, collapse = "%2C")
    unique_assigned_by <- paste0("assignedBy=",unique_assigned_by)
    ### Select fields
    unique_select_fields <-"selectedFields=geneProductId&selectedFields=symbol&selectedFields=qualifier&selectedFields=goId&selectedFields=goAspect&selectedFields=goName&selectedFields=evidenceCode&selectedFields=goEvidence&selectedFields=taxonId&selectedFields=assignedBy&selectedFields=synonyms&selectedFields=name&selectedFields=type"
    ### Gene product
    unique_uniprotid <- uniprot_id[i]
    unique_uniprotid <- paste0("geneProductId=",unique_uniprotid)
    ### Gene product type
    unique_geneProductType <- "geneProductType=protein"
    ### Gene subset
    unique_geneProductSubset <- paste(geneProductSubset, collapse = "%2C")
    unique_geneProductSubset <- paste0("geneProductSubset=",unique_geneProductSubset)
    ### Taxon
    unique_taxon <- paste(taxon, collapse = "%2C")
    unique_taxon <- paste0("taxonId=",unique_taxon)
    
    ### Taxon Usage
    unique_taxon_usage <- taxon_usage
    unique_taxon_usage <- paste0("taxonUsagee=",unique_taxon_usage)
    
    ### Evidences
    evidences <- gsub(pattern = "\\:", replacement = "%3A", x = evidences)
    unique_evidence <- paste(evidences, collapse = "%2C")
    unique_evidence <- paste0("evidenceCode=",unique_evidence)
    
    ## Evidences usage
    unique_evidence_usage <- evidence_usage
    unique_evidence_usage <- paste0("evidenceCodeUsage=",unique_evidence_usage)
    
    
    ### Generate the URL
    variable_part <- paste(unique_proteomes,unique_include_fields, unique_assigned_by, unique_select_fields, unique_uniprotid,
                           unique_geneProductType, unique_geneProductSubset,unique_evidence,unique_evidence_usage, unique_taxon, unique_taxon_usage, sep = "&")
    url <- paste(base_url,variable_part,sep = "")
    
    ### Do the API petition
    response <- httr::GET(url = url, accept("text/tsv"))
    stop_for_status(response)
    content <- httr::content(response, as = "text")
    
    ## Transform the character vector into a dataframe
    lines <- strsplit(content,"\n")[[1]]
    data <- lapply(lines, function(line) unlist(strsplit(line, "\t")))
    go_id_data <- as.data.frame(do.call(rbind, data))
    colnames(go_id_data) <- go_id_data[1, ]
    go_id_data <- go_id_data[-1, ]
    
    if (nrow(go_id_data) >= 10000){
      cat(paste0("Hey, this GO Term info might be truncated ",go_id[i],"\r"))
      cat("\n")
    }
    ## Append the list
    final_list[[i]] <- go_id_data
    cat(paste0(i,"\r"))
  }
  final_data_frame <- as.data.frame(do.call(rbind, final_list))
  return(final_data_frame)
}

# QuickGO_call <- get_go_genes(uniprot_id = "A0A075B6J1", proteome = "gcrpCan",
#                              geneProductSubset = "Swiss-Prot", assigned_by = "Ensembl",
#                              taxon = "9606", taxon_usage = "exact", evidence_usage = "exact")
 
