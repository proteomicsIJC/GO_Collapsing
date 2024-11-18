######################## 
## get_interpro_pfam ##
#######################

library(httr)

## protein_id = Uniprot protein ID

get_interpro_pfam_data <- function(protein_id,taxon){
  
  base_url <- "https://www.ebi.ac.uk:443/interpro/api/entry/pfam/protein/UniProt/"
  base_url_extra <- "?extra_fields=entry_id"
  result_list <- list()
  for (i in 1:length(protein_id)){
    ### Protein ID part
    unique_protein <- protein_id[i]
    unique_protein <- paste0(unique_protein,"/")

    ### Taxon part
    unique_taxon <- taxon
    unique_taxon <- paste0("taxonomy/uniprot/",unique_taxon,"/")

    ## variable part
    variable_part <- paste0(unique_protein,unique_taxon)
    url <- paste(base_url,variable_part,base_url_extra,sep = "")
    
    ## Do the API petition
    response <<- httr::GET(url = url, accept("application/json"))
    stop_for_status(response)
    content <- httr::content(response, as = "text", encoding = "UTF-8")
    
    ## Deparse the content of the json
    if (content != "") {"
      "} else if (content == "") {
        next }
    
    ## Put all the json data to a single dataframe

    ## Append the dataframe into the list
  }
  # final_data_frame <- as.data.frame(do.call(rbind, final_list))
  # return(final_data_frame)
  return(content)
}


a <- get_interpro_pfam_data(protein_id = c("P68871"), taxon = "9606")


################################## RESPONSE FOR 
# Load required package
library(jsonlite)
library(tidyr)

# Your API response as a JSON string
api_response <- '{"count":1,"next":null,"previous":null,"results":[{"metadata":{"accession":"PF07968","name":"Leukocidin/Hemolysin toxin family","source_database":"pfam","type":"domain","integrated":"IPR016183","member_databases":null,"go_terms":null},"proteins":[{"accession":"p09616","protein_length":319,"source_database":"reviewed","organism":"1280","in_alphafold":true,"entry_protein_locations":[{"fragments":[{"start":63,"end":313,"dc-status":"CONTINUOUS"}],"representative":true,"model":"PF07968","score":7.3e-88}]}],"taxa":[{"accession":"1280","lineage":["1","131567","2","1783272","1239","91061","1385","90964","1279","1280"],"source_database":"uniprot"}],"extra_fields":{"entry_id":null}}]}'
api_response <- "{\"count\":2,\"next\":null,\"previous\":null,\"results\":[{\"metadata\":{\"accession\":\"PF00091\",\"name\":\"Tubulin/FtsZ family, GTPase domain\",\"source_database\":\"pfam\",\"type\":\"domain\",\"integrated\":\"IPR003008\",\"member_databases\":null,\"go_terms\":null},\"proteins\":[{\"accession\":\"p0a031\",\"protein_length\":390,\"source_database\":\"reviewed\",\"organism\":\"1280\",\"in_alphafold\":true,\"entry_protein_locations\":[{\"fragments\":[{\"start\":14,\"end\":174,\"dc-status\":\"CONTINUOUS\"}],\"representative\":false,\"model\":\"PF00091\",\"score\":2.3e-50}]}],\"taxa\":[{\"accession\":\"1280\",\"lineage\":[\"1\",\"131567\",\"2\",\"1783272\",\"1239\",\"91061\",\"1385\",\"90964\",\"1279\",\"1280\"],\"source_database\":\"uniprot\"}],\"extra_fields\":{\"entry_id\":null}},{\"metadata\":{\"accession\":\"PF12327\",\"name\":\"FtsZ family, C-terminal domain\",\"source_database\":\"pfam\",\"type\":\"domain\",\"integrated\":\"IPR024757\",\"member_databases\":null,\"go_terms\":null},\"proteins\":[{\"accession\":\"p0a031\",\"protein_length\":390,\"source_database\":\"reviewed\",\"organism\":\"1280\",\"in_alphafold\":true,\"entry_protein_locations\":[{\"fragments\":[{\"start\":222,\"end\":316,\"dc-status\":\"CONTINUOUS\"}],\"representative\":false,\"model\":\"PF12327\",\"score\":7.3e-31}]}],\"taxa\":[{\"accession\":\"1280\",\"lineage\":[\"1\",\"131567\",\"2\",\"1783272\",\"1239\",\"91061\",\"1385\",\"90964\",\"1279\",\"1280\"],\"source_database\":\"uniprot\"}],\"extra_fields\":{\"entry_id\":null}}]}"
api_response <- a

# Parse the JSON
parsed_data <- fromJSON(api_response)

# Extract the "results" section
results <- parsed_data$results

# Flatten "metadata"
metadata <- as.data.frame(results$metadata)

# Flatten "proteins" (nested structure)
proteins <- do.call(rbind, lapply(results$proteins, function(x) {
  protein_locations <- do.call(rbind, lapply(x$entry_protein_locations, function(y) {
    cbind(protein_accession = x$accession,
          protein_length = x$protein_length,
          start = y$fragments[[1]]$start,
          end = y$fragments[[1]]$end,
          score = y$score)
  }))
  protein_locations
}))

# Flatten "taxa"
taxa <- do.call(rbind, lapply(results$taxa, function(x) {
  cbind(accession = x$accession, lineage = paste(x$lineage, collapse = ";"))
}))

# View the extracted dataframes
print("Metadata:")
print(metadata)

print("Proteins:")
print(proteins)

print("Taxa:")
print(taxa)


