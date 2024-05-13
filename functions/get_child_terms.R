##################### 
## get_child_terms ##
#####################

## go_id = character vector of go_ids to get childs from.

get_childs <- function(go_id){
  ##### CHILDREN
  base_url <- "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"
  child_list <- list()
  for (i in 1:length(go_id)){
    url <- paste0(base_url,go_id[i],"/children")
    response <- httr::GET(url = url)
    if (httr::http_type(response) == "application/json") {
      if (i == 1){
        children <- httr::content(response, "parsed")
        children <- children[["results"]]
        names(children) <- children[[1]]$name
        child_list <- children}
      else {
        children <- httr::content(response,"parsed")
        children <- children[["results"]]
        names(children) <- children[[1]]$name
        child_list <- list(child_list, children)}
    }}
  ## Clean children list
  final_children_list <- list()
  for (j in 1:length(child_list)){
    final_children_list[[j]] <- child_list[[j]][[1]]
    names(final_children_list)[j] <- final_children_list[[j]]$name
  }
  final_children_list <- final_children_list
  
  ## Trasnform the list into a dataset
  final_dataframe <- final_children_list %>%
    purrr::map_dfr(~ {
      parent_id <- .x$id
      parent_name <- .x$name
      children <- .x$children %>%
        purrr::map_dfr(~ as.data.frame(t(unlist(.x)), stringsAsFactors = FALSE))
      children$parent_id <- parent_id
      children$parent_name <- parent_name
      children
    })
  final_dataframe <- final_dataframe %>% 
    relocate(parent_id,parent_name, .before = id)
  
  ##### ASPECT
  base_url <- "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"
  aspect_list <- list()
  go_id <- c(unique(final_dataframe$parent_id), 
             final_dataframe$id)
  for (k in 1:length(go_id)){
    url <- paste0(base_url,go_id[k])
    response <- httr::GET(url = url)
    if (httr::http_type(response) == "application/json") {
      aspect <- httr::content(response, "parsed")
      aspect <- aspect[["results"]]
      aspect_list[[k]] <- aspect}}
  ## ids_and_aspect
  ids_and_aspect <- c()
  for (l in 1:length(aspect_list)){
    ids_and_aspect[l] <- aspect_list[[l]][[1]]$aspect
    names(ids_and_aspect)[l] <- aspect_list[[l]][[1]]$id
  }
  ## transfrom into a dataframe
  ids_and_aspect <- as.data.frame(tibble::tibble(
    id = names(ids_and_aspect),
    aspect = ids_and_aspect))
  ## separate parents and childs information
  ids_and_aspect_parents <- ids_and_aspect %>% 
    filter(id %in% final_dataframe$parent_id) %>% 
    distinct()
  colnames(ids_and_aspect_parents)[2] <- "aspect_parent"
  ids_and_aspect_childs <- ids_and_aspect %>% 
    filter(id %in% final_dataframe$id) %>% 
    distinct()
  final_dataframe <- merge(final_dataframe, ids_and_aspect_parents, 
                           by.x = "parent_id", by.y = "id")
  final_dataframe <- merge(final_dataframe, ids_and_aspect_childs,
                           by.x = "id", by.y = "id")
  final_dataframe <- final_dataframe %>% 
    select(parent_id, parent_name, aspect_parent,
           id,name,relation,hasChildren,aspect)
  return(final_dataframe)
}
