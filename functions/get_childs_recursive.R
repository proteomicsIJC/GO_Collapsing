########################## 
## get_childs_recursive ##
##########################

## childs_data = output of get_childs

get_childs <- function(go_id){
  ##### CHILDREN
  base_url <- "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"
  ##go_id <- c("GO:0048308","GO:0048313") ## <- Prova 1 funciona ! Tots tenen childs
  ##go_id <- c("GO:0048308","GO:0048313","GO:1990730","GO:0036504","GO:1990484") ## <- Prova funciona !  2 no tenen childs
  ##go_id = c("GO:1990730") ## Prova <- funciona ! no té fills i no torna res
  child_list <- list()
  ## set a counter and a counter maximum to extract the data
  counter_limit <- length(go_id)
  counter <- 1
  for (i in 1:counter_limit){
    url <- paste0(base_url,go_id[i],"/children")
    response <- httr::GET(url = url)
    if (httr::http_type(response) == "application/json") {
      children <- httr::content(response, "parsed")
      children <- children[["results"]]
      names(children) <- children[[1]]$name
      child_list <- c(child_list, list(children))
      counter <- counter +1
      ###      if (i == 1){
      ###       children <- httr::content(response, "parsed")
      ###       children <- children[["results"]]
      ###       names(children) <- children[[1]]$name
      ###       child_list <- children}
      ###      else {
      ###       children <- httr::content(response,"parsed")
      ###       children <- children[["results"]]
      ###       names(children) <- children[[1]]$name
      ###       child_list <- list(child_list, children)}
    }}
  ## Clean children list
  ## Check that any element has at least one child
  ## Define a function to it
  search_children <- function(lst) {
    for (sub_lst in lst) {
      if (is.list(sub_lst)) {
        # If the element is a list, recursively search it
        if ("children" %in% names(sub_lst)) {
          return(TRUE)  # Found "children", return TRUE
        } else {
          # Recursively search the sublist
          found_children <- search_children(sub_lst)
          if (found_children) {
            return(TRUE)  # "children" found in the sublist, return TRUE
          }
        }
      }
    }
    return(FALSE)  # "children" not found in any sublist, return FALSE
  }  
  has_children <- search_children(child_list)
  if (has_children){
    final_children_list <- list()
    for (j in 1:length(child_list)){
      final_children_list[[j]] <- child_list[[j]][[1]]
      names(final_children_list)[j] <- final_children_list[[j]]$name}
  } else {
    stop("NO CHILDRENS !")
  }
  
  final_children_list <- final_children_list
  
  ## Trasnform the list into a dataset
  final_dataframe <- final_children_list %>%
    map_dfr(~ {
      parent_id <- .x$id
      parent_name <- .x$name
      children <- .x$children %>%
        map_dfr(~ as.data.frame(t(unlist(.x)), stringsAsFactors = FALSE))
      children$parent_id <- parent_id
      children$parent_name <- parent_name
      children
    })
  final_dataframe <- final_dataframe %>% 
    relocate(parent_id,parent_name, .before = id)
  
  ## Add the names of the IDs that are not in the list cuz no childs init
  # Outersect function
  outersect <- function(x, y) {
    sort(c(setdiff(x, y),
           setdiff(y, x)))}
  # The ountersect
  out <- outersect(x = go_id,
                   y = intersect(go_id,final_dataframe$parent_id))
  
  # define the function to filter the list of results for the ones with no childs
  filtered_list <- lapply(final_children_list, function(x) {
    if (length(x) <= 2) {
      return(x)
    }
  })
  # filter the list
  filtered_list <- filtered_list[!sapply(filtered_list, is.null)]
  
  ## Extract the info
  if (length(filtered_list)>0){
    out <- c()
    for (nc in 1:length(filtered_list)){
      n_ames <- filtered_list[[nc]]$name
      go_ids <- filtered_list[[nc]]$id
      out[nc] <- go_ids
      names(out)[nc] <- n_ames
    }
    
    to_add <- as.data.frame( 
      tibble::tibble(
        parent_id = out,
        parent_name = names(out),
        id = out,
        name = rep("NO child data for this element", times = length(out)),
        relation = rep("NO child data for this element", times = length(out)),
        hasChildren = rep("NO child data for this element")))
    final_dataframe <- rbind(final_dataframe, to_add)
  }
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
get_childs_recursive <- function(childs_data) {
  ## Filter for those terms that have childs
  original_df <- childs_data
  goids_1st <- childs_data %>% 
    filter(hasChildren == "TRUE")
  goids_1st <- goids_1st$id
  
  ## If the filtering gives no chids as a result, we will terminate the operation
  if (length(goids_1st) == 0) {
    print("NO CHILD TERMS TO ADD")
    final_data <- NULL
    ## If we have some result we will continue  
  } else {
    ## Define an initial list of dataframes that we will append recursively 
    recursive_list <- list()
    recursive_list[[1]] <- original_df
    original_rec <- recursive_list
    ## Define a condition that the last element will need to have
    ## to keep with the execution of the code
    condition <- function(df) { 
      filtered_df <- subset(df, hasChildren == "TRUE")
      return(nrow(filtered_df) > 0)
    }
    
    ## The while loop on the condition for the last element of the list
    while (condition(tail(recursive_list, n=1)[[1]])) {
      # Extract the last element of the list
      recursives_goids <- tail(recursive_list, n = 1)[[1]]
      # Filter for the terms with childs
      recursives_goids <- recursives_goids %>% 
        filter(hasChildren == "TRUE")
      recursives_goids <- recursives_goids$id
      print(recursives_goids)
      
      recursive_dataframe <- get_childs(go_id = recursives_goids)
      recursive_list <- c(recursive_list, list(recursive_dataframe))
      final_dataframe <- do.call(rbind,recursive_list) 
    }
    
    final_data <- final_dataframe
  }
  return(final_data)
}
