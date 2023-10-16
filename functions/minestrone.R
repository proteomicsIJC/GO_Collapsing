################ 
## minestrone ##
###############

## collapsed_list = list of the collapsing process
## original_data = BP@results table (the filtered version)

minestrone <- function(collapsed_list, original_data){
  ### make the soup
  soup <- as.data.frame(collapsed_list$parent_paths)
  soup$child <- rownames(soup)
  colnames(soup)[1] <- "parent"
  for (i in 1:length(rownames(soup))){
    if (is.na(soup$parent[i])){
      soup$parent[i] <- soup$child[i] 
    }
  }
  ### do the soup readable
  x <- soup$parent
  soup$parent <- original_data$Description[match(x, original_data$ID)]
  y <- soup$child
  soup$child <- original_data$Description[match(y, original_data$ID)]
  
  return(soup)
}
