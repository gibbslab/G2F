.as.searchable <- function(metabolite){
  metabolite <- gsub("\\(","\\\\(",metabolite)
  metabolite <- gsub("\\)","\\\\)",metabolite)
  metabolite <- gsub("\\[","\\\\[",metabolite)
  metabolite <- gsub("\\]","\\\\]",metabolite)
  metabolite <- gsub("\\-","\\\\-",metabolite)
  metabolite <- gsub("\\+","\\\\+",metabolite)
  metabolite <- gsub("\\_","\\\\_",metabolite)
  metabolite <- paste0("[[:space:]]",metabolite,"[[:space:]]")
  return(metabolite)
}


.reactant.fill <- function(metabolite,reference){
  reference <- paste0(" ",reference," ")
  reference <- as.vector(reference)[grep(.as.searchable(metabolite),reference)]
  unique(c(as.vector(reference)[grep(paste0("<?=>?[[:print:]]+?",.as.searchable(metabolite)),reference)],
           as.vector(reference)[grep(paste0("<=>[[:print:]]+?",.as.searchable(metabolite)),reference)],
           as.vector(reference)[grep(paste0(.as.searchable(metabolite),"[[:print:]]+?<=>"),reference)]))
}

.product.fill <- function(metabolite,reference){
  reference <- paste0(" ",reference," ")
  reference <- as.vector(reference)[grep(.as.searchable(metabolite),reference)]
  unique(c(as.vector(reference)[grep(paste0(.as.searchable(metabolite),"[[:print:]]+?<?=>?"),reference)],
           as.vector(reference)[grep(paste0("<=>[[:print:]]+?",.as.searchable(metabolite)),reference)],
           as.vector(reference)[grep(paste0(.as.searchable(metabolite),"[[:print:]]+?<=>"),reference)]))
}

.addition.cost <- function(reaction,reference){
  mets <- table(metabolites(reaction)%in%metabolites(reference))
  cost <- as.vector(mets["FALSE"]/sum(mets))
  return(as.numeric(cost[!is.na(cost)]))
}

.extract <- function(reaction){
  if (grepl(";",reaction)){
    parts <- unlist(strsplit(reaction,";[[:blank:]]+"))
    return(parts[length(parts)])
  }else{ 
    if (!grepl("rn:R[[:digit:]]+[[:blank:]]",reaction)){
      return(reaction)
    } else {
      return(unlist(strsplit(reaction,"rn:R[[:digit:]]+[[:blank:]]"))[2]) 
    }
  }
}
