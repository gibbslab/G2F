additionCost <- function(reaction,reference){
  mets_r <- metabolites(reference)
  mets <- lapply(reaction,metabolites)
  mets <- lapply(mets, function(metabolites){table(metabolites%in%mets_r)})
  cost <- lapply(mets, function(metabolites){metabolites["FALSE"]/sum(metabolites)})
  cost[is.na(cost)] <- 0
  cost <- as.vector(unlist(cost))
  return(cost)
}