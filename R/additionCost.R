#' @export additionCost
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#' @title Calculate the cost of addition of a stoichiometric reaction in a metabolic network
#' @description For a given set of stoichiometric reactions this function calculates the cost of addition in a reference metabolic network. 
#' The cost is calculated by dividing the amount of non included metabolites in the reference metabolic network over the total number of metabolites involved in the reaction.
additionCost <- function(reaction,reference){
  mets_r <- metabolites(reference)
  mets <- lapply(reaction,metabolites)
  mets <- lapply(mets, function(metabolites){table(metabolites%in%mets_r)})
  cost <- lapply(mets, function(metabolites){metabolites["FALSE"]/sum(metabolites)})
  cost[is.na(cost)] <- 0
  cost <- as.vector(unlist(cost))
  return(cost)
}