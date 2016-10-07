#' @export additionCost
#' @importFrom "minval" "metabolites"
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#' @title Calculate the cost of addition of a stoichiometric reaction in a metabolic network
#' @description For a given set of stoichiometric reactions this function calculates the cost of addition in a reference metabolic network. 
#' The cost is calculated by dividing the amount of non included metabolites in the reference metabolic network over the total number of metabolites involved in the reaction.
#' @param reaction A stoichiometric reaction with the following format: 
#' 
#' \code{"H2O[c] + Urea-1-carboxylate[c] <=> 2 CO2[c] + 2 NH3[c]"} 
#' 
#' Where arrows and plus signs are surrounded by a "space character".
#' It is also expected that stoichiometry coefficients are surrounded by spaces, (nothe the "2" before the CO2[c] or the NH3[c]).
#' It also expects arrows to be in the form "\code{=>}" or "\code{<=>}". 
#' Meaning that arrows like "\code{==>}", "\code{<==>}", "\code{-->}" or "\code{->}" will not be parsed and will lead to errors.
#' @param reference A set of stoichiometric reaction with the same format of reaction.
#' @examples 
#' \dontrun{
#' # Downloading stoichiometric reactions of reference
#' hsa <- getReference(organism = "hsa")
#'
#' # Calculating cost
#' additionCost(reaction = "alpha-Amino acid + H2O + NAD+ <=> 2-Oxo acid + Ammonia + NADH + H+",
#'              reference = hsa$reaction)
#' }
additionCost <- function(reaction,reference){
  mets_r <- metabolites(reference)
  mets <- lapply(reaction,metabolites)
  mets <- lapply(mets, function(metabolites){table(metabolites%in%mets_r)})
  cost <- lapply(mets, function(metabolites){metabolites["FALSE"]/sum(metabolites)})
  cost[is.na(cost)] <- 0
  cost <- as.vector(unlist(cost))
  return(cost)
}