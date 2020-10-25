#' @export gapFind
#' @importFrom "minval" "orphanReactants" "orphanProducts"
#' @author Daniel Osorio <dcosorioh@unal.edu.co>
#' @title Find gaps in a metabolic network
#' @description This function identifies the gaps (metabolites not produced or not consumed in any other reaction or just involved in one reaction) for a set of stoichometric reactions
#' @param reactionList A set of stoichiometric reaction with the following format:
#'
#' \code{"H2O[c] + Urea-1-carboxylate[c] <=> 2 CO2[c] + 2 NH3[c]"}
#'
#' Where arrows and plus signs are surrounded by a "space character".
#' It is also expected that stoichiometry coefficients are surrounded by spaces, (nothe the "2" before the CO2[c] or the NH3[c]).
#' It also expects arrows to be in the form "\code{=>}" or "\code{<=>}".
#' Meaning that arrows like "\code{==>}", "\code{<==>}", "\code{-->}" or "\code{->}" will not be parsed and will lead to errors.
#' @param removeExchange A boolean value \code{FALSE} to ignore the exchange reactions to identify the gaps.
gapFind <- function(reactionList, removeExchange = FALSE) {
  if (class(reactionList) == "modelorg") {
    reactionList <- model2string(model = reactionList)
  }
  if (class(reactionList) == "character") {
    if (removeExchange == TRUE) {
      reactionList <- reactionList[!minval:::reactionType(reactionList = reactionList) == "Exchange reaction"]
    }
    gaps <- list()
    gaps$uptake <- orphanReactants(reactionList = reactionList)
    gaps$release <- orphanProducts(reactionList = reactionList)
    return(gaps)
  }
}
