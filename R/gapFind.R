#' @export gapFind
#' @importFrom "minval" "orphanMetabolites"
#' @author Daniel Osorio <dcosorioh@unal.edu.co>
#' @title Find gaps in a metabolic network
#' @description This function identifies the gaps (metabolites not produced or not consumed in any other reaction or just involved in one reaction) for a set of stoichometric reactions
gapFind <- function(reactionList, removeExchange = FALSE){
  if(class(reactionList) == "modelorg"){
    stoichiometricMatrix <- reactionList@S
    rownames(stoichiometricMatrix) <- reactionList@met_id
    reactionList <- minval:::rearmReactions(S = stoichiometricMatrix, reversible = reactionList@react_rev)
  }
  if (class(reactionList) == "character"){
    if(removeExchange == TRUE){
      reactionList <- reactionList[!minval:::reactionType(reactionList = reactionList) == "Exchange reaction"]
    }
    gaps <- list()
    gaps$uptake <- orphanReactants(reactionList = reactionList)
    gaps$release <- orphanProducts(reactionList = reactionList)
    return(gaps)
  }
}
