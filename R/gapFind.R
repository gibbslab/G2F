#' @export gapFind
#' @importFrom "minval" "orphanMetabolites"
#' @author Daniel Osorio <dcosorioh@unal.edu.co>
#' @title Find gaps in a metabolic network
#' @description This function identifies the gaps (metabolites not produced or not consumed in any other reaction or just involved in one reaction) for a set of stoichometric reactions
gapFind <- function(reactionList){
  orphanMetabolites(reactionList = reactionList)
}
