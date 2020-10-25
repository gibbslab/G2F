#' @export additionCost
#' @importFrom "minval" "metabolites"
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#' @title Calculate addition cost of a stoichiometric reaction in a metabolic network
#' @description This function calculates the addition-cost of a stoichiometric reaction in a metabolic network. The cost is computed for each reaction as the ratio of metabolites not included in the reference metabolic network over the total number of metabolites involved in the reaction.
#' @param reactionList A set of stoichiometric reaction with the following characteristics: \itemize{
#' \item Arrows symbols must be given in the form \code{'=>'} or \code{'<=>'}
#' \item Inverse arrow symbols \code{'<='} or other types as: \code{'-->'}, \code{'<==>'}, \code{'->'} will not be parsed and will lead to errors.
#' \item Arrow symbols and plus signs (\code{+}) must be surrounded by a space character
#' \item Stoichiometric coefficients must be surrounded by a space character and not by parentheses.
#' \item Each metabolite must have only one stoichiometric coefficient, substituents must be joined to metabolite name by a hyphen (\code{-}) symbol.
#' \item Exchange reactions have only one metabolite before arrow symbol
#' \item Compartments must be given between square brackets ([compartment]) joined at the end of metabolite name
#' }
#' Some examples of valid stoichiometric reactions are: \itemize{
#' \item \code{H2O[c] + Urea-1-Carboxylate[c] <=> 2 CO2[c] + 2 NH3[c]}
#' \item \code{ADP[c] + Phosphoenolpyruvate[c] => ATP[c] + Pyruvate[c]}
#' \item \code{CO2[c] <=> }
#' }
#' @param reference A different set of stoichiometric reactions with the same format of reactionList. This is the metabolic network expected to gap filling.
#' @examples
#' # Generate a reference, a vector of stoichiometric reactions as described above.
#'
#' example_reference <- vector()
#' example_reference[1] <- "A + B => C + D"
#' example_reference[2] <- "B + D => E"
#'
#' # Generate a set of reactions to compute the addition cost.
#'
#' example_reactionList <- vector()
#' example_reactionList[1] <- "A + E => F + G"
#' example_reactionList[2] <- "A + E <=> F"
#'
#' # The addition cost function will identify those metabolites not included in
#' # the reference (F and G) and will compute the addition cost for each reaction
#' # as the ratio of the number of metabolites not included in the reference
#' # over the total of metabolites in the reaction (2/4 and 1/3 respectively).
#'
#' additionCost(reactionList = example_reactionList, reference = example_reference)
#' # [1] 0.5000000 0.3333333
additionCost <- function(reactionList, reference) {
  if (class(reactionList) == "modelorg") {
    reactionList <- model2string(model = reactionList)
  }
  if (class(reference) == "modelorg") {
    reference <- model2string(model = reference)
  }
  mets_r <- metabolites(reference)
  mets <- lapply(reactionList, metabolites)
  cost <- unlist(lapply(mets, function(metabolites) {
    mean(!metabolites %in% mets_r)
  }))
  return(cost)
}
