#' @export gapFill
#' @importFrom "minval" "orphanReactants" "orphanProducts" "orphanMetabolites" "reactants" "products"
#' @author Kelly Botero <kjboteroo@unal.edu.co> - Maintainer: Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title Find and fill gaps in a metabolic network
#' @description This function identifies the gaps and fills it from the stoichiometric reactions of a reference metabolic reconstruction using a weighting function.
#' @seealso \code{additionCost} function documentation.
#' @param reactionList A set of stoichiometric reaction with the following format:
#'
#' \code{"H2O[c] + Urea-1-carboxylate[c] <=> 2 CO2[c] + 2 NH3[c]"}
#'
#' Where arrows and plus signs are surrounded by a "space character".
#' It is also expected that stoichiometry coefficients are surrounded by spaces, (nothe the "2" before the CO2[c] or the NH3[c]).
#' It also expects arrows to be in the form "\code{=>}" or "\code{<=>}".
#' Meaning that arrows like "\code{==>}", "\code{<==>}", "\code{-->}" or "\code{->}" will not be parsed and will lead to errors.
#' @param reference A set of stoichiometric reaction with the same format of reactionList
#' @param limit An addition cost value to be used as a limit to select reactions to be added. Is calculated as NumberNewMetabolites/NumerOfMetabolites for each reaction.
#' @param woCompartment A boolean value \code{TRUE} to define if compartment labels should be removed of the reactionList stoichiometric reactions, \code{FALSE} is used as default.
#' @param consensus A boolean value \code{TRUE} to define if reactionList and newReactions should be reported as a unique vector or \code{FALSE} if just newReactions should be reported.
#' @param nRun The number of iterations to search for new gaps and fill them according to the score, the default is five.
#'
#' @examples
#' \dontrun{
#' # Downloading stoichiometric reactions
#' all <- getReference(organism = "all", sep = ";")
#' eco <- getReference(organism = "eco", sep = ";")
#'
#' # Filtering reactions
#' all <- mapReactions(
#'   reactionList = all$reaction %in% eco$reaction,
#'   referenceData = all,
#'   by = "bool",
#'   inverse = TRUE
#' )
#'
#' # gapFill
#' gapFill(
#'   reactionList = eco$reaction,
#'   reference = all$reaction,
#'   limit = 0.25,
#'   woCompartment = TRUE,
#'   consensus = FALSE
#' )
#' }
gapFill <- function(reactionList, reference, limit = 0.25, nRun = 5, woCompartment = FALSE, consensus = FALSE) {
  reactionList <- gsub("->", "_>", reactionList)
  reference <- gsub("->", "_>", reference)
  reference_reactants <- reactants(reference)
  reference_products <- products(reference)
  newR <- data.frame("addCost" = numeric(),"react" = character())
  n <- 0
  while (n < nRun) {
    oR <- orphanReactants(reactionList)
    oP <- orphanProducts(reactionList)

    filledReactions <- fillOrphanMetabolites(reference, reactionList, limit, reference_reactants, oR, newR)
    newR <- filledReactions$new
    reactionList <- filledReactions$summary
    filledReactions <- fillOrphanMetabolites(reference, reactionList, limit, reference_products, oP, newR)
    newR <- filledReactions$new
    reactionList <- filledReactions$summary

    n <- n + 1
    message(round(mean(!unique(c(oR, oP)) %in% orphanMetabolites(reactionList)), 2) * 100, "% gaps filled in the last iteration")
  }
  if (isTRUE(consensus)) {
    return(reactionList)
  } else {
    row.names(newR) <- NULL
    return(newR)
  }
}

# hsa <- getReference(organism = "hsa")
# reactionList <- sample(hsa$reaction,100)
# reference <- hsa$reaction[!hsa$reaction %in% reactionList]
# nR2 <- gapFill(reactionList = reactionList, reference = reference, nRun = 5)
# gF <- sapply(seq(0,1,0.025),function(sLimit){
#   x <- gapFill(reactionList = reactionList, reference = reference, limit = sLimit)
#   c(mean(!orphanMetabolites(reactionList) %in% orphanMetabolites(c(reactionList,x))),length(x))
# })
# gF <- t(gF)
# png("additionCost-Threshold.png", width = 3000, height = 1500, res = 300)
# par(mfrow=c(1,2))
# plot(x = seq(0,1,0.025), y = gF[,1], type = "l", xlim = c(0,1), ylim = c(0,1), xlab="Addition Cost", ylab = "Proportion of gaps filled")
# abline(v=0.25,col="red")
# plot(seq(0,1,0.025),gF[,2],type = "l", xlab ="Addition Cost", ylab="Number of reactions added")
# abline(v=0.25,col="red")
# dev.off()
