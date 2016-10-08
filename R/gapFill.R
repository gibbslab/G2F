#' @export gapFill
#' @importFrom "minval" "orphanReactants" "orphanProducts" "reactants" "products"
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
#' 
#' @examples 
#' \dontrun{
#' # Downloading stoichiometric reactions
#' all <- getReference(organism = "all",sep = ";")
#' eco <- getReference(organism = "eco",sep = ";")
#' 
#' # Filtering reactions
#' all <- mapReactions(reactionList = all$reaction%in%eco$reaction,
#'                     referenceData = all,
#'                     by = "bool",
#'                     inverse = TRUE)
#'                     
#' # gapFill
#' gapFill(reactionList = eco$reaction,
#'         reference = all$reaction, 
#'         limit = 0.25,
#'         woCompartment = TRUE,
#'         consensus = FALSE)}

gapFill <- function(reactionList, reference, limit = 0.25, woCompartment=FALSE,consensus=FALSE){
  if(woCompartment==TRUE){
    reactionList <- gsub("\\[[[:alnum:]]*(\\_)?[[:alnum:]]*\\]$","",as.vector(reactionList))
  }
  reactions <- as.vector(unique(reactionList))
  reference <- as.vector(unique(reference))
  orphan <- unique(c(orphanReactants(reactions),orphanProducts(reactions)))
  repeat{
    orphan_r <- orphanReactants(reactions)
    orphan_r <- orphan_r[orphan_r%in%orphan]
    message(paste0(length(orphan_r)," Orphan reactants found"))
    to.add <- unique(unlist(lapply(orphan_r,function(orphan){reference[grep(orphan,reactants(reference),fixed = TRUE)]})))
    rxn <- to.add[additionCost(to.add,reactionList)<=limit]
    if(length(orphan_r)<=length(orphanReactants(c(reactions,rxn)))){
      break;
    } else {
      reactions <- unique(c(reactions,rxn))
    }
  }
  repeat{
    orphan_p <- orphanProducts(reactions)
    orphan_p <- orphan_p[orphan_p%in%orphan]
    message(paste0(length(orphan_p)," Orphan products found"))
    to.add <- unlist(lapply(orphan_p,function(orphan){reference[grep(orphan,products(reference),fixed = TRUE)]}))
    rxn <- to.add[additionCost(to.add,reactionList)<=limit]
    if(length(orphan_p)<=length(orphanProducts(c(reactions,rxn)))){
      break;
    } else{
      reactions <- unique(c(reactions,rxn))
    }
  }
  reactions <- unique(reactions)
  if(consensus == TRUE){
    return(reactions)
  } else{
    return(reactions[!reactions%in%reactionList])
  }
}
