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
#' all <- all[!(all$reaction %in% eco$reaction),]
#'                     
#' # gapFill
#' gapFill(reactionList = eco$reaction,
#'         reference = all$reaction, 
#'         limit = 0.25,
#'         woCompartment = TRUE,
#'         consensus = FALSE)}

gapFill <- function(reactionList, reference, limit = 0.25, woCompartment=FALSE,consensus=FALSE){
  
  # Setting list as vector
  reactions <- as.vector(unique(reactionList))
  reference <- as.vector(unique(reference))
  
  # Check the syntax
  if(isTRUE(any(validateSyntax(reactions)=="FALSE") | any(validateSyntax(reference)=="FALSE"))){
    stop("Syntax error found in the stoichiometric reactions")
  }

  # Remove the compartments
  if(woCompartment==TRUE){
    reactions <- gsub("\\[[[:alnum:]]*(\\_)?[[:alnum:]]*\\]$","",as.vector(reactions))
  }
  
  # Extract all orphan metabolites from reactionList (OrphanOriginal)
  orphan <- orphanMetabolites(reactionList =  reactions)
  
  # Extract all reactants
  reactantsReference <- lapply(reference, function(reaction){reactants(reaction)})
  
  # Extract all products
  productsReference <- lapply(reference, function(reaction){products(reaction)})
  
  # do
  repeat{
    # Compute the addition cost for all stoichiometric reactions from the reference
    # Select stoichiometric reactions with additionCost lower or equal than limit
    belowLimit <- additionCost(reference = reactions,reaction = reference) <= limit
    # Extract all orphan reactants from reactionList
    orphanR <- orphanReactants(reactions)
    # Count the number of orphan reactants that are in orphanOriginal
    orphanR <- orphanR[orphanR %in% orphan]
    message(paste0(length(orphanR)," Orphan reactants found"))
    # Identify the reactions that contain orphan reactants in selected stoichiometric reactions
    includeOrphan <- unlist(lapply(reactantsReference, function(reaction){any(reaction%in%orphanR)}))
    # Select reactions to be added
    toAdd <- reference[includeOrphan & belowLimit]
    # If the number of orphanOriginals \in OrphanReactant is lower than OrphanOriginals \in orphanReactant \in orphans(reactionList \cup to.add)
    if(all(orphanR %in% orphanReactants(c(reactions,toAdd)))){
      break;
    } else {
      reactions <- unique(c(reactions,toAdd))
    }
  }
  repeat{
    # Compute the addition cost for all stoichiometric reactions from the reference
    # Select stoichiometric reactions with additionCost lower or equal than limit
    belowLimit <- additionCost(reference = reactions,reaction = reference) <= limit
    # Extract all orphan products from reactionList
    orphanP <- orphanProducts(reactions)
    # Count the number of orphan reactants that are in orphanOriginal
    orphanP <- orphanP[orphanP%in%orphan]
    message(paste0(length(orphanP)," Orphan products found"))
    # Identify the reactions that contain orphan reactants in selected stoichiometric reactions
    includeOrphan <- unlist(lapply(productsReference, function(reaction){any(reaction%in%orphanP)}))
    # Select reactions to be added
    toAdd <- reference[includeOrphan & belowLimit]
    if(all(orphanP %in% orphanReactants(c(reactions,toAdd)))){
      break;
    } else {
      reactions <- unique(c(reactions,toAdd))
    }
  }
  reactions <- unique(reactions)
  if(consensus == TRUE){
    return(reactions)
  } else{
    return(reactions[!reactions%in%reactionList])
  }
}

