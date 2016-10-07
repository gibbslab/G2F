#' @export gapFill
#' @importFrom "minval" "orphanReactants" "orphanProducts" "reactants" "products"
#' @author Kelly Botero <kjboteroo@unal.edu.co> - Maintainer: Daniel Osorio <dcosorioh@unal.edu.co>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title Find and fill gaps in a metabolic network
#' @description This function identifies the gaps and fills it from the stoichiometric reactions of a reference metabolic reconstruction using a weighting function.
#' @param reactionList 
#' @param reference 
#' @param limit 
#' @param woCompartment 
#' @param consensus
#' 
#' @examples 
#' # Downloading stoichiometric reactions
#' all <- getReference(organism = "all",sep = ";")
#' hsa <- getReference(organism = "hsa",sep = ";")
#' 
#' # gapFill
#' gapFill(reactionList = hsa$reaction,
#'         reference = all$reaction, 
#'         limit = 0.25,
#'         woCompartment = TRUE,
#'         consensus = FALSE)

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
