#' @export blockedReactions
#' @author Andrés Pinzón <ampinzonv@unal.edu.co> - Mantainer: Daniel Camilo Osorio Hurtado <dcosorioh@unal.edu.co>
#' @title Identify blocked reactions in a metabolic network
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @description A blocked reaction in a metabolic network is a reaction that not participate in any optimization solution. This function set as objective function each one of the reactions (one by time) in the model, and identifies the reactions without flux under all scenarios.
#' @param model A valid model for the \code{'sybil'} package. An object of class modelorg.
#' @return A vector with the reaction ids of the blocked reactions
#' @examples
#' # Loading a model for the 'sybil' package
#' data("Ec_core")
#' 
#' # Identifying blocked reactions
#' blockedReactions(Ec_core)
#' @keywords Blocked reactions genome scale metabolic reconstruction
blockedReactions <- function(model){
  if(!is.loaded("sybil")){require("sybil")}
  locked <- NULL
  pb <- txtProgressBar(min = 1,max = model@react_num,style=3)
  for (reaction in 1:model@react_num) {
    setTxtProgressBar(pb, reaction)
    model@obj_coef <- rep(0, model@react_num)
    model@obj_coef[reaction] <- 1
    FBA <- optimizeProb(model)
    locked <- unique(c(locked, model@react_id[as.vector(FBA@fluxdist@fluxes!=0)]))
  }
  close(pb)
  locked <- model@react_id[!model@react_id%in%locked]
  return(locked)
}