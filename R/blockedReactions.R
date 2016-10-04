#' @export blockedReactions
#' @author Andrés Pinzón <ampinzonv@unal.edu.co> - Mantainer: Daniel Camilo Osorio Hurtado <dcosorioh@unal.edu.co>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

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