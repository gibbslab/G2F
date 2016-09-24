#' @export gapFill
#' @author Kelly Botero <kjboteroo@unal.edu.co> - Maintainer: Daniel Osorio <dcosorioh@unal.edu.co>
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

gapFill <- function(reactionList, reference, limit = 0.5, woCompartment=FALSE,consensus=FALSE){
  if(woCompartment==TRUE){
    reactionList <- gsub("\\[[[:alpha:]]\\]","",reactionList)
  }
  reactions <- reactionList
  orphan <- orphanReactants(reactions)
  repeat{
    orphan_r <- orphanReactants(reactions)
    orphan_r <- orphan_r[orphan_r%in%orphan]
    message(paste0(length(orphan_r)," Orphan reactants found"))
    to.add <- unique(unlist(sapply(orphan_r, function(metabolite){.reactant.fill(metabolite,reference)})))
    rxn <- unique(c(reactions,to.add[unlist(sapply(to.add, function(reaction){.addition.cost(reaction,reactions)}))<limit]))
    if(length(orphanReactants(rxn))<length(orphan_r)){
      reactions <- unique(c(reactions,rxn))
    } else {
      break
    }
  }
  oreactants <- reactions
  reactions <- reactionList
  orphan <- orphanProducts(reactions)
  repeat{
    orphan_p <- orphanProducts(reactions)
    orphan_p <- orphan_p[orphan_p%in%orphan]
    message(paste0(length(orphan_p)," Orphan products found"))
    to.add <- unique(unlist(sapply(orphan_p, function(metabolite){.product.fill(metabolite,reference)})))
    rxn <- unique(c(reactions,to.add[unlist(sapply(to.add, function(reaction){.addition.cost(reaction,reactions)}))<limit]))
    if(length(orphanProducts(rxn))<length(orphan_p)){
      reactions <- unique(c(reactions,rxn))
    } else {
      break
    }
  }
  all <-sort(unique(c(reactions,oreactants)))
  if(consensus==TRUE){
    return(all)
  }else{
    new <- all[!all%in%reactionList]
    return(new)
  }
}
