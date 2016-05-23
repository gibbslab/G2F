g2f <- function(reactionList, reference, limit = 0.5, woCompartment=FALSE,consensus=FALSE){
  if(woCompartment==TRUE){
    reactionList <- gsub("\\[[[:alpha:]]\\]","",reactionList)
  }
  reactions <- reactionList
  orphan <- orphan.reactants(reactions)
  repeat{
    orphan_r <- orphan.reactants(reactions)
    orphan_r <- orphan_r[orphan_r%in%orphan]
    message(paste0(length(orphan_r)," Orphan reactants found"))
    to.add <- unique(unlist(sapply(orphan_r, function(metabolite){.reactant.fill(metabolite,reference)})))
    rxn <- unique(c(reactions,to.add[unlist(sapply(to.add, function(reaction){.addition.cost(reaction,reactions)}))<limit]))
    if(length(orphan.reactants(rxn))<length(orphan_r)){
      reactions <- unique(c(reactions,rxn))
    } else {
      break
    }
  }
  oreactants <- reactions
  reactions <- reactionList
  orphan <- orphan.products(reactions)
  repeat{
    orphan_p <- orphan.products(reactions)
    orphan_p <- orphan_p[orphan_p%in%orphan]
    message(paste0(length(orphan_p)," Orphan products found"))
    to.add <- unique(unlist(sapply(orphan_p, function(metabolite){.product.fill(metabolite,reference)})))
    rxn <- unique(c(reactions,to.add[unlist(sapply(to.add, function(reaction){.addition.cost(reaction,reactions)}))<limit]))
    if(length(orphan.products(rxn))<length(orphan_p)){
      reactions <- unique(c(reactions,rxn))
    } else {
      break
    }
  }
  all <-sort(unique(c(reactions,oreactants)))
  all <- gsub("^[[:space:]]","",all)
  all <- gsub("[[:space:]]$","",all)
  if(consensus==TRUE){
    return(all)
  }else{
    new <- all[!all%in%reactionList]
    return(new)
  }
}
