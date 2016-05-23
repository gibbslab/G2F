get.reference<-function(organism){
  reaction_all <- data.frame(keggList("reaction"))
  id_all <- as.vector(regmatches(rownames(reaction_all),regexpr("R[[:digit:]]+",rownames(reaction_all))))
  rx_all <- as.vector(sapply(as.vector(reaction_all[,1]), .extract))
  reaction_all <- cbind(id_all,rx_all)
  if(organism == "all"){
    return (reaction_all)
  } else {
    ko_all <- keggLink("ko", "reaction")
    ko_o <- keggLink(organism, "ko")
    ko_o <- names(ko_all[ko_all%in%names(ko_o)])
    id_o <- unique(regmatches(ko_o,regexpr("R[[:digit:]]+",ko_o)))
    kegg_o <- as.data.frame(reaction_all[reaction_all[,1]%in%id_o,])
    colnames(kegg_o) <- c("id","reaction")
    return(kegg_o)
  }
}