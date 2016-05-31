associated.data <- function(reactionList,reference,reactions){
  return(reference[reference[,reactions]%in%reactionList,])
}
