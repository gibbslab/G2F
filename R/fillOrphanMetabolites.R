fillOrphanMetabolites <- function(reference, reactionList, limit, reference_metabolites, oM, newR) {
  repeat({
    aC <- additionCost(reactionList = reference, reference = reactionList)
    rA <- reference[aC <= limit]
    toAdd <- rA[unlist(lapply(reference_metabolites[aC <= limit], function(sR) {
      any(sR %in% oM)
    }))]
    if (any(!toAdd %in% newR)) {
      newR <- unique(c(newR, toAdd))
      reactionList <- unique(c(reactionList, newR))
    } else {
      break()
    }
  })
  return(list("new" = newR, "summary" = reactionList))
}