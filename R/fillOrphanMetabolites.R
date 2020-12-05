fillOrphanMetabolites <- function(reference, reactionList, limit, reference_metabolites, oM, newR) {
  repeat({
    aC <- additionCost(reactionList = reference, reference = reactionList)
    rA <- reference[aC <= limit]
    whichAdd <- unlist(lapply(reference_metabolites[aC <= limit], function(sR) {
      any(sR %in% oM)
    }))
    toAdd <- rA[whichAdd]
    aC <- aC[aC <= limit][whichAdd]
    if (any(!toAdd %in% newR$react)) {
      newR <- rbind(newR, data.frame("addCost" = aC, "react" = toAdd))
      newR <- newR[!duplicated(newR$react), ]
      reactionList <- unique(c(reactionList, newR$react))
    } else {
      break()
    }
  })
  return(list("new" = newR, "summary" = reactionList))
}
