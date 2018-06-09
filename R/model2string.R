model2string <- function(model){
  stoichiometricMatrix <- model@S
  rownames(stoichiometricMatrix) <- model@met_id
  sReactions <- minval:::rearmReactions(S = stoichiometricMatrix, reversible = model@react_rev)
  return(sReactions)
}