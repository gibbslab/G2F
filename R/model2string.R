model2string <- function(model){
  stoichiometricMatrix <- model@S
  rownames(stoichiometricMatrix) <- model@met_id
  return(model <- minval:::rearmReactions(S = stoichiometricMatrix, reversible = model@react_rev))
}