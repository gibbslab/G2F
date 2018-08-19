transferReact <- function(rxnIndex, model1, model2){
  sapply(rxnIndex, function(rxn){
    model2 <<- sybil::addReact(model = model2, 
                               id = model1@react_id[rxn], 
                               met = model1@met_id[model1@S[,rxn] != 0], 
                               Scoef = model1@S[model1@S[,rxn] != 0,rxn], 
                               reversible = model1@react_rev[rxn], 
                               lb = model1@lowbnd[rxn], 
                               ub = model1@uppbnd[rxn], 
                               subSystem = colnames(model1@subSys)[model1@subSys[rxn,]], 
                               gprAssoc = model1@gpr[rxn], 
                               reactName = model1@react_name[rxn], 
                               metName = model1@met_name[model1@S[,rxn]!=0], 
                               metComp = model1@mod_compart[model1@met_comp[model1@S[,rxn]!=0]]
    )
  })
  return(model2)
}

