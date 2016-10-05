#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

.extract <- function(reaction){
  if (grepl(";",reaction)){
    parts <- unlist(strsplit(reaction,";[[:blank:]]+"))
    return(parts[length(parts)])
  }else{ 
    if (!grepl("rn:R[[:digit:]]+[[:blank:]]",reaction)){
      return(reaction)
    } else {
      return(unlist(strsplit(reaction,"rn:R[[:digit:]]+[[:blank:]]"))[2]) 
    }
  }
}
