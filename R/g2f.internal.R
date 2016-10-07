#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#' @importFrom "utils" "download.file" "read.csv2" "setTxtProgressBar" "txtProgressBar"
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