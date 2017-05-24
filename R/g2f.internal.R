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

summary.ec <- function(id){
  data <- reaction_all[reaction_all[,"id"]%in%id,]
  data[,"ec"] <- paste0(unique(data[,"ec"]),collapse = sep)
}

summarize <- function(id, column, table, sep = "or"){
  sapply(table[,id],function(index){
    index_table <- table[table[,id] %in% index,]
    paste0(unique(index_table[,column]),collapse = paste0(" ",sep," "))
  })
}
