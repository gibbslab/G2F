#' @export gapFill
#' @importFrom "minval" "orphanMetabolites"
#' @author Daniel Osorio <dcosorioh@unal.edu.co>
gapFind <- function(reactionList){
  orphanMetabolites(reactionList = reactionList)
}