#' @export getReactionsList
#' @importFrom "sybil" "checkReactId"
#' @author Nicolas Mendoza Mejia <nimendozam@unal.edu.co>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @title extract reactions from a sybil model as a list of strings
#' @description This function is based in printReaction function from sybil, it
#' extracts all reactions from a sybil model as a list of strings with the fromat:
#' \code{"H2O[c] + Urea-1-carboxylate[c] <=> 2 CO2[c] + 2 NH3[c]"}
#'
#' Where arrows and plus signs are surrounded by a "space character".
#' It is also expected that stoichiometry coefficients are surrounded by spaces, (nothe the "2" before the CO2[c] or the NH3[c]).
#' It also expects arrows to be in the form "\code{=>}" or "\code{<=>}".
#' Meaning that arrows like "\code{==>}", "\code{<==>}", "\code{-->}" or "\code{->}" will not be parsed and will lead to errors.
#'
#' @seealso \code{additionCost} function documentation.
#' @param model A metabolic model in the form of "modelorg" class from sybil package
#'
#' @examples
#' \dontrun{
#' # Get the sybil model
#' library("sybil")
#' data("Ec_core")
#'
#' # Extract reactions
#' react <- getReactionsList(Ec_core)
#' }
getReactionsList <- function(object) {
  check <- checkReactId(object, react = object@react_id)
  cind <- react_pos(check)

  mat <- S(object)[, cind, drop = FALSE]
  nnz <- apply(mat, 2, "!=", 0)
  reaction <- character(length(cind))
  id <- character(length(cind))

  for (j in seq(along = cind)) {
    met <- met_id(object)[nnz[, j]]
    nzv <- mat[, j][nnz[, j]]

    ed <- nzv < 0
    pd <- nzv > 0

    if (sum(ed) > 0) {
      educt <- paste(met[ed], collapse = " + ")
    }
    else {
      educt <- ""
    }

    if (sum(pd) > 0) {
      product <- paste(met[pd], collapse = " + ")
    }
    else {
      product <- ""
    }

    arrow <- ifelse(react_rev(object)[cind[j]], " <=> ", " => ")

    reaction[j] <- paste(educt, product, sep = arrow)
    id[j] <- react_id(check)[j]
  }
  return(data.frame("id" = id, "react" = reaction))
}
