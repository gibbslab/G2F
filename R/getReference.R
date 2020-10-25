#' @export getReference
#' @importFrom "KEGGREST" "keggLink" "keggList" "keggGet"
#' @title Download all the set of gene-associated stoichiometric reactions for a specific organism from the KEGG database
#' @author Kelly Johana Botero <kjboteroo@unal.edu.co> - Mantainer: Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @description This function downloads all the gene-associated stoichiometric reactions for a given organism from the KEGG database. If not valid organism identifier is given, all reactions from the KEGG database are downloaded. GPR are constructed using the KEGG KO association for each enzyme in a specific organism.
#' @param organism A valid organism identifier for the KEGG database. List of valid organism identifiers are available in: http://rest.kegg.jp/list/organism. If no given, all KEGG stoichiometric reactions are downloaded.
#' @param sep A character string to separate the terms.
#' @return A data.frame with the following data associated to the stoichiometric reactions for a given organism: \itemize{
#' \item{\code{ko:}} The associated KEGG KO identifier to the reaction. In KEGG, molecular-level functions are stored in the KO (KEGG Orthology) database and associated with ortholog groups in order to enable extension of experimental evidence in a specific organism to other organisms.
#' \item{\code{id:}} The associated reaction id from the KEGG database.
#' \item{\code{reaction:}} The gene-associated stoichiometric reactions with the following format:
#'
#' \code{"H2O + Urea-1-carboxylate <=> 2 CO2 + 2 NH3"}
#'
#' Where arrows and plus signs are surrounded by a "space character", and stoichiometry coefficients are surrounded by spaces, (nothe the "2" before the CO2 or the NH3). Arrows will be in the form "\code{=>}" or "\code{<=>}". KEGG reactions are not compartmentalized.
#' \item{\code{gpr:}} The Gene-Protein-Reaction (GPR) associations for a specific organism buit from the KEGG KO identifiers.
#' }
#' @examples
#' \dontrun{
#' getReference(organism = "hsa", sep = ";")
#' }
#' @seealso The KEGG database webpage: http://www.genome.jp/kegg/
#' @keywords Download reference genome scale metabolic reconstruction
getReference <- function(organism = "all", sep = ";") {
  # Downloading organism
  if (organism != "all") {
    message("Validating selected organism ... ", appendLF = FALSE)
    kegg_download <- tempdir()
    download.file("rest.kegg.jp/list/organism", paste0(kegg_download, "organism.txt"), quiet = TRUE, method = "libcurl")
    kegg_organism <- as.data.frame.array(read.csv2(paste0(kegg_download, "organism.txt"), header = FALSE, sep = "\t"))
    organism <- match(x = organism, table = kegg_organism[, 2])
    if (is.na(organism)) {
      stop("Organism not available, choose one from http://rest.kegg.jp/list/organism")
    } else {
      organism <- kegg_organism[organism, 2]
    }
  }
  message("OK")
  # Downloading all reactions
  message("Downloading reactions ... ", appendLF = FALSE)
  reaction_all <- data.frame(keggList("reaction"))
  id <- as.vector(regmatches(rownames(reaction_all), regexpr("R[[:digit:]]+", rownames(reaction_all))))
  reaction <- as.vector(sapply(as.vector(reaction_all[, 1]), extract))
  message("DONE")
  # Downloading enzyme association
  message("Downloading enzymes ... ", appendLF = FALSE)
  ez_all <- keggLink("enzyme", "reaction")
  ez_all <- as.data.frame(cbind(id = as.vector(gsub("rn:", "", names(ez_all))), ec = as.vector(gsub("ec:", "", ez_all))))
  reaction_all <- as.data.frame(cbind(id, reaction))
  reaction_all <- merge(reaction_all, ez_all, by.x = "id", by.y = "id", all.x = TRUE)
  message("DONE")
  # Downloading KO association
  message("Downloading KO ... ", appendLF = FALSE)
  ko_all <- keggLink("ko", "reaction")
  ko_all <- cbind(id = gsub("rn:", "", names(ko_all)), ko = gsub("ko:", "", as.vector(ko_all)))
  reaction_all <- merge(reaction_all, ko_all, by.x = "id", by.y = "id", all.x = TRUE)
  message("DONE")
  # Setting up a specific organism
  if (organism != "all") {
    message("Building organism specific reference ... ")
    ko_o <- keggLink(organism, "ko")
    ko_o <- cbind(ko = as.vector(gsub("ko:", "", names(ko_o))), unigene = gsub("^[[:print:]]+:", "", ko_o))
    ko_o <- as.data.frame(cbind(ko = unique(ko_o[, "ko"]), gpr = sapply(unique(ko_o[, "ko"]), function(ko) {
      paste("(", paste0(ko_o[ko_o[, "ko"] %in% ko, 2], collapse = " and "), ")")
    })))
    reaction_all <- merge(reaction_all, ko_o, by.x = "ko", by.y = "ko", all.x = TRUE)
    reaction_all <- reaction_all[reaction_all[, "ko"] %in% ko_o[, "ko"], ]
    summary.ko <- function(id) {
      data <- reaction_all[reaction_all[, "id"] %in% id, ]
      data[, "ko"] <- paste0(unique(data[, "ko"]), collapse = sep)
    }
    summary.gpr <- function(id) {
      data <- reaction_all[reaction_all[, "id"] %in% id, ]
      data[, "gpr"] <- paste0(unique(data[, "gpr"]), collapse = " or ")
    }
    reaction_all[, "ko"] <- sapply(as.vector(reaction_all[, "id"]), summary.ko)
    reaction_all[, "gpr"] <- sapply(as.vector(reaction_all[, "id"]), summary.gpr)
    reaction_all <- unique(reaction_all)
    # Downloading directions
    message("Setting directionality ... ", appendLF = FALSE)
    direction <- lapply(names(keggList("pathway", organism)), function(x) {
      suppressMessages(keggGet(x, "kgml"))
    })
    direction <- lapply(direction, function(b) {
      b <- unlist(strsplit(b, "\n"))
      b <- b[grepl("\\<reaction[[:space:]]+", b)]
      b <- do.call(rbind.data.frame, strsplit(b, '\"'))[, c(4, 6)]
      colnames(b) <- c("rxn", "rev")
      return(b)
    })
    direction <- do.call(rbind.data.frame, direction)
    direction <- apply(direction, 1, function(x) {
      expand.grid(unlist(strsplit(x[[1]], "[[:space:]]+")), x[[2]])
    })
    direction <- unique(do.call(rbind.data.frame, direction))
    direction[, 1] <- gsub("rn:", "", direction[, 1])
    tdirection <- table(direction[, 1])
    direction[tdirection > 1, 2] <- "reversible"
    direction <- unique(direction)
    direction <- direction[direction[, 2] == "irreversible", 1]
    reaction_all <- as.data.frame.array(reaction_all)
    reaction_all$reaction[reaction_all$id %in% direction] <- gsub("<=>", "=>", reaction_all$reaction[reaction_all$id %in% direction])
    message("DONE")
  }
  summary.ec <- function(id) {
    data <- reaction_all[reaction_all[, "id"] %in% id, ]
    data[, "ec"] <- paste0(unique(data[, "ec"]), collapse = sep)
  }
  reaction_all[, "ec"] <- sapply(as.vector(reaction_all[, "id"]), summary.ec)
  # Extracting uniques
  reaction_all <- unique(reaction_all)
  rownames(reaction_all) <- NULL
  # Return
  message("DONE")
  return(reaction_all)
}
