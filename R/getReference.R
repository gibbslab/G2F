#' @export getReference
#' @importFrom "KEGGREST" "keggLink" "keggList"
#' @title Download all the set of gene-associated stoichiometric reactions for a specific organism from the KEGG database
#' @author Kelly Johana Botero <kjboteroo@unal.edu.co> - Mantainer: Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#  Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
#  Experimental and Computational Biochemistry | Pontificia Universidad Javeriana
#' @description This function downloads all the gene-associated stoichiometric reactions for a given organism from the KEGG database. If not valid organism identifier is given, all reactions from the KEGG database are downloaded. GPR are constructed using the KEGG KO association for each enzyme in a specific organism.
#' @param organism A valid organism identifier for the KEGG database. List of valid organism identifiers are available in: http://rest.kegg.jp/list/organism. If no given, all KEGG stoichiometric reactions are downloaded.
#' @param sep A character string to separate the terms.
#' @param onlyGPR A boolean value to define if just the GPR associated to each EC should be returned
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
#' getReference(organism = "hsa",sep = ";")
#' }
#' @seealso The KEGG database webpage: http://www.genome.jp/kegg/
#' @keywords Download reference genome scale metabolic reconstruction
getReference<-function(organism = "all",sep = ";", onlyGPR = FALSE){
  # Downloading organism
  kegg_download <- tempdir()
  download.file("rest.kegg.jp/list/organism",paste0(kegg_download,"organism.txt"),quiet = TRUE)
  kegg_organism <- as.data.frame.array(read.csv2(paste0(kegg_download,"organism.txt"),header = FALSE,sep ="\t"))
  organism <- match(x = organism,table = kegg_organism[,2])
  ifelse(test = is.na(organism),yes = organism <- "all",no = organism <- kegg_organism[organism,2])
  # Downloading all reactions
  reaction_all <- data.frame(keggList("reaction"))
  id <- as.vector(regmatches(rownames(reaction_all),regexpr("R[[:digit:]]+",rownames(reaction_all))))
  reaction <- as.vector(sapply(as.vector(reaction_all[,1]), .extract))
  # Downloading enzyme association
  ez_all <- keggLink("enzyme","reaction")
  ez_all <- as.data.frame(cbind(id=as.vector(gsub("rn:","",names(ez_all))),ec=as.vector(gsub("ec:","",ez_all))))
  reaction_all <- as.data.frame(cbind(id,reaction))
  reaction_all <- merge(reaction_all,ez_all,by.x="id",by.y = "id",all.x = TRUE)
  # Downloading KO association
  ko_all <- keggLink("ko", "enzyme")
  ko_all <- cbind(ec=gsub("ec:","",names(ko_all)),ko=gsub("ko:","",as.vector(ko_all)))
  reaction_all<-merge(reaction_all,ko_all,by.x = "ec",by.y = "ec",all.x = TRUE)
  # Setting up a specific organism
  if(organism != "all"){
    ko_o <- keggLink(organism, "ko")
    ko_o <- cbind(ko=as.vector(gsub("ko:","",names(ko_o))),unigene=as.vector(regmatches(ko_o,regexpr("[[:digit:]]+",ko_o))))
    ko_o <- as.data.frame(cbind(ko=unique(ko_o[,"ko"]),gpr=sapply(unique(ko_o[,"ko"]), function(ko){paste("(",paste0(ko_o[ko_o[,"ko"]%in%ko,2],collapse = " and "),")")})))
    reaction_all <- merge(reaction_all,ko_o,by.x = "ko",by.y = "ko",all.x = TRUE)
    reaction_all <- reaction_all[reaction_all[,"ko"]%in%ko_o[,"ko"],]
    reaction_all <- as.data.frame.array(reaction_all)
  if(isTRUE(onlyGPR)){
    reaction_all <- unique(reaction_all[,c("ec","gpr")])
    reaction_all[,"gpr"] <- summarize(id = "ec",column = "gpr",table = reaction_all,sep = "or")
    reaction_all <- unique(reaction_all)
    return(reaction_all)
  }
    reaction_all[,"ec"] <- summarize(id = "id",column = "ec",table = reaction_all,sep = ";")
    reaction_all[,"gpr"] <- summarize(id = "id",column = "gpr",table = reaction_all,sep = "or")
    reaction_all[,"ko"] <- summarize(id = "id",column = "ko",table = reaction_all,sep = ";")
  }
  # Extracting uniques
  reaction_all <- unique(reaction_all)
  rownames(reaction_all) <- NULL
  # Return
  return(reaction_all)
}
