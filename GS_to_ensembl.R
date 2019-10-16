##------------------------------------------
## Author: Matina Fragkogianni
## Date: 16-10-2019
##
## A function to translate Gene symbols to Ensembl IDs using the biomart database 
##
## @param subset is a string vector 
## @return shows in stdout the ensembl IDs and their corresponding gene symbols
##
##------------------------------------------

GS_to_ensembl <- function(subset, species = "Human"){
  
  ##------------------------------------------
  # Load libraries
  ##------------------------------------------
  require(biomaRt)
  
  if(species %in% "Human"){
    ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "grch37.ensembl.org")
    all_genes <- getBM(filters = "external_gene_name", 
                       attributes = c("ensembl_gene_id","external_gene_name", "description","mmusculus_homolog_ensembl_gene",
                                      "mmusculus_homolog_orthology_type", "mmusculus_homolog_perc_id"), 
                       values = subset, mart = ensembl, bmHeader = F )
    
  }else if(species %in% "Mouse"){
    ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
    all_genes <- getBM(filters = "external_gene_name", attributes = c("ensembl_gene_id", "external_gene_name", "description"),
                       values = subset, mart = ensembl)
  }
  
  return(all_genes)
}

              