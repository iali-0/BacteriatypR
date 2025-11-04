#' Print taxonomy for an unknown sequence
#'
#' The `print_taxonomy()` will output the cosensus taxonomy for an unknown
#' sequence will confidence scores for each taxonomic level and each
#' taxonomic level separated by semi-colons
#'
#' @param consensus A list object that contains two slots each with
#'                  an equal sized vector. the `taxonomy` vector contains
#'                  the classification at each taxonomic level and the
#'                  `confidence` vector contains the fraction fo bootstraps
#'                  that the specified classification
#'
#' @param n_levels An integer that indicates the number of taxonomic levels to
#'                 expect for each sequence (default = 6)
#'
#' @returns A character string indicating the classification at each taxonomic
#'          level with the corresponding confidence in parentheses.
#'          Each taxonomic level is separated by a semi-colon
#'
#' @examples
#' filtered <- list(taxonomy = c("Bacteria","Bacteroidota",
#'                  "Bacteroidia","Bacteroidales","Porphyromonadaceae"),
#'                   confidence = c(1.00,1.00,1.00,1.00,0.83))

#' tax_string <- print_taxonomy(filtered,n_levels = 6)
#'
#' @export
print_taxonomy <- function(consensus,n_levels = 6){

  original_levels <-  length(consensus$taxonomy)
  given_levels <- original_levels
  while(given_levels < n_levels){
    consensus$taxonomy[given_levels +1] <- paste(consensus$taxonomy[original_levels],
                                                 "unclassified",sep = "_")
    consensus$confidence[given_levels +1] <- consensus$confidence[original_levels]
    given_levels <- given_levels+1
  }
  pretty_confidence <- paste0("(",100*consensus$confidence,")")
  paste(consensus$taxonomy,pretty_confidence,sep = "", collapse = ";")#to rid off spaces
}
