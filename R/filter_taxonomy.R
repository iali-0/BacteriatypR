#' filter taxonomy
#'
#' The `filter_taxonomy()` function will filter a consensus taxonomy to remove
#' only taxonomic levels where the confidence score is below a `min_confidence`
#' level
#'
#' @inheritParams print_taxonomy
#' @param min_confidence The minimum fraction of bootstrap replicates that had
#'                       the same classification. Any confidence score below this
#'                       value will have the corresponding taxonomy removed
#'
#' @returns A list object containg two equally sized vectors that are filtered
#'          to remove low confidence taxonomies. One vector, `taxonomy`, contains
#'          the taxonomy at each taxonomic level and the other vector, `confidence`
#'          contains the confidence score for that taxonomy.There will be no
#'          taxonomies or confidence scores below `min_confidence`
#'
#' @examples
#' Oscillospiraceae <- list(taxonomy =c("Bacteria","Bacteroidota","Bacteroidia",
#'                                      "Bacteroidales","Porphyromonadaceae",
#'                                       "Macellibacteroides"),
#'                          confidence =c(1.00, 1.00, 1.00, 1.00, 0.83, 0.75))


#' filter_taxonomy(Oscillospiraceae,min_confidence = 0.80)

#'
#' @export

filter_taxonomy <- function(consensus,min_confidence = 0.80){

  #classification <- Bacteroidales
  high_confidence <- which(consensus$confidence >= min_confidence)
  filtered <- list()
  filtered[["taxonomy"]] <- consensus[["taxonomy"]][high_confidence]
  filtered[["confidence"]] <- consensus[["confidence"]][high_confidence]

  return(filtered)

}
