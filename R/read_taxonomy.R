#' Read in taxonomy files
#'
#'Read a [mothur-formatted taxonomy file](https://mothur.org/wiki/constaxonomy_file/)
#' into as a data frame
#'
#' @param file
#' Either a path to a file, a connection, or literal data
#' (either a single string or row vector) containing the sequence id
#' and the taxonomy information for each sequence.
#'
#'files ending in .gz,.bz2, .xz, or .zip will be automatically uncompressed.
#'files starting with http://, ftp://, or ftps:// will be automatically
#'downloaded and decompressed.
#'
#' @returns A data frame with two columns.The `id` column contains
#'          the name for each sequence and the `taxonomy` column is a series of
#'          taxonomic names separated by semi-colon.The string does not have a
#'          semi-colon at the end of the sequence.
#' @note
#' There are no checks to insure that each sequence has a unique
#' id value. it is assumed that each sequence has the same number
#' of taxonomic levelsrepresented in the second column of the input file.
#'
#'
#' @examples
#' temp <- tempfile()
#' write("seqA\tA;B;C;",file = temp)
#' write("seqB\tA;B; C;",file = temp, append = TRUE)
#' write("seqC\tA; B;C;",file = temp, append = TRUE)
#' write("seqD\tA;B;C",file = temp, append = TRUE)
#' write("seqE\tA;B; C",file = temp, append = TRUE)
#' write("seqF\tA; B;C",file = temp, append = TRUE)
#' write("seq G\tA;B;C;",file = temp, append = TRUE)

#' read_taxonomy(temp)
#'
#' @importFrom stringi  stri_replace_last_regex stri_replace_all_regex
#' @importFrom readr read_tsv cols col_character
#' @export
read_taxonomy <- function(file){

  taxonomy <- readr::read_tsv(file,
                  col_names = c("id","taxonomy"),
                  col_types = readr::cols(.default = readr::col_character()))
  taxonomy$taxonomy <- stringi::stri_replace_last_regex(taxonomy$taxonomy,";$", "")
  taxonomy$taxonomy <- stringi::stri_replace_all_regex(taxonomy$taxonomy,"; ", ";")

  as.data.frame(taxonomy)
}
