#' RDP training set v19
#'
#' The sequence and taxonomy data for the 24,642 sequences found in the Ribosomal
#' Database Project's trainset19_072023 training set for use with the naive
#' bayesian classifier as implemented in the `{bacteriatypR}` R package.
#' originally released by the RDP in July 2023
#'
#'
#' @format
#' Format data frame with 24,642 rows and three columns. Each row represents a
#' different sequence:
#' \describe{
#'   \item{id}{sequence accession identifier}
#'   \item{sequence}{DNA sequence-string}
#'   \item{taxonomy}{taxonomic string with each level separated with a `;`}
#' }
#'
#' @source
#' - [mothur_formatted files](https://mothur.org/rdp/wiki/rdp_reference_files/)
#' - [Description of how mothur-formatted files were generated](https://mothur.org/blog/2024/RDP-v19-reference-files/)
#' - [RDP sourceforge page](https://sourceforge.net/p/rdp-classifier/news/2023/08/rdp-classifier-214-august-2023-released/)

"trainset19_df"

#' trainset19_db
#'
#' The sequence and taxonomy data for the 24,642 sequences found in the Ribosomal
#' Database Project's trainset19_072023 training set for use with the naive
#' bayesian classifier as implemented in the `{bacteriatypR}` R package.
#' generated using the `bacteriatypR::build_kmer_database()` function:
#'
#' @format A bacteriatypR database with 24,642 sequences that has been trained
#'         using 8-mer kmers for the RDP taxonomy from phylum down to genus:
#' \describe{
#'   \item{conditional_prob}{the conditional probabilities for each genus and kmer}
#'   \item{genera}{taxonomic string with each level separated with a `;`}
#' }
#'
#' @source
#' - [mothur_formatted files](https://mothur.org/rdp/wiki/rdp_reference_files/)
#' - [Description of how mothur-formatted files were generated](https://mothur.org/blog/2024/RDP-v19-reference-files/)
#' - [RDP sourceforge page](https://sourceforge.net/p/rdp-classifier/news/2023/08/rdp-classifier-214-august-2023-released/)

"trainset19_db"
