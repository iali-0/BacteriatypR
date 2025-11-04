#' Read in a FASTA-formatted file containing DNA sequences
#'
#'Given a [standard FASTA-formatted file](https://en.wikipedia.org/wiki/FASTA_format),
#'read_fasta` will read in the contents of the file and create three column data
#' frame with columns for the sequence id, the sequence itself, and any comments
#' found in the header line for each sequence
#'
#'
#' @param file
#'Either a path to a file, a connection or literal data (either a single string
#'or row vector) containing DNA sequences in the standard FASTA format. There
#'are no checks to determine whether the data are DNA or amino acid sequences.
#'
#'files ending in .gz,.bz2, .xz, or .zip will be automatically uncompressed.
#'files starting with http://, ftp://, or ftps:// will be automatically
#'downloaded and decompressed.
#'
#'@note
#'The sequences in the FASTA file can have file breaks within and
#'`read_fasta()` will put those separate lines into the same sequence
#'
#' @returns A data frame object with three columns.The `id` column will contain
#'          the non-space characters following the `>` in the header line of each
#'          sequence; the `sequence` column will contain the sequence; and the
#'          `comment` column will contain any text found after the first white
#'          space character on the header line.

#'
#' @examples
#'temp <- tempfile()
#'write(">seqA\nATGCATGC\n>seqB\nATGCATGC",file = temp)
#'write(">seqC\nTCCGTACT",file = temp,append = TRUE)
#'write(">seqD B.Ceresus UM85\nTCCGTACT",file = temp,append = TRUE)
#'write(">seq4\tE.coli k12\tBacteria;proteobacteria;\nTCCGTACT",file = temp,append = TRUE)
#'write(">seq_4\tSalmonella LT2\tBacteria;proteobacteria;\nTCCGTACT",file = temp,append = TRUE)
#'write(">seqE B.Ceresus UM123\nTCCGTACT\nTCCGTACT",file = temp,append = TRUE)
#'
#'sequence_df <- read_fasta(temp)
#'
#' @importFrom readr read_lines
#' @importFrom stringi stri_startswith_fixed stri_replace_first_regex
#' @export
read_fasta <- function(file){
  fasta_data <- readr::read_lines(file)
  is_header <- stringi::stri_startswith_fixed(fasta_data,">")
  header_lines <- fasta_data[is_header]

  id <- header_lines|>
    stringi::stri_replace_first_regex("^>(\\w*)\\s*.*","$1")

  comment <- header_lines |>
    stringi::stri_replace_first_regex("^>\\w*\\s*(.*)","$1")


  number <- cumsum(is_header)
  seq_lines <- fasta_data[!is_header]
  seq_number <- number[!is_header]
  sequence <- tapply(seq_lines,seq_number,stringi::stri_c, collapse = "")|> unname()


  data.frame(id = id,
             sequence = sequence,
             comment = comment)

}
