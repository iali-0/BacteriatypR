#'Build kmer database
#'
#' Build kmer database for classifiying 16s rRNA and other gene sequences to
#' a genus when a kmer size is provided.
#'
#' @param sequences A vector of reference sequences for which we have
#'                  genus-level taxonomic information in the same order as the
#'                  value for genera
#' @param genera A vector of genus-level taxonomic information for reference
#'               sequences in the same order as the value sequences.
#'               Ideally, taxonomic information will be provided back to the
#'               domain level with each level separated by semicolons and no
#'               spaces.
#'
#' @param kmer_size The length of the nucleotide word to base our classification
#'                  on (default = 8).
#'
#' @returns A list object containing the genus level conditional probability
#'          (`conditional_prob`) of seeing each kmer in a given genus as well
#'          as genus names (`genera`).
#'
#' @references
#' Wang Q, Garrity GM, Tiedje JM, Cole JR. Naive Bayesian classifier for rapid
#' assignment of rRNA sequences into the new bacterial taxonomy. Appl Environ
#' Microbiol. 2007 Aug;73(16):5261-7. [doi: 10.1128/AEM.00062-07](https://dx.doi.org/10.1128/AEM.00062-07).
#'
#'
#' @examples
#' kmer_size <- 3
#' sequences <- c("ATGCGCTA","ATGCGCTC","ATGCGCTC")
#' genera <- c("A", "B", "B")
#' build_kmer_database(sequences, genera, kmer_size)
#'
#' @export
build_kmer_database <- function(sequences, genera, kmer_size = 8){

  genera_indices <- genera_str_to_index(genera)

  detected_kmers <- detect_kmers_across_sequences(sequences, kmer_size = kmer_size)


  priors <- calcul_word_specific_priors(detected_kmers,kmer_size)

  cond_prob <- calc_genus_conditional_prob(detected_kmers,genera_indices, priors)
  genera_names <- get_unique_genera(genera)
  return(list(conditional_prob = cond_prob, genera = genera_names))
}

#' Classify 16s rRNA gene sequence fragment
#'
#' The `classify_seq()` function implements the Wong et al. naive Bayesian
#' classification algorithm for 16s rRNA gene sequence.
#'
#' @param unknown_sequence A DNA sequence that needs to be classified
#' @param database A kmer database generated using `build_kmer_database`
#' @param kmer_size An integer value (default of 8) indicating the size of kmer
#'                  to use for classifying sequences.higher value more RAM
#'                  with potentially more specificity lower value use less RAM.
#'                   Benchmarking has found that the default of 8 provides the
#'                   best specificity with the lowest possible memory requirement
#'                   and fastest execution time.
#'
#' @param num_bootstraps An integer value (default of 100). The value of
#'                       `num_bootstraps` is the number of randomization
#'                       to perform where `1/kmer_size` of all kmers are sampled
#'                       (without replacement) from `unknown sequence`. Higher
#'                       value will provide greater precision on the confidence
#'                       score.
#'
#'
#' @returns A list object of two vectors. one vector `taxonomy`
#'          is the taxonomic assignment for each level. the secod
#'          vector `confidence` is the fraction of `num_bootstraps`
#'          that the classifier gave the same classification at that level.
#'
#' @inherit build_kmer_database references
#'
#'
#' @export
#'
#' @examples
#' kmer_size <- 3
#' sequences <- c("ATGCGCTA","ATGCGCTC","ATGCGCTC")
#' genera <- c("A", "B", "B")
#' db <- build_kmer_database(sequences, genera, kmer_size)
#' unknown_sequence <- "ATGCGCTC"

#' classify_sequence(unknown_sequence,
#'                   database = db,
#'                   kmer_size = kmer_size)
#'
#'

classify_sequence <- function(unknown_sequence,database,
                              kmer_size = 8, num_bootstraps = 100){

  kmers <- detect_kmers(sequence = unknown_sequence,kmer_size)
  bs_class <- numeric(length = num_bootstraps)
  for(i in 1:num_bootstraps){
    bs_kmers <- bootstrap_kmers(kmers,kmer_size)
    bs_class[[i]] <- classify_bs(bs_kmers,database)
  }
  consensus_bs_class( bs_class,database)
}



#' @noRd
#' @importFrom stringi stri_length stri_sub
get_all_kmers <- function(x,kmer_size = 8){
  seq_lenght <- stringi::stri_length(x)
  n_kmers <- seq_lenght - kmer_size + 1
  seq_kmers <- stringi::stri_sub(x, 1:n_kmers,kmer_size:seq_lenght)
  return(seq_kmers)
}

#' @noRd
#' @importFrom stringi stri_trans_toupper stri_replace_all_charclass stri_trans_char
seq_to_base4 <- function(sequence){
  stringi::stri_trans_toupper(sequence)|>
  stringi::stri_replace_all_charclass(str = _,
                                      pattern = "[^ACGT]",
                                      replacement = "N")|>
    stringi::stri_trans_char(str = _, pattern = "ACGT",replacement = "0123")
}

#' @noRd
#' @importFrom stats na.omit
base4_to_index <- function(base4_str){
  # I need output to be index to start at position 1 rather than 0
  # therefore we add 1 to all base10 values
  stats::na.omit(strtoi(base4_str, base = 4) +1) |> as.numeric()

}

#' @noRd
detect_kmers <- function(sequence, kmer_size = 8){
  seq_to_base4(sequence) |>
    get_all_kmers(kmer_size) |>
    base4_to_index()
}

#' @noRd
detect_kmers_across_sequences <- function(sequences, kmer_size = 8){
  n_sequences <- length(sequences)
  kmer_list <- vector(mode = "list", length = n_sequences)
  for(i in seq_along(1:n_sequences)){
    kmer_list[[i]] <- detect_kmers(sequences[[i]],kmer_size = kmer_size)
  }
  return(kmer_list)

}

#' @noRd
calcul_word_specific_priors <- function(detect_list, kmer_size){
  #n_seqs <- length(detect_list) #optimise the code
  #n_kmers <- 4^kmer_size # optimise the code
  priors <- detect_list |> unlist() |> tabulate(bin = _, nbins = 4^kmer_size)
  (priors +0.5)/(length(detect_list) + 1)
}

#' @noRd
calc_genus_conditional_prob <- function(detect_list,genera,#need to be an interger
                                        calcul_word_specific_priors){
  genus_counts <- tabulate(genera)#|> as.vector()
  n_genera <- length(genus_counts)
  n_sequences <- length(genera)
  n_kmers <- length(calcul_word_specific_priors)
  kmer_genus_count <- matrix(0,nrow = n_kmers,ncol= n_genera)
  for(i in 1:n_sequences){
    kmer_genus_count[detect_list[[i]],genera[i]] <- kmer_genus_count[detect_list[[i]],genera[i]]+1

  }
  #log(t(t(kmer_genus_count + calcul_word_specific_priors)/(genus_counts + 1)))# 10.2 sec
  #log((kmer_genus_count + calcul_word_specific_priors) %*% diag(1/(genus_counts + 1)))# forever
  #log(sweep(kmer_genus_count + calcul_word_specific_priors,2,genus_counts + 1,"/"))#10.3sec

  log((kmer_genus_count + calcul_word_specific_priors)/rep(genus_counts + 1,each = n_kmers))#6.2sec
  #log((kmer_genus_count + calcul_word_specific_priors)/t(replace(genus_counts),
                                                            #TRUE,genus_counts +1))# 10.8sec
  #(kmer_genus_count + calcul_word_specific_priors)/(genus_counts +1)[col(genus_counts)]# 8.9sec
  #calculate_log_probability(kmer_genus_count,calcul_word_specific_priors,genus_counts)
}


#' @noRd
genera_str_to_index <- function(string){
  factor(string) |> as.numeric()

}

#' @noRd
get_unique_genera <- function(string){
  factor(string)|> levels()

}
#' @noRd
bootstrap_kmers <- function(kmers,kmer_size = 8){
  n_kmers <- as.integer(length(kmers)/kmer_size)
  sample(kmers,n_kmers, replace = TRUE)
}
#' @noRd
#' @importFrom Rfast colsums
classify_bs <- function(unknown_kmers,db){
  probabilities <- Rfast::colsums(db$conditional_prob[unknown_kmers,])
  which.max(probabilities)
}
#' @noRd
consensus_bs_class <- function(bs_class,db){
  taxonomy <- db[["genera"]][bs_class]
  taxonomy_split <- stringi::stri_split_fixed(taxonomy,pattern = ";")
  n_levels <- length(taxonomy_split[[1]])
  consensus_list <- lapply(1:n_levels,
                           \(i) sapply(taxonomy_split,
                                       \(p)paste(p[1:i],collapse = ";"))|>
                             get_consensus()
    )
  list(taxonomy = stringi::stri_split_fixed(consensus_list[[n_levels]][["id"]],
                                            pattern = ";")|>
         unlist(),
         confidence = sapply(consensus_list,'[[',"frac"))
}

#' @noRd
get_consensus <- function(taxonomy){
  n_bs <- length(taxonomy)
  taxonomy_table <- table(taxonomy)
  max_index <- which.max(taxonomy_table)
  list(frac = taxonomy_table[[max_index]]/n_bs,
       id = names(max_index))
}

