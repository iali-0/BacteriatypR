test_that("can extract all possible 8-mers from a sequence",{
  x <- "ATGCGCTAGGCATGCTG"

  all_kmers <- get_all_kmers(x,kmer_size = 8)

  expected_kmers <- c("ATGCGCTA","TGCGCTAG","GCGCTAGG","CGCTAGGC",
                      "GCTAGGCA","CTAGGCAT","TAGGCATG","AGGCATGC",
                      "GGCATGCT","GCATGCTG")

  all_kmers <- get_all_kmers(x)
  expect_equal(all_kmers,expected_kmers)
  expect_length(all_kmers,length(expected_kmers))


    all_kmers <- get_all_kmers(x, kmer_size = 9)
  expected_kmers <- c("ATGCGCTAG","TGCGCTAGG","GCGCTAGGC","CGCTAGGCA",
                      "GCTAGGCAT","CTAGGCATG","TAGGCATGC","AGGCATGCT",
                      "GGCATGCTG")
  expect_equal(all_kmers,expected_kmers)
  expect_length(all_kmers, length(expected_kmers))

})
test_that("can extract speicific kmer from a starting position and size",{
  x <- "ATGCGCTAGGCATGCTG"
  kmer <- get_kmer(x,1,kmer_size = 8)
  expect_equal(kmer,"ATGCGCTA")
  kmer <- get_kmer(x,5,kmer_size = 8)
  expect_equal(kmer,"GCTAGGCA")
  kmer <- get_kmer(x,10,kmer_size = 8)
  expect_equal(kmer,"GCATGCTG")
  expect_error(get_kmer(x,11,kmer_size=8))
  kmer <- get_kmer(x,5)
  expect_equal(kmer,"GCTAGGCA")
})
test_that("conversion between DNA sequence and quanternary",{
  x <-      "ATGCGCTAGGCATGCTG"
  actual <- "03212130221032132"
 base4_seq <- seq_to_base4(x)
 expect_equal(base4_seq,actual)

 x <-      "ATGCGCTRGGCATGCTG"
 actual <- "0321213N221032132"
 base4_seq <- seq_to_base4(x)
 expect_equal(base4_seq,actual)

 x <-      tolower("ATGCGCTRGGCATGCTG")
 actual <- "0321213N221032132"
 base4_seq <- seq_to_base4(x)
 expect_equal(base4_seq,actual)
})
test_that("generate base10 value from kmers in base4",{
  x <- "0000"
  expected <- 1
  actual <- base4_to_index(x)
  expect_equal(expected, actual)

  x <- "1000"
  expected <- 65
  actual <- base4_to_index(x)
  expect_equal(expected, actual)

  x <- "0123"
  expected <- 28
  actual <- base4_to_index(x)
  expect_equal(expected, actual)
  x <- c("0123","1000","0000")
  expected <- c(28,65,1)
  actual <- base4_to_index(x)
  expect_equal(expected,actual)

  x <- c("0123","1000","0000", "000N")
  expected <- c(28,65,1)
  actual <- base4_to_index(x)
  expect_equal(expected,actual)

})
test_that("correctly detect kmers from a sequence",{
  sequence <- "03212130221032132"
  kmers <- get_all_kmers(sequence)
  indices <- base4_to_index(kmers)
  detected <- detect_kmers(sequence)

  expect_equal(length(detected[indices]),length(indices))
  expect_equal(length(detected[-indices]),4^8 - length(indices))
  expect_equal(sum(detected[indices]),length(indices))

  sequence <- "03212130221032132N"
  kmers <- get_all_kmers(sequence)
  indices <- base4_to_index(kmers)
  detected <- detect_kmers(sequence)

  expect_equal(length(detected[indices]),length(indices))
  expect_equal(length(detected[-indices]),4^8 - length(indices))
  expect_equal(sum(detected[indices]),length(indices))

  sequence <- "03212130221032132N"
  kmers <- get_all_kmers(sequence, kmer_size = 7)
  indices <- base4_to_index(kmers)
  detected <- detect_kmers(sequence, kmer_size = 7)

  expect_equal(length(detected[indices]),length(indices))
  expect_equal(length(detected[-indices]),4^7 - length(indices))
  expect_equal(sum(detected[indices]),length(indices))
})
test_that("correctly detect kmers across multiple sequence",{
  kmer_size <- 3
  sequences <- c("03212130","03212131")
  expected <- matrix(0,nrow = 4^kmer_size, ncol = 2)
  expected[base4_to_index(get_all_kmers(sequences[1],kmer_size)),1] <- 1
  expected[base4_to_index(get_all_kmers(sequences[2],kmer_size)),2] <- 1

  detect_matrix <- detect_kmers_across_sequences(sequences, kmer_size)
  expect_equal(detect_matrix, expected)

})

test_that("calculate word specific priors probabilities",{
  kmer_size <- 3
  sequences <- c("03212130","03212131","03212131")
  detect_matrix <- detect_kmers_across_sequences(sequences, kmer_size)

  expected <- (apply(detect_matrix,1,sum) + 0.5)/(1 + length(sequences))
  priors <- calcul_word_specific_priors(detect_matrix)
  expect_equal(priors,expected)
  expect_equal(priors[26],0.875)
  expect_equal(priors[29],0.375)
  expect_equal(priors[30],0.625)
  expect_equal(priors[64],0.125)


  })

test_that("calculate genus-specific conditional probablities",{
  kmer_size <- 3
  sequences <- c("03212130","03212131","03212131")
  genera <- c(1, 2, 2)
  detect_matrix <- detect_kmers_across_sequences(sequences, kmer_size)
  priors <- calcul_word_specific_priors(detect_matrix)

  conditional_prob <- calc_genus_conditional_prob(detect_matrix, genera, priors)

  expect_equal(conditional_prob[26,],(c(1,2)+0.875)/(c(1,2)+1))
  expect_equal(conditional_prob[29,],(c(1,0)+0.375)/(c(1,2)+1))
  expect_equal(conditional_prob[30,],(c(0,2)+0.625)/(c(1,2)+1))
  expect_equal(conditional_prob[64,],(c(0,0)+0.125)/(c(1,2)+1))
})

test_that("Create kmer database from sequences,taxonomy, and kmer_size ",{
  kmer_size <- 3
  sequences <- c("ATGCGCTA","ATGCGCTC","ATGCGCTC")
  genera <- c("A", "B", "B")
  db <- build_kmer_database(sequences, genera, kmer_size)

  expect_equal(db[["conditional_prob"]][26,],(c(1,2)+0.875)/(c(1,2)+1))
  expect_equal(db[["conditional_prob"]][29,],(c(1,0)+0.375)/(c(1,2)+1))
  expect_equal(db[["conditional_prob"]][30,],(c(0,2)+0.625)/(c(1,2)+1))
  expect_equal(db[["conditional_prob"]][64,],(c(0,0)+0.125)/(c(1,2)+1))

  expect_equal(db[["genera"]][[1]],"A")
  expect_equal(db[["genera"]][[2]],"B")
})
test_that("convert back and forth between genus and indices",{
  genera_str <- c("A","B","B")
  genera_index <- c(1,2,2)
  expect_equal(genera_str_to_index(genera_str),genera_index)
  expect_equal(get_unique_genera(genera_str),c("A","B"))
})
