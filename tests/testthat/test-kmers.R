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
test_that("can extract speicifickmer from a starting position and size",{
  x <- "ATGCGCTAGGCATGCTG"
  kmer <- get_kmer(x,1,kmer_size = 8)
  expect_equal(kmer,"ATGCGCTA")
  kmer <- get_kmer(x,5,kmer_size = 8)
  expect_equal(kmer,"GCTAGGCA")
  kmer <- get_kmer(x,10,kmer_size = 8)
  expect_equal(kmer,"GCATGCTG")
  expect_error(get_kmer(x,11,kmer_size))
  kmer <- get_kmer(x,5)
  expect_equal(kmer,"GCTAGGCA")
})
