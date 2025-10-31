test_that("read in fasta formatted data generates dataframe", {
  temp <- tempfile()
  write(">seqA\nATGCATGC\n>seqB\nATGCATGC",file = temp)
  write(">seqC\nTCCGTACT",file = temp,append = TRUE)
  write(">seqD B.Ceresus UM85\nTCCGTACT",file = temp,append = TRUE)
  write(">seq4\tE.coli k12\tBacteria;proteobacteria;\nTCCGTACT",file = temp,append = TRUE)
  write(">seq_4\tSalmonella LT2\tBacteria;proteobacteria;\nTCCGTACT",file = temp,append = TRUE)
  write(">seqE B.Ceresus UM123\nTCCGTACT\nTCCGTACT",file = temp,append = TRUE)
  sequence_df <- read_fasta(temp)
  expected <- data.frame(id = c("seqA","seqB","seqC","seqD","seq4","seq_4","seqE"),
                         sequence = c("ATGCATGC","ATGCATGC","TCCGTACT","TCCGTACT","TCCGTACT","TCCGTACT","TCCGTACTTCCGTACT" ),
                         comment = c("","","","B.Ceresus UM85","E.coli k12\tBacteria;proteobacteria;","Salmonella LT2\tBacteria;proteobacteria;","B.Ceresus UM123"))
  expect_equal(sequence_df,expected)


})
