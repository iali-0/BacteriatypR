test_that("apply filter to confidence score",{
  Oscillospiraceae <- list(taxonomy =c("Bacteria","Bacteroidota",
                                       "Bacteroidia","Bacteroidales",
                                       "Porphyromonadaceae", "Macellibacteroides"),
                           confidence =c(1.00, 1.00, 1.00, 1.00, 0.83, 0.75))

  filtered <- list(taxonomy =c("Bacteria","Bacteroidota",
                               "Bacteroidia","Bacteroidales",
                               "Porphyromonadaceae"),
                   confidence =c(1.00, 1.00, 1.00, 1.00, 0.83))

  observed <- filter_taxonomy(Oscillospiraceae,min_confidence = 0.80)
  expect_equal(observed,filtered)

})
