test_that("print out consensus taxonomy",{
  filtered <- list(taxonomy = c("Bacteria","Bacteroidota",
                                "Bacteroidia","Bacteroidales","Porphyromonadaceae"),
                   confidence = c(1.00,1.00,1.00,1.00,0.83))
  expected <- "Bacteria(100);Bacteroidota(100);Bacteroidia(100);Bacteroidales(100);Porphyromonadaceae(83);Porphyromonadaceae_unclassified(83)"
  tax_string <- print_taxonomy(filtered,n_levels = 6)
  expect_equal(tax_string,expected)


  oscillospiraceae <- list(taxonomy = c("Bacteria","Bacteroidota",
                                        "Bacteroidia","Bacteroidales","Porphyromonadaceae"),
                           confidence = c(1.00,1.00,1.00,1.00,0.83))
  expected <- "Bacteria(100);Bacteroidota(100);Bacteroidia(100);Bacteroidales(100);Porphyromonadaceae(83);Porphyromonadaceae_unclassified(83)"
  tax_string <- print_taxonomy(oscillospiraceae)
  expect_equal(tax_string,expected)

  Bacteroidales <- list()
  Bacteroidales[["taxonomy"]] <- c("Bacteria","Bacteroidota","Bacteroidia","Bacteroidales","Porphyromonadaceae")
  Bacteroidales[["confidence"]] <- c(1.00, 1.00, 1.00, 1.00,0.82)

  expected <- "Bacteria(100);Bacteroidota(100);Bacteroidia(100);Bacteroidales(100);Porphyromonadaceae(82);Porphyromonadaceae_unclassified(82)"
  tax_string <- print_taxonomy(Bacteroidales,n_levels = 6)
  expect_equal(tax_string, expected)

})
