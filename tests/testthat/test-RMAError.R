

test_that("Combining duplicate AffyBatch objects doesn't add new samples.", {
  # Super-suppressed because of warnings, status messages, etc.
  suppressPackageStartupMessages({
    suppressWarnings({
      library(affy)
      library(affydata)
    })
  })
  data(Dilution, verbose = FALSE)
  rep <- combine_and_deduplicate(Dilution, Dilution)
  expect_equal(sampleNames(rep), sampleNames(Dilution))
  expect_identical(Biobase::exprs(rep), Biobase::exprs(Dilution))
})

test_that("Error in incrementalRMA is zero when it is not incremental.", {
  # Super-suppressed because of warnings, status messages, etc.
  suppressPackageStartupMessages({
    suppressWarnings({
      library(affy)
      library(affydata)
    })
  })
  data(Dilution, verbose = FALSE)
  suppressWarnings({
    incrma <- incrementalRMA(Dilution, parameterizeRMA(Dilution), calculate_error = TRUE )
  })
  expect_equal(Biobase::assayData(incrma)[["se.exprs"]],
               matrix(0, ncol = ncol(incrma), nrow = nrow(incrma), dimnames=dimnames(incrma)),
               tolerance = 1e-10)
})

test_that("Error in incrementalRMA is non-zero when data is slightly different.", {
  # Super-suppressed because of warnings, status messages, etc.
  suppressPackageStartupMessages({
    suppressWarnings({
      library(affy)
      library(affydata)
    })
  })
  data(Dilution, verbose = FALSE)
  suppressWarnings({
    incrma <- incrementalRMA(Dilution[,4], parameterizeRMA(Dilution[,1:3]), calculate_error = TRUE )
  })
  expect_equal(quantile(Biobase::assayData(incrma)[["se.exprs"]], c(0, 0.25, 0.5, 0.75, 1)),
               c("0%"=-0.3789974, "25%"=0.1437935, "50%"=0.1844053, "75%"=0.2323840, "100%"=0.8848735),
               tolerance = 1e-6
  )
})
