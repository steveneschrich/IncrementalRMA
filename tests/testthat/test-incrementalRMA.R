test_that("Incremental RMA provides same answers as RMA.", {
    # Super-suppressed because of warnings, status messages, etc.
    suppressPackageStartupMessages({
      suppressWarnings({
        library(affy)
        library(affydata)
      })
    })
    data(Dilution, verbose = FALSE)
    suppressWarnings({
      incrma <- incrementalRMA(Dilution, parameterizeRMA(Dilution))
      canrma <- affy::rma(Dilution)
    })
    expect_equal(
      exprs(incrma),
      exprs(canrma)

  )
})
