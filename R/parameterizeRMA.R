#' Calculate parameters for incremental RMA
#'
#' @description Calculate and store parameters for future use in incrementalRMA.
#'
#' @details
#' IncrementalRMA works by applying previous computed parameters from a full RMA process. This
#' function calculates the necessary parameters for future \code{\link{incrementalRMA}} use.
#'
#' RMA consists of three steps: background correction, quantile normalization and median polish.
#' Background correction is done per-chip, so there is no need for pre-calculating parameters.
#' Quantile normalization computes values for each quantile in expression data and then sets the
#' intensity distribution equal to these quantiles. Median polish sweeps out row-wise and column-wise
#' medians from probe expression; column-wise effects are used as expression estimates.
#'
#' Given the above description, quantiles (from quantile normalization) and row-wise effects
#' (from median polish) can be computed on a set of samples and stored. With this information,
#' \code{\link{incrementalRMA}} can then calculate gene expression from a new sample without
#' re-computing these parameters.
#'
#' This function computes these parameters and stores them in a list, such that the list can then
#' be used. The elements of the list currently include:
#' \itemize{
#'   \item{\code{probeEffects}}{Row-wise effects from median polish step.}
#'   \item{\code{normalizationVector}}{Quantile values from quantile normalization step.}
#'   \item{\code{referenceCELFiles}}{An \code{\link[affy]{AffyBatch-class}} containing the reference cel
#'   files (for error calculations).}
#' }
#' A simple test of incrementalRMA is to parameterizeRMA, then apply incrementalRMA and compare with
#' using rma on the set. See the example for code to test this.
#'
#' @param abatch An \code{\link[affy]{AffyBatch-class}} to calculate parameters on.
#'
#' @return A list containing
#' @export
#'
#' @examples
#' \dontrun{
#' params <- parameterizeRMA(abatch)
#'
#' all.equal(exprs(incrementalRMA(abatch, parameterizeRMA(abatch))), exprs(affy::rma(abatch)))
#' }
parameterizeRMA <- function(abatch) {
  stopifnot("AffyBatch is required for input." = class(abatch)=="AffyBatch")

  params <- list(referenceCELFiles = abatch)

  ## Note: The trick with doing this work is that we have to essentially calculate RMA on
  ## the reference set, just storing things along the way.

  ## RMA Step 1: Background correction - does not require parameterization.
  message("Background Correcting ...")
  abatch <- affy::bg.correct.rma(abatch)

  ## RMA Step 2: Quantile normalization - requires parameterization
  message("Normalizing ...")
  params$normalizationVector <- preprocessCore::normalize.quantiles.determine.target(affy::pm(abatch))
  names(params$normalizationVector) <- affy::probeNames(abatch)

  pms <- preprocessCore::normalize.quantiles.use.target(affy::pm(abatch), params$normalizationVector)

  ## RMA Step 3: Summarization with median polish - requires parameterization
  message("Summarizing ...")
  mp <- preprocessCore::subrcModelMedianPolish(log2(pms), affy::probeNames(abatch))
  # Take out estimates: the first ncol(pms) are column effects so ignore these.
  estimates <- sapply(mp, function(tmp){tmp$Estimates[(ncol(pms)+1):length(tmp$Estimates)]})
  # Set the names next
  params$probeEffects<-stats::setNames(unlist(estimates),
                rep(names(estimates), sapply(estimates, length)))

  params
}
