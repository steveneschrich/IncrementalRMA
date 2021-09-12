
#' Extract value or return default
#'
#' A pretty clever approach for accessing part of a list,
#'  returning a default if there is nothing of that name in the list.
#' This was taken from https://stackoverflow.com/questions/26561518/built-in-for-getting-list-element-or-default
#'
#'
#' @param a Value to use, unless empty.
#' @param b Alternative value to use, if a is empty.
#'
#' @return a unless empty, in which case b is returned.
#' @export
#'
#' @examples
#'
#' NULL %or% "default"
#'
`%or%` <- function(a, b) {
  cmp = function(a,b) if (identical(a, FALSE) ||
                          is.null(a) ||
                          is.na(a) ||
                          is.nan(a) ||
                          length(a) == 0) b else a

  if (length(a) > 1)
    mapply(cmp, a, b)
  else
    cmp(a, b)
}
