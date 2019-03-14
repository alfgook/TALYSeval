
#' Normalize TALYS Input
#'
#' @param input list with TALYS parameters
#'
#' @return normalized list with TALYS parameters
#' @export
#'
normalizeTalysInput <- function(input) {

  names(inp) <- tolower(names(inp))
  lapply(inp, function(x) { if (is.character(x)) tolower(x) else x })
}



#' Define an energy grid
#'
#' @param energies a sorted experimental energy grid
#' @param compGrid a sorted computational grid
#' @param enPolicy see result section.
#'
#' @return
#' The returned energy grid depends on \code{enPolicy}
#' \tabular{ll}{
#'   \code{expgrid} \tab returns the experimental grid with duplicates removed \cr
#'   \code{compgrid} \tab returns the computational grid \cr
#'   \code{thinning} \tab returns a subset of the computational grid which encloses the energies in the experimental grid
#' }
#' @export
#'
defineEnergyGrid <- function(energies, compGrid=NULL, enPolicy="thinning") {

  if (is.unsorted(compGrid))
    stop("the provided energy grid must be sorted")

  if (!enPolicy %in% c("thinning","expgrid","compgrid"))
    stop("enPolicy must be 'thinning', 'expgrid' or 'compgrid'")

  if (enPolicy %in% c("thinning", "compgrid") && is.null(compGrid))
    stop("computational grid must be specified")

  if (enPolicy=="thinning")
  {
    enIdx <- findInterval(energies, compGrid)
    liesOnGrid <- energies == compGrid[enIdx]
    enIdxNotGrid <- enIdx[!liesOnGrid]

    if (any(enIdx==0) || any(enIdxNotGrid==length(compGrid)))
      stop("Some experimental energy outside computational grid")

    enIdx <- sort(unique(c(enIdxNotGrid, enIdxNotGrid+1,
                           enIdx[liesOnGrid])))
    enGrid <- compGrid[enIdx]
  }
  else if (enPolicy=="compgrid")
  {
    stopifnot(min(energies) >= min(compGrid) && max(energies) <= max(compGrid))
    enGrid <- compGrid
  }
  else if (enPolicy=="expgrid")
  {
    enGrid <- sort(unique(energies))
  }
  else
  {
    stop("should not happen by construction!")
  }

  enGrid
}


approxWithSnap <- function(x, y = NULL, xout, ..., digits=6) {

  if (length(x) == 1 && length(y)==1) {
    if (!all(signif(xout, digits) == signif(x, digits)))
      stop("only one x-value provided, all xout must be identical to x")
    list(x = rep(x, length(xout)),
         y = rep(y, length(y)))
  }
  else
    approx(signif(x, digits), y, signif(xout, digits), ...)
}





