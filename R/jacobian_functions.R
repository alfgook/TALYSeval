#' Create Default Mask
#'
#' @param paramDt datatable with parameter specifications
#' @param needsDt datatable with output specifications
#'
#' @return
#' datatable with columnns \code{SRCIDX, DSTIDX} telling which parameters associated with a \code{SRCIDX}
#' have an impact on which observables associated with \code{DSTIDX}
#' @export
#'
createDefaultMask <- function(paramDt, needsDt) {

  if (! "IDX" %in% names(paramDt) || ! "IDX" %in% names(needsDt))
    stop("both paramDt and needsDt must contain a column IDX")
  reqColsInNeedsDt <- c("PROJECTILE", "ELEMENT", "MASS", "REAC", "L1", "L2", "L3")
  if (! all(reqColsInNeedsDt %in% names(needsDt)))
    stop(paste0("needsDt must contain the columns ", paste0(reqColsInNeeds, collapse=", ")))
  reqColsInParamDt <- c("PROJECTILE", "ELEMENT", "MASS")
  if (! all(reqColsInParamDt %in% names(paramDt)))
    stop(paste0("paramDt must contain the columns ", paste0(reqColsInParamDt, collapse=", ")))


  setkey(needsDt, PROJECTILE, ELEMENT, MASS, REAC, L1, L2, L3)
  setkey(paramDt, PROJECTILE, ELEMENT, MASS)
  paramDt[,{
    srcIdx <- IDX[ADJUSTABLE]
    dstIdx <- needsDt[.BY, IDX]
    stopifnot(!any(is.na(dstIdx)))
    expand.grid(SRCIDX = srcIdx, DSTIDX = dstIdx)
  },by=c("PROJECTILE", "ELEMENT", "MASS")][, list(SRCIDX, DSTIDX)]
}


#' Create Inputs To Calculate Jacobian
#'
#' @param paramDt datatable with parameter specifications
#' @param needsDt datatable with output specifications
#' @param mask datatable which determines the dependence structure between parameters and observables
#' @param eps difference to employ in the finite difference approximation to the gradient
#' @param trafo either NULL or a list containing the functions \code{fun} and \code{invfun}
#'
#' @return
#' datatable with list columns \code{inputs} and \code{outspecs}
#' @export
#'
createInputsForJacobian <- function(paramDt, needsDt, mask=NULL, eps=1e-3, trafo = NULL) {

  if (!is.null(trafo)) {
    stopifnot(is.function(trafo$fun))
    stopifnot(is.function(trafo$invfun))
  }

  if (! "IDX" %in% names(paramDt) || ! "IDX" %in% names(needsDt))
    stop("both paramDt and needsDt must contain a column IDX")

  if (is.null(mask))
    mask <- createDefaultMask(paramDt, needsDt)

  setkey(paramDt, IDX)
  setkey(needsDt, IDX)
  setkey(mask, SRCIDX, DSTIDX)

  resDt <- paramDt[, {

    adjustableIdcs <- mask[J(IDX[ADJUSTABLE]), list(SRCIDX, DSTIDX), mult="first", nomatch=NA]
    existsIdx <- which(!is.na(adjustableIdcs$DSTIDX))
    adjustableIdcs <- adjustableIdcs[existsIdx, SRCIDX]
    subAdjIdcs <- which(ADJUSTABLE)[existsIdx]
    numVars <- length(adjustableIdcs)

    curList <- list(PARNAME = NA_character_,
                    ADJIDX = c(0, adjustableIdcs),
                    EPS = eps,
                    inputs = replicate(numVars+1, NULL, simplify = FALSE),
                    outspecs = replicate(numVars+1, NULL, simplify = FALSE))

    curNeedsDt <- copy(needsDt)
    refInp <- PARVAL
    names(refInp) <- tolower(PARNAME)

    sysInp <- copy(.BY)
    names(sysInp) <- tolower(names(sysInp))

    refInp <- c(sysInp, refInp)
    refInp <- refInp[!duplicated(names(refInp))]

    curList$PARNAME[1] <- NA_character_
    curList$ADJIDX[1] <- 0
    curList$inputs[[1]] <- refInp
    curList$outspecs[[1]] <- curNeedsDt

    # create variations
    for (j in seq_along(adjustableIdcs)) {

      # print(paste0("j ", j, " out of ", length(adjustableIdcs)))
      subIdx <- subAdjIdcs[j]
      curIdx <- adjustableIdcs[j]
      curParname <- tolower(PARNAME[subIdx])

      uniqueDstIdcs <- mask[J(curIdx), DSTIDX]
      curNeedsDt <- needsDt[J(uniqueDstIdcs),]
      stopifnot(!any(is.na(curNeedsDt$L1)))

      origVal <- PARVAL[[subIdx]]
      variation <- {
        if (is.null(trafo))
          origVal + eps
        else
          trafo$fun(trafo$invfun(origVal) + eps)
      }

      refInp[[curParname]] <- variation
      refInp$energy <- sort(unique(curNeedsDt$L1))

      curList$PARNAME[j+1] <- curParname
      curList$ADJIDX[j+1] <- curIdx
      curList$inputs[[j+1]] <- refInp
      curList$outspecs[[j+1]] <- curNeedsDt

      refInp[[curParname]] <- origVal
    }

    curList
  }, by=c("PROJECTILE", "ELEMENT", "MASS")]

  resDt[]
}


#' Compute Jacobian
#'
#' @param variationDt a datatable of the same structure as returned by \link{createInputsForJacobian} and with
#'                    \code{outspecs} with an additional column \code{V1} containg the predictions of TALYS.
#' @param drop0 logical. Should rows corresponding to zero elements be removed?
#'
#' @return
#' datatable with columns \code{IDX1, IDX2, X} indicating values and positions of non-zero elements of the
#' Jacobian matrix.
#' @export
#'
computeJacobian <- function(variationDt, drop0 = FALSE) {

  setkey(variationDt, ADJIDX)
  curRef <- NULL
  resultDt <- variationDt[,{

    if (ADJIDX==0) {
      curRef <<- outspecs[[1]]
      setkey(curRef, IDX)
      NULL
    }
    else
    {
      curOut <- outspecs[[1]]
      sensElements <- (curOut$V1 - curRef[J(curOut$IDX),V1]) / EPS
      list(IDX1 = curOut$IDX, IDX2 = ADJIDX, X = sensElements)
    }
  }, by="ADJIDX"]
  resultDt[, IDX1 := as.integer(IDX1)]
  resultDt[, IDX2 := as.integer(IDX2)]
  resultDt[, ADJIDX := NULL]
  if (drop0) resultDt <- resultDt[X != 0,]
  resultDt[]
}

