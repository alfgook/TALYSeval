#  TALYSeval - working with TALYS from R
#  Copyright (C) 2019  Georg Schnabel
#  
#  TALYSeval is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  TALYSeval is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>

#' Convert Extended Output Specification
#'
#' @param extSpec extended output specification
#' @param compGrid incident energy grid for the computation
#' @param enPolicy1 policy to reduce energy grid passed to TALYS
#' @param enPolicy2 policy to construct energy grid to which TALYS results will be interpolated
#'
#' @return A list with elements \code{inputs} \code{outspecs}
#'         \code{ACCNUM, SUBACCNUM, REAC}
#'
#' @note TODO: energy grid reduction for DA/DE data
#'
#' @export
#' @import data.table
#'
convertToInput <- function(paramDt, needsDt) {

  reqCols <- c("PROJECTILE", "ELEMENT", "MASS", "PARNAME", "PARVAL")
  if (!all(reqCols %in% names(paramDt)))
    stop(paste0("paramDt must contain columns ", paste0(reqCols, collapse=", ")))
  forbiddenParnames <- c("PROJECTILE", "MASS", "ELEMENT")
  if (any(toupper(paramDt$PARNAME) %in% forbiddenParnames))
    stop(paste0("paramDt$PARNAME must not contain ", paste0(forbiddenParnames, collapse=", ")))

  setkey(paramDt, PROJECTILE, ELEMENT, MASS)
  needsDt[,{
    curInp <- copy(.BY)
    names(curInp) <- tolower(names(curInp))
    curInpRows <- paramDt[.BY,,nomatch=0]
    if (nrow(curInpRows) > 0)
    {
      if (! "energy" %in% curInpRows$PARNAME)
        stop("energy must be present in each parameter list")
      additionalParams <- curInpRows$PARVAL
      names(additionalParams) <- tolower(curInpRows$PARNAME)
      curInp <- modifyList(curInp, additionalParams)
    }
    else
      stop("each parameter set needs at least the energy keyword")

    list(inputs=list(curInp), outspecs=.(copy(.SD)))
  }, by=c("PROJECTILE", "ELEMENT", "MASS")]
}


#' Convert Calculation Results
#'
#' @param results list with calculation results
#' @param errorPolicy default 'skip', skips incomplete results, 'abort' aborts
#'
#' @return An extended output specification
#' @export
#'
convertResultsToNeeds <- function(results, errorPolicy="skip") {

  stopifnot(is.list(results))
  rbindlist(lapply(results, function(el) {

    if (!is.null(el$error)) {
      if (!isTRUE(errorPolicy=="skip"))
        stop("incomplete result found")
      else
        NULL
    } else {
      if (is.null(el$result))
        stop("element 'result' does not exist in one list element, should not happen!")

      data.table(PROJECTILE = toupper(el$input$projectile),
                 ELEMENT = toupper(el$input$element),
                 MASS = el$input$mass,
                 el$result)
    }
  }))
}


#' Compactify Extended Output Specification
#'
#' @param needsDt datatable with columns \code{PROJECTILE,ELEMENT,MASS,REAC,L1,L2,L3}
#' @param compGrid computational grid for the incident energies
#' @param enPolicy how to reduce the energies, choices are 'thinning', 'expgrid', and 'compgrid'
#'
#' @return datatable of the same form as \code{needsDt} but potentially less rows
#' @export
#'
compactifyNeeds <- function(needsDt, compGrid, enPolicy = "thinning") {

  resDt <- needsDt[,{

    if (grepl("^CS/", .BY[["REAC"]])) {
      enGrid <- defineEnergyGrid(L1, compGrid, enPolicy)
      list(L1=enGrid, L2=0, L3=0)
    }
    else
      list(L1=L1, L2=L2, L3=L3)

  }, by=c("PROJECTILE","ELEMENT","MASS", "REAC")]
  resDt[, IDX:=seq_len(.N)]
  resDt[]
}



#' Add Energies To ParamDt
#'
#' @param paramDt a datatable with parameter specifications
#' @param needsDt a datatable with output specifications
#' @param compGrid the energy grid used for computation
#' @param enPolicy how to reduce the energy, i.e. \code{'expgrid', 'compgrid', 'thinning'}
#'
#' @return a paramDt with the energy keyword added to the parameter lists
#' @export
#'
addEnergiesToParamDt <- function(paramDt, needsDt, compGrid, enPolicy="thinning") {

  reqCols <- c("PROJECTILE", "ELEMENT", "MASS", "PARNAME", "PARVAL", "ADJUSTABLE")
  if (!all(reqCols %in% names(paramDt)))
    stop(paste0("paramDt must contain columns ", paste0(reqCols, collapse=", ")))
  forbiddenParnames <- c("PROJECTILE", "MASS", "ELEMENT")
  if (any(toupper(paramDt$PARNAME) %in% forbiddenParnames))
    stop(paste0("paramDt$PARNAME must not contain ", paste0(forbiddenParnames, collapse=", ")))

  setkeyv(needsDt, c("PROJECTILE", "ELEMENT", "MASS", "REAC", "L1", "L2", "L3"))
  paramDt[,{

    origEnGrid <- needsDt[.BY, L1]
    enGrid <- defineEnergyGrid(origEnGrid, compGrid, enPolicy)
    enIdx <- which("energy" == .SD$PARNAME)
    if (length(enIdx) > 1)
      stop("parameter name 'energy' must not occur more than once")
    else if (length(enIdx) == 1) {
      PARVAL[[enIdx]] <- enGrid
    } else {
      PARNAME <- c(PARNAME, "energy")
      PARVAL <- c(PARVAL, list(enGrid))
      ADJUSTABLE <- c(ADJUSTABLE, FALSE)
    }
    list(PARNAME = PARNAME, PARVAL = PARVAL, ADJUSTABLE=ADJUSTABLE)
  }, by=c("PROJECTILE", "ELEMENT", "MASS")]
}



