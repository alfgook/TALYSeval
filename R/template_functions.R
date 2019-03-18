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


#' Read TALYS Input Template
#'
#' @param filename path to TALYS input file
#'
#' @return character vector with lines in TALYS input file
#' @export
#'
readTalysTemplate <- function(filename) {

  scan(filename, character(0), sep="\n", quiet=TRUE)
}


#' Remove Lines From TALYS Template
#'
#' @param regPats character vector with regular expressions
#' @param template character vector with a TALYS template
#'
#' @return Modified template, i.e. with lines removed that matched some patterin in \code{regPats}
#' @export
#'
removeLinesFromTalysTemplate <- function(regPats, template) {

  for (curPat in regPats)
    template <- template[!grepl(curPat, template, ignore.case=TRUE)]
  template
}


#' Update TALYS Template
#'
#' @param input list with TALYS parameters
#' @param template character vector with a TALYS template
#'
#' @return A modified TALYS template including the parameter specifications in \code{input}
#' @export
#'
updateTalysTemplate <- function(input, template) {

  parNames <- names(input)
  stopifnot(!anyDuplicated(parNames))
  parNames <- ifelse(grepl("_", parNames),
                     parNames,
                     paste0(parNames, "_"))
  parNames <- gsub("_", " _ ", parNames)
  parNames <- gsub(" +", " ", parNames)
  parNames <- gsub("(^ *| *$)", "", parNames)
  escParNames <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1",
                      parNames)
  regPats <- gsub(" ", " +", escParNames)
  regPats <- gsub("_", ")([^ ]+)(", regPats)
  regPats <- paste0("^( *", regPats, " *)$")

  templateAnnex <- character(0)
  for (i in seq_along(input)) {

    curVal <- input[[i]]
    idcs <- grep(regPats[i], template, ignore.case=TRUE)
    if (length(idcs)>0) {
      template[idcs] <- gsub(regPats[i],
                             paste0("\\1",as.character(curVal),"\\3"),
                             template[idcs])
    } else {
      curParName <- tolower(parNames[i])
      curLine <- gsub("_", as.character(curVal), curParName)
      templateAnnex <- c(templateAnnex, curLine)
    }
  }
  c(template, templateAnnex)
}


#' Extract Adjusted TALYS Parameters
#'
#' Extract parameters of the form *adjust from a TALYS input file.
#'
#' @param template character vector with a TALYS template
#' @param last in the case some parameters are defined multiple times with last=TRUE
#'             the last specifications are returned and with last=FALSE the first one.
#'
#' @return a named list with the parameters and their values
#' @export
#'
extractAdjustedTalysParameters <- function(template, last=TRUE) {

  parList <- list()
  parNames <- character(0)
  lines <- grep("^( |[^#])*[0-9a-zA-Z]+adjust", template, value=TRUE, ignore.case=TRUE)
  lines <- gsub("(^ *| *$)", "", lines)
  splitRes <- strsplit(lines, "( |\\t)+")
  for (i in seq_along(lines)) {
    curSplitRes <- splitRes[[i]]
    idx <- grep("[0-9]+\\.[0-9]+", curSplitRes)
    stopifnot(length(idx)==1)
    curParName <- paste(curSplitRes[1:(idx-1)], collapse=" ")
    if (idx < length(curSplitRes))
      curParName <- paste0(curParName, "_",
                           paste(curSplitRes[(idx+1):length(curSplitRes)], collapse=" "))
    parNames[i] <- curParName
    parList[[i]] <- as.numeric(curSplitRes[idx])
  }
  parSel <- !duplicated(parNames, fromLast=last)
  parList <- parList[parSel]
  parNames <- parNames[parSel]
  names(parList) <- tolower(parNames)
  parList
}
