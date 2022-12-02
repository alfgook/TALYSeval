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

#' Create a TALYS model class
#'
#' \code{createModelTALYS} returns a list which contains functions to
#' \itemize{
#'   \item{setup TALYS model calculations}
#'   \item{provide information to run the code}
#'   \item{read the results}
#'   }
#'
#' @param initList not used at the moment. Future: a list of additional configuration variables.
#'
#' @return Returns a list which provides the following functions:
#'
#' \item{\code{getName()}}{returns a model identification string.}
#' \item{\code{getArgs()}}{returns a string with the command line arguments to be used in call of the executable.}
#' \item{\code{prepare(dir, parset)}}{sets up the files in directory \code{dir}
#'       required to start a model calculation. Parameters for the calculation are passed in the list \code{parset}.}
#' \item{\code{finalize(dir)}}{has to be called after the completion of the model calculation.
#'       This function can do some postprocessing of the results. Returns \code{TRUE} if calculation finished successfully,
#'       otherwise \code{FALSE}.}
#' \item{\code{read(dir, tbl)}}{reads the results of a model calculation performed in directory \code{dir}.
#'       This function assumes that \code{finalize} was called before.}
#' @export
#' @import data.table
createModelTALYS <- function(initList = NULL) {

  model <- list()

  # parameters (only internally used)

  model[["initList"]] <- initList

  # functions (called from the outside)

  model[["getName"]] <- function() {

    return("TALYS")
  }

  model[["getArgs"]] <- function() {

    paste0("< ",model$intern$inputFile," > ",model$intern$outputFile)
  }

  model[["getScript"]] <- function() {

    file.path(path.package("TALYSeval"),"exec/talys.sh")
  }


  model[["prepare"]] <- function(dir, parset, clusterOpts=NULL) {

    model$intern$prepareCalc(dir,parset)
  }

  model[["finalize"]] <- function(dir, ..., clusterOpts=NULL) {

    # some additional information
    extraInfo <- Sys.time()
    extraFile <- file.path(dir,"extrainfo")
    write(extraInfo, extraFile)
    # check if successful
    outputPath <- file.path(dir,model$intern$outputFile)
    if (file.exists(outputPath) &&
        grepl('TALYS team congratulates',system(paste0("tail -n1 ",outputPath),intern=TRUE),fixed=TRUE))
      return(TRUE)
    else
      return(FALSE)
  }

  model[["read"]] <- function(dir, ...) {

    return(model$intern$read(dir, ...))
  }

  # functions only internally used
  model$intern$inputFile <- "input"
  model$intern$outputFile <- "output"

  model$intern$prepareCalc <- function(dir,parset) {

    inputFile <- model$intern$inputFile
    absInputFile <- file.path(dir,inputFile)
    # remove the template from the parameter list
    template <- parset$template
    parset$template <- NULL

    # special treatment for some input keywords

    # special treatment of energy dependent parameters
    endepParPat <- "^ *([0-9A-Za-z]+adjust)\\(([^)]+)\\) *([npdtha]) *$"
    endepParSel <- grepl(endepParPat, names(parset))
    endepPars <- names(parset)[endepParSel]
    if (length(endepPars) > 0) {
      endepParLines <- regmatches(endepPars, regexec(endepParPat, endepPars))
      dt <- as.data.table(t(matrix(unlist(endepParLines), nrow=4)))
      setnames(dt, c("line", "parname", "energy", "projectile"))
      dt[, energy := as.numeric(energy)]
      dt[, multfact := as.numeric(unlist(parset[endepParSel]))]
      dt[, lineidx := which(endepParSel)]
      setkey(dt, parname, projectile, energy)
      dt[, {
        curFilename <- paste0(parname, "-", projectile, ".table")
        write.table(.SD[, list(sprintf("%.8f", energy),
                               sprintf("%.8f", multfact))],
                    file = file.path(dir, curFilename), row.names = FALSE,
                    col.names = FALSE, quote = FALSE)
        names(parset)[lineidx[1]] <<- paste0(.BY[["parname"]], " ", .BY[["projectile"]], " 1.")
        parset[[lineidx[1]]] <<- curFilename
        NULL
      }, by=c("parname", "projectile")]
      
      deleteLines <- dt[, list(lineidx = lineidx[-1]), by=c("parname", "projectile")]
      if (length(deleteLines$lineidx) <= 0)
          stop("at least two energy points for each energy-dependent parameter needed") 
      parset <- parset[-deleteLines$lineidx]
    }

    # special treatment of energy
    stopifnot(sum("energy" == names(parset))==1)
    energyValues <- parset[["energy"]]
    parset[["energy"]] <- "energies"
    write(signif(energyValues, digits = 6),  # round to precision of TALYS
          file.path(dir,"energies"),ncolumns=1)

    # sepcial treatment of ompenergy
    if ("ompenergy" %in% names(parset)) {
      stopifnot(sum(names(parset) %in% "ompenergy")==1)
      ompenergyValues <- parset[["ompenergy"]]
      names(parset)[names(parset)=="ompenergy"] <- "ompenergyfile"
      parset[["ompenergyfile"]] <- "ompenergies"
      write(signif(ompenergyValues, digits = 6), # round to precision of TALYS
            file.path(dir, "ompenergies"), ncolumns=1)
    }

    # generic treatment for the rest
    parlines <- updateTalysTemplate(parset, template)

    # write the input file and the energy file
    write(parlines,absInputFile)
  }

  model$intern$read <- function(calcdir,spec,keep.temp=FALSE, packed=TRUE) {

    if (!length(calcdir)==1 || !is.character(calcdir))
      stop("calcdir must be a character vector of length one")
    if (!dir.exists(calcdir))
      stop(paste0("directory ", calcdir, " does not exist"))
    if (!is.data.table(spec))
      stop("spec must be a data.table")
    if (!all(c("REAC","L1","L2","L3") %in% colnames(spec)))
      stop("spec must contain columns REAC, L1, L2, L3")

    originalCalcdir <- calcdir
    if (!isTRUE(packed))
      unpackedCalcdir <- calcdir
    else
      unpackedCalcdir <- file.path(calcdir,"unpacked")

    # prepare data.table
    setkey(spec,REAC,L1,L2,L3)
    if ("V1" %in% names(spec))
      spec[,V1:=NULL]

    reacFun <- function(calcdir, reacstr, data, retfiles) {
      if (grepl("^CS/TOT$",reacstr))
        model$intern$getTotalXS(calcdir,reacstr,data,retfiles)
      else if (grepl("^CS/EL$",reacstr))
        model$intern$getElasticXS(calcdir,reacstr,data,retfiles)
      else if (grepl("^CS/NONEL$",reacstr))
        model$intern$getNonelasticXS(calcdir,reacstr,data,retfiles)
      else if (grepl("^CS/REAC/[[:digit:]]{6}/TOT$",reacstr))
        model$intern$getReacXS(calcdir,reacstr,data,retfiles)
      else if (grepl("^CS/REAC/([[:digit:]]{6})/L([[:digit:]]+)$",reacstr))
        model$intern$getReacLev(calcdir,reacstr,data,retfiles)
      else if (grepl("^CS/LEV/(n[npdtha])/L([[:digit:]]+)$",reacstr))
        model$intern$getLevelXS(calcdir,reacstr,data,retfiles)
      else if (grepl("^SP/PART/([npdtha])$",reacstr))
        model$intern$getPartSpec(calcdir,reacstr,data,retfiles)
      else if (grepl("^CS/BINARY/([npdtha])/(L[0-9][0-9]|CON)$",reacstr))
        model$intern$getBinaryXS(calcdir,reacstr,data,retfiles)
      else if (grepl("^CS/PROD/[[:digit:]]{6}/TOT$",reacstr))
        model$intern$getProdXS(calcdir,reacstr,data,retfiles)
      else
        stop(paste0("cannot handle REAC = ", reacstr))
    }

    files <- spec[,list(files=unique(reacFun(originalCalcdir,.BY[[1]],.SD,retfiles=TRUE))),by=REAC][,files]
    files <- files[!file.exists(file.path(unpackedCalcdir,files))]

    # make the files accessible in a temporary directory
    if (length(files) > 0 && isTRUE(packed)) {

      dir.create(unpackedCalcdir,recursive=FALSE,showWarnings=FALSE)
      cmdStr <- paste("tar","-xf",file.path(originalCalcdir,"archive.tgz"),
                      "-C",unpackedCalcdir,sep=" ")
      print(cmdStr)
      cmdStr <- paste(cmdStr,paste0(files,collapse=" "),sep=" ")
      system(cmdStr)
    }

    result <- spec[, V1:=reacFun(unpackedCalcdir, .BY[[1]],
                                 .SD, retfiles=FALSE), by=REAC]


    # remove the unpacked folder
    if (!isTRUE(keep.temp) && isTRUE(packed))
      unlink(unpackedCalcdir,recursive=TRUE)

    result[]
  }

  # SPECIALIZED HANDLERS
  # each handler retrieves a certain type of cross section
  # linear/bilinear interpolation is performed if required

  # specialized handler to get total xs
  model$intern$getTotalXS <- function(calcdir,reacstr,spec,retfiles=FALSE) {

    xsfile <- "totalxs.tot"
    if (isTRUE(retfiles))
      return(xsfile)

    # provide full path
    xsfile <- file.path(calcdir,xsfile)

    xsdata <- model$intern$readfixed(xsfile,c(1,13),c(12,24))
    approxWithSnap(xsdata[,1],xsdata[,2],xout=spec[,L1])$y
  }

  # specialized handler to get elastic xs
  model$intern$getElasticXS <- function(calcdir,reacstr,spec,retfiles=FALSE) {

    xsfile <- "elastic.tot"
    if (isTRUE(retfiles))
      return(xsfile)

    # provide full path
    xsfile <- file.path(calcdir,xsfile)

    xsdata <- model$intern$readfixed(xsfile,c(1,13),c(12,24))
    approxWithSnap(xsdata[,1],xsdata[,2],xout=spec[,L1])$y
  }

  # specialized handler to get nonelastic xs
  model$intern$getNonelasticXS <- function(calcdir,reacstr,spec,retfiles=FALSE) {

    xsfile <- "nonelastic.tot"
    if (isTRUE(retfiles))
      return(xsfile)

    # provide full path
    xsfile <- file.path(calcdir,xsfile)

    xsdata <- model$intern$readfixed(xsfile,c(1,13),c(12,24))
    approxWithSnap(xsdata[,1],xsdata[,2],xout=spec[,L1])$y
  }


  # specialized handler to read reaction xs
  model$intern$getReacXS <- function(calcdir,reacstr,spec,retfiles=FALSE) {

    particleCode <- regmatches(reacstr,regexec("^CS/REAC/([[:digit:]]{6})/TOT$",reacstr))[[1]][2]
    stopifnot(!is.na(particleCode))

    # construct filename
    xsfile <- paste0("xs",particleCode,".tot")
    if (isTRUE(retfiles))
      return(xsfile)
    xsfile <- file.path(calcdir,xsfile)
    # retrieve data
    xsdata <- model$intern$readfixed(xsfile,c(1,13),c(12,24))
    approxWithSnap(xsdata[,1],xsdata[,2],xout=spec[,L1])$y
  }


  # specialized handler to read binary xs
  model$intern$getBinaryXS <- function(calcdir, reacstr, spec, retfiles = FALSE) {

    regPat <- "^CS/BINARY/([npdtha])/(L[0-9][0-9]|CON|TOT)$"
    regRes <- regexec(regPat, reacstr)
    particleCode <- regmatches(reacstr, regRes)[[1]][2]
    levelCode <- regmatches(reacstr, regRes)[[1]][3]

    if (levelCode %in% c("TOT","CON")) levelCode <- tolower(levelCode)
    stopifnot(!is.na(particleCode))
    stopifnot(!is.na(levelCode))

    # construct filename
    xsfile <- paste0("n",particleCode,".",levelCode)
    print(xsfile)
    if (isTRUE(retfiles))
      return(xsfile)
    xsfile <- file.path(calcdir,xsfile)
    # retrieve data
    xsdata <- model$intern$readfixed(xsfile,c(1,13),c(12,24))
    approxWithSnap(xsdata[,1],xsdata[,2],xout=spec[,L1])$y
  }


  # specialized handler to read reaction xs scattered at some level
  model$intern$getReacLev <- function(calcdir,reacstr,spec,retfiles=FALSE) {

    particleCode <- regmatches(reacstr,regexec("^CS/REAC/([[:digit:]]{6})/L([[:digit:]]+)$",reacstr))[[1]][-1]
    stopifnot(!any(is.na(particleCode)))

    # construct filename
    xsfile <- sprintf("xs%s.L%02d",particleCode[1],as.integer(particleCode[2]))
    if (isTRUE(retfiles))
      return(xsfile)
    xsfile <- file.path(calcdir,xsfile)

    # retrieve data
    xsdata <- model$intern$readfixed(xsfile,c(1,13),c(12,24))
    approxWithSnap(xsdata[,1],xsdata[,2],xout=spec[,L1])$y
  }




  # specialized handler to read scattering at levels
  model$intern$getLevelXS <- function(calcdir,reacstr,spec,retfiles=FALSE) {

    particleCode <- regmatches(reacstr,regexec("^CS/LEV/(n[npdtha])/L([[:digit:]]+)$",
                                               reacstr))[[1]][2:3]
    stopifnot(!any(is.na(particleCode)))

    # construct filename
    xsfile <- sprintf("%s.L%02d",particleCode[1],as.integer(particleCode[2]))
    if (isTRUE(retfiles))
      return(xsfile)
    xsfile <- file.path(calcdir,xsfile)

    # retrieve data
    xsdata <- model$intern$readfixed(xsfile,c(1,11),c(10,22))
    approxWithSnap(xsdata[,1],xsdata[,2],xout=spec[,L1])$y
  }

  # specialized handler to read spectra
  model$intern$getPartSpec <- function(calcdir,reacstr,spec,retfiles=FALSE) {

    particleCode <- regmatches(reacstr,regexec("^SP/PART/([npdtha])$",
                                               reacstr))[[1]][2]
    stopifnot(!any(is.na(particleCode)))

    # get relevant filenames
    hypxsfiles <- sprintf("%sspec%07.3f.tot",particleCode,spec[,L1])
    filelistFile <- file.path(sub("/unpacked$","",calcdir),"filelist")
    filelist <- fread(filelistFile,header=FALSE,col.names="file")
    setattr(filelist, "sorted", "file")
    idx <- filelist[J(hypxsfiles), list(idx=.I), roll = TRUE, by = .EACHI][,idx]
    regex <- paste0("^",particleCode,"spec([[:digit:]]{3}\\.[[:digit:]]{3})")
    idxLo <- ifelse(idx>0 & grepl(regex,filelist[idx,file]),idx,NA)
    idxHi <- ifelse(idx<nrow(filelist) & grepl(regex,filelist[idx+1,file]),idx+1,NA)

    if (isTRUE(retfiles)) {
      validIdx <- c(idxLo[!is.na(idxLo)],idxHi[!is.na(idxHi)])
      return(unique(filelist[validIdx,file]))
    }

    # retrieve the data
    xsfilesLo <- filelist[idxLo,file]
    xsfilesHi <- filelist[idxHi,file]
    foundLo <- !is.na(xsfilesLo)
    foundHi <- !is.na(xsfilesHi)
    EincLo <- rep(NA,length(xsfilesLo))
    EincHi <- rep(NA,length(xsfilesHi))
    EincLo[foundLo] <- as.numeric(unlist(regmatches(xsfilesLo[foundLo],regexec(regex,xsfilesLo[foundLo])))[c(FALSE,TRUE)])
    EincHi[foundHi] <- as.numeric(unlist(regmatches(xsfilesHi[foundHi],regexec(regex,xsfilesLo[foundHi])))[c(FALSE,TRUE)])
    xsfilesLo[foundLo] <- file.path(calcdir,xsfilesLo[foundLo])
    xsfilesHi[foundHi] <- file.path(calcdir,xsfilesHi[foundHi])
    Ediff <- EincHi - EincLo

    # prepare a data table for manipulation
    tmpspec <- copy(spec)
    tmpspec[,FILEIDX:=idx]
    tmpspec[,c("EincLo","EincHi","Ediff"):=list(EincLo,EincHi,Ediff)]
    tmpspec[,c("xsfilesLo","xsfilesHi"):=list(xsfilesLo,xsfilesHi)]
    tmpspec[,V1:={
      if (!is.na(xsfilesLo[1]) && !is.na(xsfilesHi[1])) {

        tbl1 <- model$intern$readfixed(xsfilesLo[1],c(1,8),c(7,19))
        tbl2 <- model$intern$readfixed(xsfilesHi[1],c(1,8),c(7,19))

        yLo <- approxWithSnap(tbl1[,1],tbl1[,2],xout=L2)$y
        yHi <- approxWithSnap(tbl2[,1],tbl2[,2],xout=L2)$y
        res <- (EincHi-L1)/Ediff * yLo + (L1-EincLo)/Ediff* yHi
        # assume missing values (interpolation outside the table) to be zero
        res[is.na(res)] <- 0
        res
      }
      else {
        0
      }
    },by=FILEIDX]

    return(tmpspec[,V1])
  }


  # LOW LEVEL CONVENIENCE FUNCTIONS

  model$intern$readfixed <- function(filename,first,last) {

    rawdata <- try(scan(filename,"character",sep="\n",quiet=TRUE))
    if ("try-error" %in% class(rawdata))
      stop(paste0("Cannot open file ", filename))

    iscomment <- grepl("^[[:blank:]]*#",rawdata)
    xsdata <- t(sapply(rawdata[!iscomment],substring,first,last,USE.NAMES=FALSE))
    xsdata <- apply(xsdata,2,as.numeric)
    if (is.null(dim(xsdata)))
      dim(xsdata) <- c(1, length(xsdata))
    return(xsdata)
  }

  model$intern$getProdXS <- function(calcdir,reacstr,spec,retfiles=FALSE) {

      residualCode <- regmatches(reacstr,regexec("^CS/PROD/([[:digit:]]{6})/TOT$",reacstr))[[1]][2]
      stopifnot(!is.na(residualCode))

      # construct filename for nuclide specification
      xsfile <- paste0("rp",residualCode,".tot")

      # if the residual is one of p,n,d,t,h,a
      # we need special cases
      if(residualCode=="001001") # proton
        xsfile <- "pprod.tot"
      if(residualCode=="000001") # neutron
        xsfile <- "nprod.tot"
      if(residualCode=="001002") # deuteron
        xsfile <- "dprod.tot"
      if(residualCode=="001003") # trition
        xsfile <- "tprod.tot"
      if(residualCode=="002003") # He-3
        xsfile <- "hprod.tot"
      if(residualCode=="002004") # alpha
        xsfile <- "aprod.tot"

      if (isTRUE(retfiles))
        return(xsfile)
      xsfile <- file.path(calcdir,xsfile)
      # retrieve data
      xsdata <- model$intern$readfixed(xsfile,c(1,13),c(12,24))
      approxWithSnap(xsdata[,1],xsdata[,2],xout=spec[,L1])$y
    }


  class(model) <- append(class(model),"model.TALYS")

  return(model)
}

