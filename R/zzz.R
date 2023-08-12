#' @importFrom utils URLencode capture.output data flush.console 
#'   packageVersion read.csv read.csv2
#'   setTxtProgressBar str txtProgressBar type.convert write.csv write.table
#' @import digest
#' @import readr
#' @importFrom assertthat assert_that has_args
#' @import tibble
#' @import glue
#' @importFrom MSnbase writeMgfContent


# Auxiliary for getSplash.R so we can use the original file and don't have to change anything there

.onLoad <- function(libname, pkgname) {
  RMassBank.env <<- new.env()
  RMassBank.env$ReadAnnotation <- FALSE
  RMassBank.env$testnumber <- 1
  ## new variables
  RMassBank.env$verbose.output <- FALSE
  RMassBank.env$export.invalid <- FALSE
  RMassBank.env$export.molfiles <- TRUE
  RMassBank.env$strictMsMsSpectraSelection <- FALSE
  
  mb <- list()
  attach(RMassBank.env)
}

utils::globalVariables(c("cpdID", "isotopes","mzCalc"))


