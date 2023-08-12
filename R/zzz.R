# Central import section

#' @importFrom assertthat assert_that has_args
#' @importFrom Biobase isVersioned isCurrent classVersion<- classVersion
#' @importFrom ChemmineR smiles2sdf validSDF write.SDF
#' @importFrom data.table fread fwrite
#' @import digest
#' @import glue
#' @import httr
#' @import logger
#' @import MSnbase
#' @import mzR
#' @import rcdk
#' @import Rcpp
#' @import readJDX
#' @import readr
#' @import rjson
#' @importFrom stats lm loess median predict smooth.spline
#' @import S4Vectors
#' @import tibble
#' @importFrom utils URLencode capture.output data flush.console 
#' @importFrom utils packageVersion read.csv read.csv2 setTxtProgressBar
#' @importFrom utils str txtProgressBar type.convert write.csv write.table
#' @importFrom utils globalVariables
#' @importFrom webchem cir_query
#' @import XML
#' @import yaml


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


