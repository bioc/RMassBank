
library(RMassBank)
library(RMassBankData)

loadList(system.file("list/NarcoticsDataset.csv", package="RMassBankData"))
storedW <- loadMsmsWorkspace(system.file("results/pH_narcotics_RF.RData",
package="RMassBankData"))
w <- storedW
loadRmbSettings("inst/RMB_options.ini")
mb <- newMbWorkspace(w)
mb <- loadInfolists(mb, system.file("infolists", package="RMassBankData"))
expect_silent(mb <- mbWorkflow(mb))

o <- getOption("RMassBank")
o$accessionBuilder <- "MASBNK-{contributor_prefix}-{entry_prefix}{compound_id(4)}{scan_id(2)}"
options(RMassBank = o)
expect_error(mb <- mbWorkflow(mb))

o$accessionValidate <- FALSE
options(RMassBank = o)
expect_silent(mb <- mbWorkflow(mb))


notaBuilder <- function(a) {
  stop("not a builder")
}

# notaBuilder should not be acceptable ACCESSION builder
expect_error(setAccessionBuilder(notaBuilder))

library(glue)
alternativeBuilder <- function(cpd, spectrum, subscan) {
  glue("MSBNK-ABCDE-12345")
}
# alternativeBuilder should be an acceptable ACCESSION builder
expect_silent(setAccessionBuilder(alternativeBuilder))
expect_silent(mb <- mbWorkflow(mb))
