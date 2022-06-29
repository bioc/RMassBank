library(testthat)
library(RMassBank)
library(RMassBankData)

loadList(system.file("list/NarcoticsDataset.csv", package="RMassBankData"))
storedW <- loadMsmsWorkspace(system.file("results/pH_narcotics_RF.RData",
package="RMassBankData"))
w <- storedW
w@spectra <- w@spectra[1:3]
loadRmbSettings("inst/RMB_options.ini")
mb <- newMbWorkspace(w)
mb <- loadInfolists(mb, system.file("infolists", package="RMassBankData"))
expect_silent(mb <- mbWorkflow(mb))

# MASBNK instead of MSBNK should fail the validator
setAccessionBuilder("MASBNK-{contributor_prefix}-{entry_prefix}{compound_id(4)}{scan_id(2)}")
expect_error(mb <- mbWorkflow(mb))

# Ignoring the validator should make the workflow work
o <- getOption("RMassBank")
o$accessionValidate <- FALSE
options(RMassBank = o)
expect_silent(mb <- mbWorkflow(mb))

o$accessionValidate <- TRUE
options(RMassBank = o)
setAccessionBuilder("MSBNK-{contributor_prefix}-{entry_prefix}{compound_id(4)}{scan_id(2)}")

# Test some accession_builder variants with 3 spectra only

mb <- loadInfolists(mb, system.file("infolists", package="RMassBankData"))
expect_silent(mb <- mbWorkflow(mb))


expect_silent(mb <- mbWorkflow(mb))
setAccessionBuilder("MSBNK-{contributor_prefix}-{info('INCHIKEY2D')}_{info('INCHIKEY')}{scan_id(2)}")
expect_silent(mb <- mbWorkflow(mb))
setAccessionBuilder("MSBNK-{contributor_prefix}-{compound_id(4)}_{mode}_{mode_hash}_{collision_energy_raw}")
expect_silent(mb <- mbWorkflow(mb))
setAccessionBuilder("MSBNK-{contributor_prefix}-{compound_id(4)}_{mode}_{mode_hash}_{condition_hash}_{polarity()}{polarity(5)}")
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
