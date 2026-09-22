# Ag-depth models with rbacon for sediment cores of lakes Castor, Gagnon, and Clo
# UTF-8
# Author: Dorian Gaboriau - dorian.gaboriau@uqat.ca
# Last update : September 2026

# Builds Bayesian age-depth models with the Bacon algorithm
# (Blaauw & Christen 2011) via the rbacon package

# Deleting variables from the environment, and closing any open graphics device

rm(list = ls())
if (!is.null(grDevices::dev.list())) dev.off()

####################
# Import libraries
####################
## ---- 1. Packages ---------------------------------------------------
required_pkgs <- c("rbacon", "ggplot2", "patchwork")
missing_pkgs  <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

library(rbacon)
library(ggplot2)
library(patchwork)

## ---- 2. Reproducibility & folders -----------------------------------
set.seed(785876478)

project_dir <- "D:/PROJETS_POSTDOC/Aiguebelle/Datasets/Files"   # <-- adjust if needed
setwd(project_dir)

# rbacon expects: <coredir>/<CoreName>/<CoreName>.csv
coredir <- file.path(project_dir, "Bacon_runs")
dir.create(coredir, showWarnings = FALSE)

## ---- 3. Helper: prepare input + run Bacon for one core --------------
run_bacon_core <- function(core_name, source_csv, d.min, d.max, d.by = 0.5) {
  
  core_folder <- file.path(coredir, core_name)
  dir.create(core_folder, showWarnings = FALSE)
  
  target_csv <- file.path(core_folder, paste0(core_name, ".csv"))
  if (!file.exists(target_csv) && file.exists(source_csv)) {
    file.copy(source_csv, target_csv, overwrite = TRUE)
  }
  
  # Sanity check: print the dates table before modelling
  cat("\n----", core_name, "input dates ----\n")
  print(read.csv(target_csv))
  
  Bacon(
    core_name,
    coredir     = coredir,
    d.min       = d.min,
    d.max       = d.max,
    d.by        = d.by,
    runname     = core_name,
    prob        = 0.95,
    rev.yr      = FALSE,
    rev.d       = TRUE,
    unit        = "cm",
    rotate.axes = TRUE,
    normal      = TRUE,
    ask         = FALSE      # run unattended, no y/n prompts
  )
  
  invisible(core_folder)
}


## ---- 4. Run the three cores ------------------------------------------
run_bacon_core("Castor", "Castor_Bacon.csv", d.min = 0, d.max = 510)
run_bacon_core("Gagnon", "Gagnon_Bacon.csv", d.min = 0, d.max = 380)
run_bacon_core("Clo",    "Clo_Bacon.csv",    d.min = 0, d.max = 207)
y



