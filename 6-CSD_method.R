# Charcoal Size Distribution (CSD) method to separate local from regional fires
# Author: Dorian Gaboriau - dorian.gaboriau@uqat.ca
# Last update : September 2026

## This script classifies charcoal peaks detected by CharAnalysis as either
## "Local" or "Regional" fires, using the Charcoal Size Distribution (CSD)
## method described in Asselin et Payette (2005) / Oris et al. (2014).

# Deleting variables from the environment, and closing any open graphics device
rm(list = ls())
if (!is.null(grDevices::dev.list())) dev.off()

####################
# Import libraries
####################
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(dplyr)
library(readxl)
library(tidyr)

## ---- 2. CONFIGURATION -------------------------------------------------------

# - charcoal_file : WinSeedle compiled charcoal measurements (semicolon-separated csv)
# - n_rows        : number of data rows to keep from charcoal_file (removes header
#                    artifact / trailing non-sediment sample, as in original script)
# - char_xls      : CharAnalysis output workbook (.xls), sheet "CharResults"
# - depth_xlsx    : correspondence table between master-core (MC) depth and
#                    intra-core ("DEPTH CORE") depth
# - lake_label    : short name appended to sample IDs (e.g. "Castors", "Gagnon")
#                   -> IMPORTANT: this must match exactly the suffix used in the
#                   SampleId column of the charcoal_file (e.g. "Castors", not "Castor")

lake_configs <- list(
  Castor = list(
    charcoal_file = "D:/PROJETS_POSTDOC/Aiguebelle/Datasets/Files/Lake_Castor/Castor_compil_CSD.csv",
    n_rows        = 7481,
    char_xls      = "D:/PROJETS_POSTDOC/Aiguebelle/Datasets/Files/Lake_Castor/CharAnalysis_Outputs/Castor.xls",
    depth_xlsx    = "D:/PROJETS_POSTDOC/Aiguebelle/Datasets/Files/Lake_Castor/Lac_Castor.xlsx",
    lake_label    = "Castors"
  ),
  Gagnon = list(
    charcoal_file = "D:/PROJETS_POSTDOC/Aiguebelle/Datasets/Files/Lake_Gagnon/Gagnon_compil_CSD.csv",
    n_rows        = 10000,
    char_xls      = "D:/PROJETS_POSTDOC/Aiguebelle/Datasets/Files/Lake_Gagnon/CharAnalysis_Outputs/Gagnon.xls",
    depth_xlsx    = "D:/PROJETS_POSTDOC/Aiguebelle/Datasets/Files/Lake_Gagnon/Lac_Gagnon.xlsx",
    lake_label    = "Gagnon"
  )
)


# Analysis parameters
MIN_PROJ_AREA   <- 0.025   # minimum particle projected area (WinSeedle protocol filter)
SLOPE_THRESHOLD <- -1.77   # CSD regression slope threshold: >= this value => "Local"
MIN_LOG_SIZE    <- -0.3    # a sample must contain >= 1 particle at/above this
# log10(geometric diameter, mm) to be eligible as "Local"

# Output directory for exported CSVs (adjust as needed)
OUTPUT_DIR <- "D:/PROJETS_POSTDOC/Aiguebelle/Datasets/"


## ---- 3. FUNCTIONS ------------------------------------------------------------

#' Compute the CSD regression slope for one charcoal sample
#' Bins log10(geometric diameter) values into six fixed size classes,
#' computes the proportion of particles in each class, and fits a linear
#' regression (proportion ~ class center). Empty edge/middle classes are
#' handled with small fallback values (0.01) or dropped, exactly as in the
#' original protocol, so that the regression can still be fit.
#'
#' @param sample_charpart numeric vector of log10(geometric diameter, mm)
#'   for all particles in one sample.
#' @return list(slope = <regression slope>, has_large_particle = <logical>)
compute_sample_slope <- function(sample_charpart) {
  
  # Fixed log-scale size classes (mm), each 0.2 log-units wide
  sc09 <- sample_charpart[sample_charpart >= -0.8 & sample_charpart < -0.6]
  sc07 <- sample_charpart[sample_charpart >= -0.6 & sample_charpart < -0.4]
  sc05 <- sample_charpart[sample_charpart >= -0.4 & sample_charpart < -0.2]
  sc03 <- sample_charpart[sample_charpart >= -0.2 & sample_charpart <  0.0]
  sc01 <- sample_charpart[sample_charpart >=  0.0 & sample_charpart <  0.2]
  sc02 <- sample_charpart[sample_charpart >=  0.2 & sample_charpart <  0.4]
  
  N <- length(sample_charpart)
  p09 <- length(sc09) / N
  p07 <- length(sc07) / N
  p05 <- length(sc05) / N
  p03 <- length(sc03) / N
  p01 <- length(sc01) / N
  p02 <- length(sc02) / N
  
  prop     <- c(p09, p07, p05, p03, p01, p02)
  CenterSC <- c(-0.7, -0.5, -0.3, -0.1, 0.1, 0.3)
  
  # --- Edge-case handling (kept identical to the original protocol) ---
  if (p09 == 0 & p07 != 0) {
    prop <- c(p07, p05, p03, p01, p02)
    CenterSC <- c(-0.5, -0.3, -0.1, 0.1, 0.3)
  }
  if (p07 == 0) {
    p07 <- 0.01
    prop <- c(p09, p07, p05, p03, p01, p02)
    CenterSC <- c(-0.7, -0.5, -0.3, -0.1, 0.1, 0.3)
  }
  if (p05 == 0 & p03 != 0) {
    p05 <- 0.01
    prop <- c(p09, p07, p05, p03, p01, p02)
  }
  if (p01 == 0 & p02 == 0) {
    prop <- c(p09, p07, p05, p03)
    CenterSC <- c(-0.7, -0.5, -0.3, -0.1)
  }
  if (p03 == 0 & p02 != 0 & p01 != 0) {
    p03 <- 0.01
    prop <- c(p09, p07, p05, p03, p01, p02)
  }
  if (p03 == 0 & p02 == 0 & p01 != 0) {
    p03 <- 0.01
    prop <- c(p09, p07, p05, p03, p01)
    CenterSC <- c(-0.7, -0.5, -0.3, -0.1, 0.1)
  }
  if (p02 == 0 & p01 != 0) {
    prop <- c(p09, p07, p05, p03, p01)
    CenterSC <- c(-0.7, -0.5, -0.3, -0.1, 0.1)
  }
  if (p02 != 0 & p01 == 0) {
    p01 <- 0.01
    prop <- c(p09, p07, p05, p03, p01, p02)
    CenterSC <- c(-0.7, -0.5, -0.3, -0.1, 0.1, 0.3)
  }
  if (p03 == 0 & p02 == 0 & p01 == 0) {
    prop <- c(p09, p07, p05)
    CenterSC <- c(-0.7, -0.5, -0.3)
  }
  if (p03 == 0 & p02 == 0 & p01 == 0 & p05 == 0) {
    prop <- c(p09, p07)
    CenterSC <- c(-0.7, -0.5)
  }
  
  reg   <- lm(prop ~ CenterSC)
  slope <- reg$coefficients[[2]]
  
  list(
    slope              = slope,
    has_large_particle = any(sample_charpart >= MIN_LOG_SIZE)
  )
}

#' Compute CSD slope and fire-type classification for every sample in a lake
#'
#' @param charcoal_file path to the WinSeedle compiled csv (";"-separated)
#' @param n_rows number of rows to keep (drops header artifact / trailing rows)
#' @return data.frame with columns SampleID, SlopeCoeff, MinCharcoalSizeOK, FireType
compute_csd_slopes <- function(charcoal_file, n_rows) {
  
  char_results <- read.csv(charcoal_file, sep = ";")
  char_results <- char_results[1:n_rows, ]
  
  # Keep only particles large enough per the WinSeedle measurement protocol
  seedle <- subset(char_results, ProjArea > MIN_PROJ_AREA)
  
  # id column name is inconsistent between input files in the original data
  # ("Id" for Castor, "id" for Gagnon) -- detect it automatically here.
  id_col <- if ("Id" %in% names(seedle)) "Id" else "id"
  ids    <- seedle[[id_col]]
  
  n_samples <- nrow(seedle)
  results <- data.frame(
    SampleID          = character(n_samples),
    SlopeCoeff        = numeric(n_samples),
    MinCharcoalSizeOK = logical(n_samples),
    stringsAsFactors  = FALSE
  )
  n_filled <- 0
  
  # Walk through rows; each block of consecutive rows sharing the same
  # SampleId is one sediment sample containing several particle measurements.
  i <- 1
  while (i <= max(ids)) {
    if (i < nrow(seedle) && seedle$SampleId[i] == seedle$SampleId[i + 1]) {
      sample <- subset(seedle, SampleId == seedle$SampleId[i])
      sample$ProjArea <- as.numeric(sample$ProjArea)
      
      # Geometric diameter (Oris et al., 2014): sqrt(area), then log10.
      # Uses the ProjArea column BY NAME (not by position) so this stays
      # correct regardless of column order in a given file.
      sample_charpart <- log10(sqrt(sample$ProjArea))
      
      csd <- compute_sample_slope(sample_charpart)
      
      n_filled <- n_filled + 1
      results$SampleID[n_filled]          <- as.character(seedle$SampleId[i])
      results$SlopeCoeff[n_filled]        <- csd$slope
      results$MinCharcoalSizeOK[n_filled] <- csd$has_large_particle
    }
    i <- i + 1
  }
  
  results <- results[seq_len(n_filled), ]
  
  # Classification rule: steep slope + at least one large particle => Local;
  # otherwise Regional.
  results$FireType <- ifelse(
    results$SlopeCoeff >= SLOPE_THRESHOLD & results$MinCharcoalSizeOK,
    "Local", "Regional"
  )
  
  results
}

#' Match CharAnalysis-detected fire peaks to master-core depths, then to
#' intra-core (sub-core) sample IDs, and label each with its CSD fire type.
#'
#' @param char_xls path to CharAnalysis output workbook
#' @param depth_xlsx path to MC <-> intra-core depth correspondence table
#' @param csd_results output of compute_csd_slopes() for the same lake
#' @param lake_label short lake name used to build sample IDs
#' @return list(top = data.frame, bottom = data.frame) each with columns
#'   MC_depth_cm, SampleID, FireType
match_fires_to_samples <- function(char_xls, depth_xlsx, csd_results, lake_label) {
  
  # Only keep confirmed fire peaks (Final peak == 1) from CharAnalysis
  char_fire <- read_excel(char_xls, sheet = "CharResults")
  char_fire <- subset(char_fire, char_fire$`peaks         Final` == 1)
  
  top_depth <- char_fire$`cm Top_i  (cm)`
  
  # Each detected fire's depth is rounded to the nearest 0.5 cm sampling
  # interval to find the bracketing top/bottom samples in the master core.
  round_to_half <- function(x) {
    rem <- x %% 1
    if (rem >= 0 & rem < 0.5) x - rem else x - rem + 0.5
  }
  top_sample_depth    <- vapply(top_depth, round_to_half, numeric(1))
  bottom_sample_depth <- top_sample_depth + 0.5
  
  # Correspondence table: master-core (MC) depth <-> sub-core ("DEPTH CORE")
  # depth, used to rebuild each sample's original SampleID.
  prof <- read_excel(depth_xlsx)
  
  # Formats a depth value the way it appears in the real SampleId strings:
  # a whole number must come out as "11", never "11.0" -- otherwise the
  # rebuilt ID ("A1_11_0_...") never matches the real one ("A1_11_...")
  # and the fire silently falls back to "Regional". This was the source of
  # the Castor misclassification.
  format_depth <- function(x) {
    if (is.na(x)) return(NA_character_)
    if (x == round(x)) as.character(as.integer(round(x))) else as.character(x)
  }
  
  build_sample_id_lookup <- function(mc_depths) {
    vapply(mc_depths, function(d) {
      match_row <- which(prof$DEPTH_SEDIMENT == d)
      if (length(match_row) == 0) return(NA_character_)
      depth_str <- format_depth(prof$`DEPTH CORE`[match_row[1]])
      id <- paste(prof$CORE[match_row[1]], depth_str, lake_label, sep = "_")
      sub("\\.", "_", id)
    }, character(1))
  }
  
  make_matched_df <- function(mc_depths) {
    sample_id <- build_sample_id_lookup(mc_depths)
    fire_type <- csd_results$FireType[match(sample_id, csd_results$SampleID)]
    fire_type <- tidyr::replace_na(fire_type, "Regional")
    data.frame(
      MC_depth_cm = mc_depths,
      SampleID    = sample_id,
      FireType    = fire_type,
      stringsAsFactors = FALSE
    )
  }
  
  list(
    top    = make_matched_df(top_sample_depth),
    bottom = make_matched_df(bottom_sample_depth)
  )
}

#' Plot Local vs Regional fires along core depth
plot_fire_types <- function(df, position, lake_label) {
  ggplot(df, aes(x = MC_depth_cm, y = FireType, group = FireType)) +
    geom_point(aes(shape = FireType, color = FireType)) +
    scale_color_manual(values = c('#440154FF', '#73D055FF')) +
    rotate() +
    ggtitle(paste0("Local and regional fires with CSD method \napplied to the ",
                   position, " sample - ", lake_label)) +
    scale_x_reverse(limits = c(500, 0), breaks = seq(500, 0, by = -100)) +
    theme_bw()
}

#' Run the full CSD workflow for one lake and return/export the results
run_csd_analysis <- function(cfg) {
  
  message("Processing lake: ", cfg$lake_label)
  
  csd_results <- compute_csd_slopes(cfg$charcoal_file, cfg$n_rows)
  
  matched <- match_fires_to_samples(
    char_xls    = cfg$char_xls,
    depth_xlsx  = cfg$depth_xlsx,
    csd_results = csd_results,
    lake_label  = cfg$lake_label
  )
  
  message("Fire type counts (bottom sample):")
  print(matched$bottom %>% group_by(FireType) %>% count())
  message("Fire type counts (top sample):")
  print(matched$top %>% group_by(FireType) %>% count())
  
  print(plot_fire_types(matched$bottom, "bottom", cfg$lake_label))
  print(plot_fire_types(matched$top, "top", cfg$lake_label))
  
  # Uncomment to export:
  # write.csv(matched$top,
  #           file = file.path(OUTPUT_DIR, paste0("csd_", cfg$lake_label, ".csv")),
  #           row.names = FALSE)
  
  invisible(matched)
}

## ---- 4. RUN FOR ALL LAKES ---------------------------------------------------
# Results are stored per lake in `all_results`, e.g. all_results$Castor$top
all_results <- lapply(lake_configs, run_csd_analysis)
