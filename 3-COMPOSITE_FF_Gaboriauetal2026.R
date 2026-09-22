# Fire regime Composite (for Fire Frequency))
# UTF-8
# Author: Dorian Gaboriau - dorian.gaboriau@uqat.ca
# Last update : September 2026

# Deleting variables from the environment, and closing any open graphics device

rm(list = ls())
if (!is.null(grDevices::dev.list())) dev.off()

#Launch Rscripts from 'Scripts_Paleofire' folder

####################
# Import libraries
####################

library(devtools)
library(ggplot2)
library(MASS)
library(cowplot)      # plot_grid(), get_legend()
library(boot)          # bootstrap confidence intervals
library(tidyverse)   # dplyr, readr, etc.
library(nlme)          # gls(), corAR1()
library(mgcv)          # gam()
library(scales)        # alpha() for transparent colors

#Working directory
setwd("D:/PROJETS_POSTDOC/Aiguebelle/Datasets/Files/")

## ============================================================================
## IMPORT DATA
## ============================================================================
# Each file contains, for one region (lowlands / hills), one column per lake
# with the dates (cal yr BP) of the significant fire peaks identified by
# CharAnalysis.

lakes_fires_low   <- read.csv("Lake_fire_frequency_lowv2.csv",   header = TRUE, sep = ";")
lakes_fires_hills <- read.csv("Lake_fire_frequency_hillsv2.csv", header = TRUE, sep = ";")

## ============================================================================
## FIRE FREQUENCY RECONSTRUCTION PER LAKE
## ============================================================================
# Fire frequency (number of fires per millennium) is estimated from a set of
# fire-event dates using a Gaussian kernel density estimation procedure with
# a fixed bandwidth (see Mudelsee 2004 for methodological details).
# pseudo = TRUE enables pseudo-replicated values to correct for edge bias,
# equivalent to the "minimum slope" correction described in Mann (2004).

#' Reconstruct fire frequency for a single lake and build its plot.
#'
#' @param fire_events numeric vector of fire event dates (cal yr BP)
#' @param lake_name    character, lake name (used as plot title)
#' @param line_color   color used for the frequency curve
#' @param up, lo       age bounds passed to kdffreq()
#' @param bandwidth, interval, nbboot, alpha  kdffreq() parameters
#' @param show_y_title logical, whether to display the y-axis title
#' @return a list with the raw kdffreq object, the tidy table, the mean
#'         frequency, and the ggplot object
compute_ff <- function(fire_events, lake_name, line_color,
                       up = -72, lo = 10000, bandwidth = 500, interval = 10,
                       nbboot = 1000, alpha = 0.01, show_y_title = FALSE) {
  
  ff <- kdffreq(fire_events, up = up, lo = lo, bandwidth = bandwidth,
                interval = interval, nbboot = nbboot, alpha = alpha, pseudo = TRUE)
  
  ff_table <- data.frame(
    Age  = ff$age,
    Freq = ff$ff * 1000,   # convert to "fires per millennium"
    Lo   = ff$lo * 1000,
    Up   = ff$up * 1000
  )
  mean_freq <- mean(ff_table$Freq)
  
  p <- ggplot(ff_table) +
    geom_line(aes(x = Age, y = Freq), colour = line_color, size = 1.3) +
    geom_ribbon(aes(x = Age, ymin = Lo, ymax = Up), linetype = 2, alpha = 0.3, fill = "grey") +
    geom_hline(yintercept = mean_freq, colour = "black", linetype = "dashed", size = 1) +
    scale_x_continuous(trans = "reverse", breaks = seq(0, 10000, 2000), limits = c(10000, 0)) +
    scale_y_continuous(breaks = c(0, 5, 10, 15), limits = c(0, 15)) +
    theme_classic() +
    labs(title = lake_name, y = if (show_y_title) "Number of fires per millennium" else NULL) +
    theme(
      panel.grid.minor.y = element_line(colour = "black"),
      legend.position     = "none",
      legend.title        = element_blank(),
      axis.title.x        = element_blank(),
      axis.title.y        = if (show_y_title) element_text(size = 10, face = "bold") else element_blank(),
      axis.text.x         = element_text(color = "black", size = 11, angle = 30),
      axis.text.y         = element_text(color = "black", size = 11, angle = 30),
      axis.line           = element_line(color = "black", size = 0.6, linetype = "solid")
    )
  
  list(raw = ff, table = ff_table, mean = mean_freq, plot = p)
}

# --- Hills lakes --------------------------------------------------------------
# The first hills lake keeps the y-axis title; the rest omit it (as in the
# original multi-panel figure).
hills_specs <- list(
  Castor     = list(events = lakes_fires_hills$Castor,     max_age = 9630,  min_age = 0,    y_title = TRUE),
  Gagnon     = list(events = lakes_fires_hills$Gagnon,     max_age = 9910,  min_age = 0,    y_title = FALSE),
  Perche     = list(events = lakes_fires_hills$Perche,     max_age = 9630,  min_age = 6500, y_title = FALSE),
  Desperiers = list(events = lakes_fires_hills$Desperiers, max_age = 9930,  min_age = 6670, y_title = FALSE),
  Labelle    = list(events = lakes_fires_hills$Labelle,    max_age = 10370, min_age = 0,    y_title = FALSE)
)

hills_results <- lapply(names(hills_specs), function(nm) {
  spec <- hills_specs[[nm]]
  compute_ff(spec$events, lake_name = nm, line_color = "brown", show_y_title = spec$y_title)
})
names(hills_results) <- names(hills_specs)

# --- Lowlands lakes ------------------------------------------------------------
low_specs <- list(
  Lili      = list(events = lakes_fires_low$Lili,      max_age = 8450, min_age = 0, interp_rule = 2,  y_title = TRUE),
  Pasdefond = list(events = lakes_fires_low$Pasdefond, max_age = 7700, min_age = 0, interp_rule = 2,  y_title = FALSE),
  Francis   = list(events = lakes_fires_low$Francis,   max_age = 7080, min_age = 0, interp_rule = 2,  y_title = FALSE),
  Clo       = list(events = lakes_fires_low$Clo,       max_age = 8500, min_age = 0, interp_rule = 10, y_title = FALSE)  # NB: rule=10 kept as in original script
)

low_results <- lapply(names(low_specs), function(nm) {
  spec <- low_specs[[nm]]
  compute_ff(spec$events, lake_name = nm, line_color = "steelblue", show_y_title = spec$y_title)
})
names(low_results) <- names(low_specs)

# --- Multi-panel figure: one panel per lake -----------------------------------
graph1 <- plot_grid(
  hills_results$Castor$plot, hills_results$Gagnon$plot, hills_results$Perche$plot,
  hills_results$Desperiers$plot, hills_results$Labelle$plot,
  low_results$Lili$plot, low_results$Pasdefond$plot, low_results$Francis$plot,
  low_results$Clo$plot,
  nrow = 2, ncol = 5
)
graph1

## ============================================================================
## LINEAR INTERPOLATION TO REGULAR TIME STEPS
## ============================================================================
# Each lake's fire-frequency series is interpolated onto a regular time axis
# (10-year steps) over its own valid age range, so that series from different
# lakes can later be merged and averaged.

interpolate_ff <- function(ff_result, min_age, max_age, rule = 2) {
  interp <- approx(ff_result$raw$age, ff_result$raw$ff, method = "linear",
                   xout = seq(min_age, max_age, by = 10), rule = rule)
  data.frame(Age = interp$x, FF = interp$y * 1000)
}

ff_hills_interp <- Map(function(res, spec) interpolate_ff(res, spec$min_age, spec$max_age),
                       hills_results, hills_specs)

ff_low_interp <- Map(function(res, spec) interpolate_ff(res, spec$min_age, spec$max_age, spec$interp_rule),
                     low_results, low_specs)

## ============================================================================
## 6. REGIONAL COMPOSITES (HILLS AND LOWLANDS)
## ============================================================================
# Merge the interpolated series of all lakes within a region and compute the
# regional mean fire frequency (row-wise mean across lakes, ignoring NAs).

build_composite <- function(interp_list) {
  composite <- Reduce(function(a, b) merge(a, b, by = "Age", all = TRUE),
                      lapply(names(interp_list), function(nm) {
                        df <- interp_list[[nm]]
                        colnames(df)[2] <- nm
                        df
                      }))
  composite$FFmean <- apply(composite[, -1, drop = FALSE], 1, mean, na.rm = TRUE)
  composite
}

FFhills <- build_composite(ff_hills_interp)
FFlow   <- build_composite(ff_low_interp)

head(FFhills)
head(FFlow)

# Quick overview plots of individual lake series within each region
plot(FFhills$Age, FFhills$Castor, xlim = rev(range(FFhills$Age)), type = "l",
     ylim = c(0, 15), lwd = 2, ylab = "", col = "blue", cex.axis = 1.3)
lines(FFhills$Age, FFhills$Gagnon,     col = "forestgreen", lwd = 2)
lines(FFhills$Age, FFhills$Perche,     col = "red",         lwd = 2)
lines(FFhills$Age, FFhills$Desperiers, col = "black",       lwd = 2)
lines(FFhills$Age, FFhills$Labelle,    col = "orange",      lwd = 2)

plot(FFlow$Age, FFlow$Lili, xlim = rev(range(FFlow$Age)), type = "l",
     ylim = c(0, 15), lwd = 2, ylab = "", col = "blue", cex.axis = 1.3)
lines(FFlow$Age, FFlow$Pasdefond, col = "forestgreen", lwd = 2)
lines(FFlow$Age, FFlow$Francis,   col = "red",         lwd = 2)
lines(FFlow$Age, FFlow$Clo,       col = "black",       lwd = 2)


## ============================================================================
## 7. BOOTSTRAP CONFIDENCE INTERVALS FOR THE REGIONAL MEAN
## ============================================================================
# For each time step, a 90% confidence interval around the mean fire
# frequency across lakes is estimated by non-parametric bootstrap.

#' @param lake_matrix data.frame with one column per lake (no Age/mean column),
#'        restricted to the time steps where at least two lakes have data.
#' @param conf confidence level
#' @param R number of bootstrap resamples
compute_boot_ci <- function(lake_matrix, conf = 0.90, R = 999) {
  boot_mean <- function(y, idx) mean(y[idx])
  transposed <- as.data.frame(t(lake_matrix))  # required layout for boot()
  
  ci_rows <- lapply(transposed, function(y) {
    y <- na.omit(y)
    boot_obj <- boot(y, boot_mean, R)
    boot.ci(boot_obj, conf = conf, type = "norm")$normal
  })
  as.data.frame(do.call(rbind, ci_rows))
}

# --- Hills: keep only time steps with at least two lakes (rows 1:994) --------
reg_ff_hills <- FFhills[1:994, c("Castor", "Gagnon", "Perche", "Desperiers", "Labelle")]
boot_ci_hills <- compute_boot_ci(reg_ff_hills)

# --- Lowlands: keep only time steps with at least two lakes (rows 1:846) ----
reg_ff_low <- FFlow[1:846, c("Lili", "Pasdefond", "Francis", "Clo")]
boot_ci_low <- compute_boot_ci(reg_ff_low)

# --- Attach CI bounds to the regional composite tables ------------------------
FFBootCIhills <- cbind(FFhills$Age[1:994], FFhills$FFmean[1:994], boot_ci_hills)
colnames(FFBootCIhills) <- c("Age", "Mean", "CI", "Infe", "Supe")

FFBootCIlow <- cbind(FFlow$Age[1:846], FFlow$FFmean[1:846], boot_ci_low)
colnames(FFBootCIlow) <- c("Age", "Mean", "CI", "Infe", "Supe")

# --- Smoothing splines (mean and CI bounds), for a cleaner regional curve ----
FFmean_smoothhills <- smooth.spline(FFBootCIhills$Age, FFBootCIhills$Mean, spar = 0.4)
FFmean_infhill      <- smooth.spline(FFBootCIhills$Age, FFBootCIhills$Infe, spar = 0.4)
FFmean_suphill      <- smooth.spline(FFBootCIhills$Age, FFBootCIhills$Supe, spar = 0.4)

FFmean_smoothlow <- smooth.spline(FFBootCIlow$Age, FFBootCIlow$Mean, spar = 0.4)
FFmean_inflow     <- smooth.spline(FFBootCIlow$Age, FFBootCIlow$Infe, spar = 0.4)
FFmean_suplow     <- smooth.spline(FFBootCIlow$Age, FFBootCIlow$Supe, spar = 0.4)

# Optionally export the regional composites for reuse without re-running
# the reconstruction above:
# write.csv(FFhills, file = "D:/FF_hills.csv", row.names = FALSE)
# write.csv(FFlow,   file = "D:/FF_low.csv",   row.names = FALSE)

## ============================================================================
## WILCOXON TESTS: HILLS vs LOWLANDS, BY PERIOD
## ============================================================================
# Non-parametric comparison of regional mean fire frequency between hills and
# lowlands, for the full Holocene and two sub-periods.

#' Build a boxplot comparing hills vs lowlands FFmean for a given row range.
make_comparison_boxplot <- function(row_range, title, show_p_annotation = FALSE,
                                    y_limits = NULL, keep_y_title = FALSE) {
  df <- data.frame(
    valeur = c(FFhills$FFmean[row_range], FFlow$FFmean[row_range]),
    groupe = factor(c(rep("FFhills", length(row_range)), rep("FFlow", length(row_range))))
  )
  
  p <- ggplot(df, aes(x = groupe, y = valeur, fill = groupe)) +
    geom_boxplot(alpha = 0.4) +
    scale_fill_manual(values = c(FFhills = "brown", FFlow = "steelblue"),
                      labels = c(FFhills = "Hills", FFlow = "Lowlands"), name = NULL) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title.x    = element_blank(),
      axis.title.y    = if (keep_y_title) element_text() else element_blank(),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      plot.title      = element_text(size = 12)
    ) +
    labs(y = "FFmean", title = title)
  
  if (!is.null(y_limits)) p <- p + scale_y_continuous(limits = y_limits)
  if (show_p_annotation) {
    p <- p + annotate("text", x = 1.5, y = 11, label = "p < 0.001", size = 4.2, fontface = "bold")
  }
  
  list(plot = p, test = wilcox.test(df$valeur[df$groupe == "FFhills"],
                                    df$valeur[df$groupe == "FFlow"], paired = FALSE))
}

# Full Holocene (8,450-0 cal yr BP)
holocene_cmp <- make_comparison_boxplot(1:4226, "Holocene (8450 - 0)", keep_y_title = TRUE)
box1 <- holocene_cmp$plot
w1   <- holocene_cmp$test
box1; w1

# Late Holocene (4,200-0 cal yr BP)
late_cmp <- make_comparison_boxplot(1:421, "Late Holocene (4,200-0 cal. yr BP)",
                                    show_p_annotation = TRUE, y_limits = c(0, 12))
box2 <- late_cmp$plot
w2   <- late_cmp$test
box2; w2

# Mid-Holocene (8,200-4,200 cal yr BP)
mid_cmp <- make_comparison_boxplot(421:821, "Mid-Holocene (8,200-4,200 cal. yr BP)",
                                   show_p_annotation = TRUE, y_limits = c(0, 12))
box3 <- mid_cmp$plot
w3   <- mid_cmp$test
box3; w3


## ============================================================================
## REGIONAL FIRE-FREQUENCY CURVE (HILLS vs LOWLANDS, WITH CI RIBBONS)
## ============================================================================

df_hill <- data.frame(
  Age   = FFmean_smoothhills$x,
  Mean  = FFmean_smoothhills$y,
  Lower = approx(FFmean_infhill$x, FFmean_infhill$y, xout = FFmean_smoothhills$x, rule = 2)$y,
  Upper = approx(FFmean_suphill$x, FFmean_suphill$y, xout = FFmean_smoothhills$x, rule = 2)$y,
  zone  = "Hill"
)

df_low <- data.frame(
  Age   = FFmean_smoothlow$x,
  Mean  = FFmean_smoothlow$y,
  Lower = approx(FFmean_inflow$x, FFmean_inflow$y, xout = FFmean_smoothlow$x, rule = 2)$y,
  Upper = approx(FFmean_suplow$x, FFmean_suplow$y, xout = FFmean_smoothlow$x, rule = 2)$y,
  zone  = "Low"
)

linemeanFFhill <- mean(df_hill$Mean, na.rm = TRUE)
linemeanFFlow  <- mean(df_low$Mean,  na.rm = TRUE)

df_all <- rbind(df_hill, df_low)

cols_line <- c(Hill = "brown",      Low = "steelblue")
cols_fill <- c(Hill = alpha("brown", 0.20), Low = alpha("steelblue", 0.20))

# Restrict each zone to its valid age range
df_all_filt <- df_all |>
  dplyr::filter(
    (zone == "Hill" & Age >= 0 & Age <= 9900) |
      (zone == "Low"  & Age >= 0 & Age <= 8450)
  )

p <- ggplot(df_all_filt, aes(x = Age, y = Mean, color = zone, fill = zone)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), color = NA) +
  geom_line(linewidth = 1.3) +
  geom_hline(yintercept = linemeanFFhill, color = "brown",     linetype = "dashed") +
  geom_hline(yintercept = linemeanFFlow,  color = "steelblue", linetype = "dashed") +
  scale_x_reverse(limits = c(10000, 0), breaks = seq(0, 10000, by = 1000), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 11), expand = c(0, 0)) +
  scale_color_manual(values = cols_line, labels = c(Hill = "Hills", Low = "Lowlands"), name = NULL) +
  scale_fill_manual(values = cols_fill,  labels = c(Hill = "Hills", Low = "Lowlands"), name = NULL) +
  labs(x = "Age (cal. yr BP)", y = "Number of fires per millennium", color = NULL, fill = NULL) +
  theme_classic() +
  theme(
    axis.title.x    = element_text(size = 13),
    axis.title.y    = element_text(size = 12, face = "bold"),
    axis.text       = element_text(size = 12, color = "black"),
    axis.ticks      = element_line(color = "black"),
    axis.line       = element_line(linewidth = 1),
    plot.margin     = margin(10, 10, 10, 10),
    legend.position = "top",
    legend.text     = element_text(size = 12)
  )
p


## ============================================================================
## FINAL COMPOSITE FIGURE
## ============================================================================

right_boxes <- plot_grid(box3, box2, ncol = 1, rel_heights = c(1, 1), align = "v")
bottom_row  <- plot_grid(p, right_boxes, ncol = 2, rel_widths = c(3, 1), align = "h")
final_fig   <- plot_grid(graph1, bottom_row, ncol = 1, rel_heights = c(1, 1.3))
final_fig

## ============================================================================
## BINNED CORRELATIONS (500-YEAR BINS): HILLS vs LOWLANDS
## ============================================================================
# Regional FF series are averaged into 500-year bins and compared with 95%
# quantile-based intervals, for the full record and two sub-periods.

#' Bin a regional composite table into fixed-width age classes and summarize
#' the mean fire frequency with 95% quantile bounds per bin.
bin_ff_summary <- function(df, age_min, age_max, bin_width = 500) {
  df %>%
    dplyr::filter(Age >= age_min & Age <= age_max) %>%
    dplyr::mutate(bin = cut(Age, breaks = seq(age_min, age_max, bin_width),
                            include.lowest = TRUE, right = FALSE)) %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(
      mean  = mean(FFmean, na.rm = TRUE),
      lower = quantile(FFmean, 0.025, na.rm = TRUE),
      upper = quantile(FFmean, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

#' Build a Hills-vs-Lowlands biplot with error bars from two binned summaries.
make_biplot <- function(ic_hills, ic_low, title) {
  df_bi <- ic_hills %>%
    dplyr::rename(FFmean_Hills = mean, Hills_low = lower, Hills_up = upper) %>%
    dplyr::inner_join(
      ic_low %>% dplyr::rename(FFmean_Lowlands = mean, Low_low = lower, Low_up = upper),
      by = "bin"
    )
  
  ggplot(df_bi, aes(x = FFmean_Hills, y = FFmean_Lowlands)) +
    geom_errorbarh(aes(xmin = Hills_low, xmax = Hills_up), height = 0.05, color = "grey40") +
    geom_errorbar(aes(ymin = Low_low, ymax = Low_up), width = 0.05, color = "grey40") +
    geom_point(size = 3, color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
    theme_minimal(base_size = 10) +
    labs(x = "FFmean Hills", y = "FFmean Lowlands", title = title)
}

# Full record (0-8,500 cal yr BP)
ic_hills_full <- bin_ff_summary(FFhills, 0, 8500)
ic_low_full   <- bin_ff_summary(FFlow,   0, 8500)
biplotff1_IC  <- make_biplot(ic_hills_full, ic_low_full, "FFmean (0-8,500 cal. yr BP)")
biplotff1_IC

# 0-4,000 cal yr BP
ic_hills_0_4000 <- bin_ff_summary(FFhills, 0, 4000)
ic_low_0_4000   <- bin_ff_summary(FFlow,   0, 4000)
biplotff2_IC    <- make_biplot(ic_hills_0_4000, ic_low_0_4000, "FFmean (0-4,000 cal. yr BP)")
biplotff2_IC

# 4,000-8,000 cal yr BP
ic_hills_4000_8000 <- bin_ff_summary(FFhills, 4000, 8000)
ic_low_4000_8000   <- bin_ff_summary(FFlow,   4000, 8000)
biplotff3_IC       <- make_biplot(ic_hills_4000_8000, ic_low_4000_8000, "FFmean (4,000-8,000 cal. yr BP)")
biplotff3_IC

final_biplot <- plot_grid(biplotff1_IC, biplotff2_IC, biplotff3_IC,
                          ncol = 3, labels = c("A", "B", "C"), label_size = 12)
final_biplot

# --- Same full-record biplot, but colored by age class (4,200-0 vs 8,200-4,200) ---
# This is the version used for the final "RegFF" figure: each 500-year bin is
# colored according to whether it falls in the late-Holocene (blue) or
# mid-Holocene (red) window, instead of using a single color/trend line.
df_bi_full <- ic_hills_full %>%
  dplyr::rename(FFmean_Hills = mean, Hills_low = lower, Hills_up = upper) %>%
  dplyr::inner_join(
    ic_low_full %>% dplyr::rename(FFmean_Lowlands = mean, Low_low = lower, Low_up = upper),
    by = "bin"
  ) %>%
  dplyr::mutate(
    bin_start = readr::parse_number(as.character(bin)),
    AgeClass  = dplyr::case_when(
      bin_start < 4200 ~ "4200-0",
      TRUE              ~ "8200-4200"
    )
  )

biplotff1_byclass <- ggplot(df_bi_full, aes(x = FFmean_Hills, y = FFmean_Lowlands, color = AgeClass)) +
  geom_errorbarh(aes(xmin = Hills_low, xmax = Hills_up), height = 0.05) +
  geom_errorbar(aes(ymin = Low_low, ymax = Low_up), width = 0.05) +
  geom_point(size = 3) +
  scale_color_manual(values = c("8200-4200" = "red", "4200-0" = "blue"), name = " ") +
  theme_minimal(base_size = 14) +
  labs(x = "Hills", y = "Lowlands", title = "RegFF")
biplotff1_byclass

## ============================================================================
## GLS + AR1 CORRELATIONS ON DETRENDED (GAM RESIDUAL) SERIES
## ============================================================================
# The long-term trend is removed from each regional series with a GAM, and
# the residuals are compared using a generalized least-squares model with an
# AR1 correlation structure (to account for temporal autocorrelation), plus
# a Spearman correlation as a robustness check.

fit_hills <- gam(FFhills$FFmean ~ s(FFhills$Age, k = 50))
fit_low   <- gam(FFlow$FFmean   ~ s(FFlow$Age,   k = 50))

res_hills <- resid(fit_hills)
res_low   <- resid(fit_low)

pred_hills <- predict(fit_hills, se.fit = TRUE)
pred_low   <- predict(fit_low,   se.fit = TRUE)

#' Run the GLS+AR1 model, Spearman test, and residual scatterplot for a given
#' row range (time period) of the detrended series.
run_period_analysis <- function(row_range, title) {
  df <- data.frame(
    age      = FFhills$Age[row_range],
    hills    = res_hills[row_range],
    lowlands = res_low[row_range],
    hills_se = pred_hills$se.fit[row_range],
    low_se   = pred_low$se.fit[row_range]
  )
  
  model_det <- gls(hills ~ lowlands, correlation = corAR1(form = ~ age), data = df)
  spearman  <- cor.test(df$hills, df$lowlands, method = "spearman")
  
  p <- ggplot(df, aes(x = lowlands, y = hills)) +
    geom_errorbar(aes(ymin = hills - 1.96 * hills_se, ymax = hills + 1.96 * hills_se),
                  width = 0, alpha = 0.3) +
    geom_errorbarh(aes(xmin = lowlands - 1.96 * low_se, xmax = lowlands + 1.96 * low_se),
                   height = 0, alpha = 0.3) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10)) +
    labs(x = "Lowlands (residuals from GAM)", y = "Hills (residuals from GAM)", title = title)
  
  list(model = model_det, spearman = spearman, plot = p)
}

full_period <- run_period_analysis(1:846,   "RegFF - 8,450 - 0 cal. BP")
mid_period  <- run_period_analysis(421:821, "RegFF - 8,200 - 4,200 cal. BP")
late_period <- run_period_analysis(1:421,   "RegFF - 4,200 - 0 cal. BP")

summary(full_period$model); full_period$spearman
summary(mid_period$model);  mid_period$spearman
summary(late_period$model); late_period$spearman

A_plot <- full_period$plot
B_plot <- mid_period$plot
C_plot <- late_period$plot

combined_plot <- plot_grid(A_plot, B_plot, C_plot, labels = c("A", "B", "C"),
                           label_size = 12, ncol = 3, align = "hv")
combined_plot

## ============================================================================
## SIMPLE SPEARMAN CORRELATIONS (NON-DETRENDED SERIES)
## ============================================================================

# Full Holocene
cor(FFhills$FFmean[1:846], FFlow$FFmean[1:846], method = "spearman", use = "complete.obs")
cor.test(FFhills$FFmean[1:846], FFlow$FFmean[1:846], method = "spearman")

# Mid-Holocene
cor(FFhills$FFmean[421:821], FFlow$FFmean[421:821], method = "spearman", use = "complete.obs")
cor.test(FFhills$FFmean[421:821], FFlow$FFmean[421:821], method = "spearman")

# Late Holocene
cor(FFhills$FFmean[1:421], FFlow$FFmean[1:421], method = "spearman", use = "complete.obs")
cor.test(FFhills$FFmean[1:421], FFlow$FFmean[1:421], method = "spearman")

