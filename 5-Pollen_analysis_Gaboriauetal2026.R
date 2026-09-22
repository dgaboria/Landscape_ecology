# ============================================================================
# Pollen diagrams for lakes on hilltops vs. lowlands
# Author: Dorian Gaboriau - dorian.gaboriau@uqat.ca
# Last update: September 2026
# ============================================================================

# What this script does, for each of 9 taxa (Abies, Betula, Picea, Pinus
# banksiana, Pinus strobus, Cupressaceae, Alnus, Populus, Quercus):
#   1. Plots the pollen (%) time series for hilltop vs. lowland lakes
#   2. Builds a hilltop-vs-lowland biplot with error bars, coloured by
#      time period (8200-4200 vs 4200-0 cal. yr BP)
#   3. Detrends both series with a GAM, then tests whether the residuals
#      are correlated (GLS with an AR1 correlation structure)
#   4. Runs a Wilcoxon (Mann-Whitney) test + boxplot comparing hilltop vs.
#      lowland values
#   5. Assembles everything into combined figures
#

# Deleting variables from the environment, and closing any open graphics device

rm(list = ls())
if (!is.null(grDevices::dev.list())) dev.off()


## ---- 1. Libraries ---------------------------------------------------------
library(tidyverse)   # dplyr, tidyr, readr, ggplot2
library(nlme)        # gls(), corAR1()
library(mgcv)        # gam()
library(cowplot)     # plot_grid()
library(patchwork)   # wrap_plots(), "/" and "&" layout operators


## ---- 2. File locations -----------------------------------------------------
dir <- "D:/PROJETS_POSTDOC/Aiguebelle/Datasets/Files/Pollen"


## ---- 3. Per-taxon configuration --------------------------------------------
# main_path : file used for the time-series plot and the GAM detrending
#             (semicolon-separated, comma decimal)
# bi_path   : file used only for the biplot
# bi_reader : "readr" (default comma/dot format) or "semicolon_dot"
# bi_x/bi_y : column names for hilltop/lowland values in bi_path
# n         : number of rows (from the top) used in the GAM detrending
# show_xlab : TRUE for taxa in the middle column of the 3x3 grid, which
#             carry the shared "cal. yr BP" x-axis title
# xlim_detrend : optional x-axis limit override for the residuals plot
species_list <- list(
  list(letter = "A", gam_title = "Abies balsamea",
       ts_title = expression(italic("A. Abies balsamea")),
       main_path = file.path(dir, "Balsamfir.csv"),
       bi_path   = file.path(dir, "Abies.csv"), bi_reader = "readr",
       bi_x = "High", bi_y = "Low", bi_title = "Abies",
       n = 17, show_xlab = FALSE),
  
  list(letter = "B", gam_title = "Betula spp.",
       ts_title = expression(italic("B. Betula") ~ " spp."),
       main_path = file.path(dir, "Betula.csv"),
       bi_path   = file.path(dir, "Birch.csv"), bi_reader = "readr",
       bi_x = "High", bi_y = "Low", bi_title = "Birch",
       n = 17, show_xlab = TRUE),
  
  list(letter = "C", gam_title = "Picea spp.",
       ts_title = expression(italic("C. Picea") ~ " spp."),
       main_path = file.path(dir, "Picea.csv"),
       bi_path   = file.path(dir, "Spruce.csv"), bi_reader = "readr",
       bi_x = "High", bi_y = "Low", bi_title = "Spruce",
       n = 17, show_xlab = FALSE),
  
  list(letter = "D", gam_title = "Pinus banksiana",
       ts_title = expression(italic("D. Pinus banksiana")),
       main_path = file.path(dir, "Pinusbanksiana.csv"),
       bi_path   = file.path(dir, "Jackpine.csv"), bi_reader = "readr",
       bi_x = "High", bi_y = "Low", bi_title = "Jack pine",
       n = 17, show_xlab = FALSE),
  
  list(letter = "E", gam_title = "Pinus strobus",
       ts_title = expression(italic("E. Pinus strobus")),
       main_path = file.path(dir, "Pinusstrobus.csv"),
       bi_path   = file.path(dir, "Whitepine.csv"), bi_reader = "readr",
       bi_x = "High", bi_y = "Low", bi_title = "White pine",
       n = 17, show_xlab = TRUE),
  
  list(letter = "F", gam_title = "Cupressaceae",
       ts_title = expression(italic("F. Cupressaceae")),
       main_path = file.path(dir, "Thuja.csv"),
       bi_path   = file.path(dir, "Cedar.csv"), bi_reader = "semicolon_dot",
       bi_x = "High", bi_y = "Low", bi_title = "Cupressaceae",
       n = 18, show_xlab = FALSE),
  
  list(letter = "G", gam_title = "Alnus incana",
       ts_title = expression(italic("G. Alnus incana") ~ " subsp. " ~ italic("rugosa")),
       main_path = file.path(dir, "Alder.csv"),
       bi_path   = file.path(dir, "Alnus.csv"), bi_reader = "readr",
       bi_x = "Hilltops", bi_y = "Lowlands", bi_title = "Alnus",
       n = 18, show_xlab = FALSE, xlim_detrend = c(NA, 10)),
  
  list(letter = "H", gam_title = "Populus tremuloides",
       ts_title = expression(italic("H. Populus tremuloides")),
       main_path = file.path(dir, "Aspen.csv"),
       bi_path   = file.path(dir, "Populus.csv"), bi_reader = "readr",
       bi_x = "Hilltops", bi_y = "Lowlands", bi_title = "Trembling aspen",
       n = 18, show_xlab = TRUE),
  
  list(letter = "I", gam_title = "Quercus spp.",
       ts_title = expression(italic("I. Quercus") ~ " spp."),
       main_path = file.path(dir, "Quercus.csv"),
       bi_path   = file.path(dir, "Oak.csv"), bi_reader = "readr",
       bi_x = "Hilltops", bi_y = "Lowlands", bi_title = "Quercus",
       n = 17, show_xlab = FALSE)
)
names(species_list) <- vapply(species_list, `[[`, character(1), "letter")



## ---- 4. Wilcoxon test datasets (hilltop vs. lowland) -----------------------

wilcoxon_data <- list(
  A = read.table(header = TRUE, text = "
X A_High B_Low
0 2.513641242 1.563507364
500 2.297010566 2.298662026
1000 2.771204767 2.326814426
1500 3.599757374 2.237197685
2000 3.21089901 1.536924496
2500 1.380422533 1.959496695
3000 1.371120367 1.937831757
3500 1.729662555 1.245657664
4000 1.300593516 1.625247995
4500 0.786810706 1.199783153
5000 1.714634328 1.207627418
5500 0.674267328 0.900001522
6000 0.670555793 1.19677098
6500 0.984616892 0.59455915
7000 0.837615087 0.60861365
7500 1.191886642 0.406836237
8000 1.73705115 0.659955078"),
  
  B = read.table(header = TRUE, text = "
X A_High B_Low
0 28.4087045 17.59529554
500 25.2954113 29.90382742
1000 30.20449255 28.64795518
1500 27.82150767 25.5384135
2000 27.77671176 28.80598278
2500 28.18394488 31.69360695
3000 26.21348337 33.7277448
3500 28.99837565 34.93956472
4000 29.63972699 31.49759241
4500 27.76354647 29.5856302
5000 17.74934613 22.69513433
5500 15.48105951 19.00716224
6000 9.905423207 18.27490754
6500 15.30721533 23.41639187
7000 16.85642363 23.51633479
7500 22.08455839 19.3437611
8000 30.75597023 21.35770858
8500 14.96896249 30.4281457"),
  
  C = read.table(header = TRUE, text = "
X A_High B_Low
0 28.04099194 31.47097785
500 24.23164817 25.90496824
1000 23.6846839 24.42771257
1500 19.84395239 23.30841213
2000 18.45030452 25.04387686
2500 18.78269099 21.876522
3000 17.36997407 17.58839693
3500 13.18544609 14.90006627
4000 12.6964916 14.97380208
4500 11.15879415 11.52101096
5000 10.25004599 11.80978161
5500 8.727767834 11.19208532
6000 11.02193377 8.631433815
6500 9.666927652 8.757876225
7000 13.5518834 9.029602352
7500 11.07998253 7.977971478
8000 8.490234047 9.635006297
8500 13.46953317 9.923509934"),
  
  D = read.table(header = TRUE, text = "
X A_High B_Low
0 20.08550276 17.25781998
500 21.5116998 8.436561431
1000 17.79647773 9.308064745
1500 15.46090051 7.769259795
2000 13.89465869 9.317345377
2500 15.16598459 7.685692467
3000 15.22792739 8.689703027
3500 17.1683304 5.052212823
4000 13.79521879 4.573380934
4500 17.6164554 8.198419193
5000 18.12787251 7.787996622
5500 18.26040038 6.142354421
6000 15.47108487 7.744934025
6500 23.16105601 10.01418915
7000 28.88278054 18.48142262
7500 26.69733427 18.47594308
8000 23.61505685 23.58017917
8500 36.2604545 22.68741722"),
  
  E = read.table(header = TRUE, text = "
X A_High B_Low
0 5.068718487 4.946180899
500 8.938587922 5.359108787
1000 8.854168169 7.494627266
1500 15.3874972 7.210056083
2000 17.01499559 7.851944064
2500 19.9037516 9.567147076
3000 22.376121 9.872831788
3500 19.69427994 10.91432391
4000 21.36794568 12.97608959
4500 22.86048037 12.75167657
5000 29.12188264 15.75667233
5500 28.05379916 15.11015116
6000 39.56453605 22.10951449
6500 31.12282512 17.43291389
7000 20.62034285 13.54176736
7500 16.01861546 10.48892636
8000 4.174725225 4.173963627
8500 3.41261739 4.22781457"),
  
  F = read.table(header = TRUE, text = "
X A_High B_Low
0 0.422158729 1.798626627
500 0.328975858 2.277778642
1000 0.669063443 2.628521104
1500 0.954308334 3.132929903
2000 1.067213248 3.159258117
2500 0.606419709 3.360419285
3000 0.838673414 7.461788042
3500 0.701076243 6.324185663
4000 3.9306775 13.40017964
4500 2.245823671 10.67274016
5000 0.907132838 11.67290294
5500 1.330677561 16.3488249
6000 1.313737141 6.030366267
6500 1.532037427 2.360708507
7000 0.787837768 2.423091064
7500 0.780037313 2.479664209
8000 0.915716786 3.159348267
8500 0.937383163 5.570198676"),
  
  G = read.table(header = TRUE, text = "
X A_High B_Low
0 2.478759381 8.200145966
500 3.299767797 5.194860743
1000 2.53754871 6.478350543
1500 2.493218104 9.016276065
2000 3.162984461 7.962985828
2500 2.26546906 8.276479579
3000 3.185974632 6.640267549
3500 2.112836961 6.994488452
4000 1.951700033 9.083409203
4500 2.788686211 7.930426943
5000 3.365800376 7.553268875
5500 1.906488531 7.452559767
6000 0.577321067 5.939496904
6500 0.62589739 4.357553859
7000 1.198823091 1.944242988
7500 1.901447955 1.583300527
8000 3.309451333 0.62205692
8500 2.576133297 2.246688742"),
  
  H = read.table(header = TRUE, text = "
X A_High B_Low
0 0.912996539 1.655422
500 0.268629307 0.663566808
1000 0.621391634 0.423464294
1500 0.498933315 0.766328056
2000 0.428340616 1.060645767
2500 0.397459054 0.614396915
3000 0.539484593 0.932137138
3500 0.983424156 1.766052985
4000 1.224742337 2.331923372
4500 1.05437057 1.496848912
5000 1.438874118 1.353050837
5500 1.012002311 3.144495581
6000 1.318919126 5.099654636
6500 2.346140341 6.439966025
7000 2.656105161 5.885159216
7500 4.154522957 7.18139566
8000 5.77597956 7.704700503
8500 2.497683179 8.989072848"),
  
  I = read.table(header = TRUE, text = "
X A_High B_Low
0 1.223393139 1.700549058
500 1.29176578 1.184923324
1000 1.452792868 2.009416145
1500 1.010812628 1.772113982
2000 1.231718772 1.331370228
2500 0.930425546 1.327938242
3000 0.815784283 1.788256106
3500 1.425123412 1.890146473
4000 1.508951825 1.433845603
4500 1.210989863 1.113324697
5000 0.904772154 1.454285286
5500 1.324683802 1.240470567
6000 1.015573043 1.806940321
6500 1.801538035 1.791408068
7000 2.745903571 1.859238057
7500 2.074234473 2.125298424
8000 2.435744761 2.471001974")
)


## ---- 5. Helper functions ----------------------------------------------------

# Read the "main" file and standardise its 7 columns regardless of whether
# the original header said High/Low or Hilltops/Lowlands.
read_main_data <- function(path) {
  df <- read.csv(path, sep = ";", dec = ",", stringsAsFactors = FALSE)[, 1:7]
  colnames(df) <- c("Age", "High", "Low", "MinL", "MaxL", "MinH", "MaxH")
  df
}

# Read the (separate) file used only for the biplot.
read_biplot_data <- function(path, reader = c("readr", "semicolon_dot")) {
  reader <- match.arg(reader)
  if (reader == "readr") {
    readr::read_csv(path, show_col_types = FALSE)
  } else {
    read.csv(path, sep = ";", dec = ".", stringsAsFactors = FALSE)
  }
}

# Time-series plot: hilltop vs. lowland pollen (%) through time.
build_timeseries_plot <- function(df, title, show_xlab) {
  ggplot(df) +
    geom_line(aes(Age, High), colour = "brown", linewidth = 1.3) +
    geom_point(aes(Age, High), colour = "brown", size = 2) +
    geom_line(aes(Age, Low), colour = "steelblue", linewidth = 1.3) +
    geom_point(aes(Age, Low), colour = "steelblue", size = 2) +
    geom_ribbon(aes(Age, ymin = MinL, ymax = MaxL), linetype = 2, alpha = 0.3, fill = "steelblue") +
    geom_ribbon(aes(Age, ymin = MinH, ymax = MaxH), linetype = 2, alpha = 0.3, fill = "brown") +
    ggtitle(title) +
    scale_x_continuous(trans = "reverse", breaks = seq(9000, 0, -1000),
                       name = if (show_xlab) "cal. yr BP" else NULL) +
    scale_y_continuous(name = "%") +
    theme_classic() +
    theme(
      plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x     = element_text(color = "black", size = 11, angle = 30, vjust = 0.5),
      axis.text.y     = element_text(face = "bold", color = "black", size = 14, angle = 30),
      axis.title.x    = if (show_xlab) element_text(color = "black", size = 14, face = "bold") else element_blank(),
      axis.title.y    = element_text(color = "black", size = 14, face = "bold.italic"),
      axis.line       = element_line(color = "black", linewidth = 0.6),
      panel.grid.minor.y = element_line(colour = "black"),
      legend.position = "none"
    )
}

# Hilltop-vs-lowland biplot with error bars, coloured by time period.
# (Only the final, "kept" version from the original is reproduced.)
build_biplot <- function(df, x_col, y_col, title, xlim = NULL) {
  df <- df %>%
    dplyr::filter(Age <= 8500) %>%
    tidyr::drop_na(dplyr::all_of(y_col)) %>%
    dplyr::mutate(AgeClass = dplyr::case_when(
      Age >= 4200 & Age <= 8200 ~ "8200\u20134200",
      Age < 4200                ~ "4200\u20130",
      TRUE                      ~ NA_character_
    )) %>%
    dplyr::filter(!is.na(AgeClass), !is.na(MinH), !is.na(MaxH), !is.na(MinL), !is.na(MaxL))
  
  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color = AgeClass)) +
    geom_errorbarh(aes(xmin = MinH, xmax = MaxH), height = 0.05) +
    geom_errorbar(aes(ymin = MinL, ymax = MaxL), width = 0.05) +
    geom_point(size = 3) +
    scale_color_manual(values = c("8200\u20134200" = "red", "4200\u20130" = "blue"), name = NULL) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(size = 10)) +
    labs(x = "Hills", y = "Lowlands", title = title)
  
  if (!is.null(xlim)) p <- p + scale_x_continuous(limits = xlim)
  p
}


# GAM-detrend both series, then GLS(AR1) on the residuals + residual plot.
build_detrend_plot <- function(main_df, n, title, xlim = NULL) {
  age <- main_df$Age[1:n]
  fit_hills <- gam(main_df$High[1:n] ~ s(age, k = 5))
  fit_low   <- gam(main_df$Low[1:n]  ~ s(age, k = 5))
  
  pred_hills <- predict(fit_hills, se.fit = TRUE)
  pred_low   <- predict(fit_low,   se.fit = TRUE)
  
  df <- data.frame(
    age      = age,
    hills    = resid(fit_hills),
    lowlands = resid(fit_low),
    hills_se = pred_hills$se.fit,
    low_se   = pred_low$se.fit
  )
  
  model <- gls(hills ~ lowlands, correlation = corAR1(form = ~ age), data = df)
  
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
  
  if (!is.null(xlim)) p <- p + scale_x_continuous(limits = xlim)
  list(model = model, plot = p)
}

# Wilcoxon (Mann-Whitney) test + boxplot comparing hilltop vs. lowland values.
# Rationale (unpaired series measured at matching depths -> Mann-Whitney U).
build_wilcoxon_boxplot <- function(wide_df, letter, p_cutoff = 0.001) {
  test <- wilcox.test(wide_df$A_High, wide_df$B_Low, paired = FALSE)
  label <- if (test$p.value < p_cutoff) "p < 0.001" else paste0("p = ", signif(test$p.value, 3))
  
  long_df <- data.frame(
    Group  = factor(rep(c("High", "Low"), each = nrow(wide_df)), levels = c("High", "Low")),
    Value  = c(wide_df$A_High, wide_df$B_Low)
  )
  long_df$Group2 <- ifelse(long_df$Group == "High", 1, 1.15)
  ylim_sup <- max(long_df$Value) * 1.25
  ylim_inf <- min(long_df$Value)
  
  p <- ggplot(long_df, aes(x = Group2, y = Value, fill = Group)) +
    geom_boxplot(width = 0.12, alpha = 0.3, lwd = 0.5, color = "black") +
    scale_fill_manual(values = c(High = "brown", Low = "steelblue")) +
    coord_cartesian(ylim = c(ylim_inf, ylim_sup), clip = "off") +
    annotate("text", x = 0.95,  y = ylim_sup,        label = paste0(letter, "."), size = 4, fontface = "bold") +
    annotate("text", x = 1.075, y = ylim_sup * 0.90,  label = label,              size = 4, fontface = "bold") +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.title.y = element_blank(), axis.text.y = element_text(size = 12, margin = margin(r = 5)),
      axis.ticks.y = element_line(),
      plot.margin = margin(1, 4, 1, 4)
    )
  
  list(test = test, plot = p)
}

## ---- 6. Run the pipeline for every taxon -----------------------------------
results <- list()

for (sp in species_list) {
  cat("\n=== Processing:", sp$gam_title, "(", sp$letter, ") ===\n")
  
  res <- tryCatch({
    main_df <- read_main_data(sp$main_path)
    bi_df   <- read_biplot_data(sp$bi_path, sp$bi_reader)
    
    detrend <- build_detrend_plot(main_df, sp$n, sp$gam_title, sp$xlim_detrend)
    wilcox  <- build_wilcoxon_boxplot(wilcoxon_data[[sp$letter]], sp$letter)
    
    cat("--- GLS(AR1) summary:", sp$gam_title, "---\n")
    print(summary(detrend$model))
    
    list(
      ts        = build_timeseries_plot(main_df, sp$ts_title, sp$show_xlab),
      bi        = build_biplot(bi_df, sp$bi_x, sp$bi_y, sp$bi_title),
      gam_plot  = detrend$plot,
      gam_model = detrend$model,
      box       = wilcox$plot,
      wilcox    = wilcox$test
    )
  }, error = function(e) {
    warning(sprintf("Skipping '%s' (%s): %s", sp$gam_title, sp$letter, conditionMessage(e)), call. = FALSE)
    NULL
  })
  
  if (!is.null(res)) results[[sp$letter]] <- res
}

# Report any taxa that failed and were skipped, so the combined figures below
# (built only from the taxa that succeeded) make sense.
missing_letters <- setdiff(names(species_list), names(results))
if (length(missing_letters) > 0) {
  cat("\nNOTE: the following taxa were skipped due to a read/processing error and are excluded from the figures below:",
      paste(missing_letters, collapse = ", "), "\n")
}

## ---- 7. Combined figures ----------------------------------------------------
letters_order <- names(results)  # "A".."I", in the order defined above

fig_timeseries <- plot_grid(plotlist = lapply(results, `[[`, "ts"),
                            labels = letters_order, label_size = 12, ncol = 3, align = "hv")

fig_biplots <- plot_grid(plotlist = lapply(results, `[[`, "bi"),
                         labels = letters_order, label_size = 12, ncol = 3, align = "hv")

# Bonus: a combined grid of the 9 GAM-residual correlation plots
# (computed, but never assembled into a figure, in the original script).
fig_gam <- plot_grid(plotlist = lapply(results, `[[`, "gam_plot"),
                     labels = letters_order, label_size = 12, ncol = 3, align = "hv")

# Shared legend (hilltop = brown, lowland = steelblue)
legend_df <- data.frame(Group = factor(c("High", "Low"), levels = c("High", "Low")), x = 1:2, y = 1)
legend_plot <- ggplot(legend_df, aes(x, y, fill = Group)) +
  geom_point(shape = 22, size = 4, stroke = 0.4, colour = "black", alpha = 0) +
  scale_fill_manual(
    values = c(High = "brown", Low = "steelblue"), labels = c("Hills", "Lowlands"), name = NULL,
    guide = guide_legend(override.aes = list(shape = 22, size = 4, stroke = 0.4, alpha = 0.3),
                         nrow = 1, byrow = TRUE)
  ) +
  theme_void() +
  theme(
    legend.position = "top", legend.direction = "horizontal",
    legend.key.size = unit(0.35, "cm"), legend.text = element_text(size = 16),
    legend.spacing.x = unit(0.15, "cm"),
    legend.box.margin = margin(0, 0, 0, 0), plot.margin = margin(0, 0, 0, 0)
  )

top_grid    <- wrap_plots(lapply(results, `[[`, "ts"),  ncol = 3)
bottom_grid <- wrap_plots(lapply(results, `[[`, "box"), ncol = 9)

fig_final <- (legend_plot / top_grid / bottom_grid) +
  plot_layout(heights = c(0.10, 4, 1)) &
  theme(plot.margin = margin(3, 3, 3, 3))

## ---- 8. Display --------------------------------------------------------------
fig_timeseries
fig_biplots
fig_gam
fig_final
