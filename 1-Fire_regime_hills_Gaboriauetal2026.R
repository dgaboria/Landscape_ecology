# Fire regime (occurrence and biomass burned) reconstruction for each lake (hill sites)
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

library(ggplot2)
library(MASS)
library('ggpubr')
library('vegan')
library(cowplot)
library(BINCOR)
library(pracma)

####################
# File paths
####################

base_dir       <- "D:/PROJETS_POSTDOC/Aiguebelle/Datasets/Files"
gagnon_dir     <- file.path(base_dir, "Lake_Gagnon")
castor_dir     <- file.path(base_dir, "Lake_Castor")
bincor_outfile <- "D:/PROJETS_POSTDOC/Aiguebelle/Datasets/Files"

##########################################################################################################
# 1 - Fire Frequency + confidence interval (CI)
##########################################################################################################

# For each lake: read the dated, significant charcoal peaks (fire events)
# identified by CharAnalysis, then estimate a smoothed fire
# frequency curve (kernel-density fire frequency, kdffreq) with a
# bootstrapped confidence interval.

########################################
#LAKE GAGNON
########################################

Lakes_fires_gagnon<-read.csv(file.path(gagnon_dir, "Lake_fire_frequency.csv"), header = TRUE, sep = ";")
Gagnon<-Lakes_fires_gagnon[, 1]
Gagnon<-Gagnon[!is.na(Gagnon)]
par(mar = rep(2, 4))

# Kernel-density fire frequency estimate
fevent  <- c(Gagnon)
ffGagnon <- kdffreq(fevent, up = -72, lo = 10250, bandwidth = 500, interval = 10,
                    nbboot = 1000, alpha = 0.01, pseudo = TRUE)

mat_ffGagnon <- matrix(cbind(ffGagnon$age, ffGagnon$ff), ncol = 2)

# Base-R diagnostic plot
plot.kdffreq(ffGagnon, ylim = c(0, 0.015), xlim = c(10000, 0), colour = "red",
             bty = "n", ylab = "Fire.millennium", main = "Fire occurrence Gagnon")

GagnonFF1 <- as.data.frame(cbind(ffGagnon$age, ffGagnon$ff, ffGagnon$lo, ffGagnon$up))
colnames(GagnonFF1) <- c("Age", "Freq", "Lo", "Up")

moyGagnon1 <- mean(ffGagnon$ff)
abline(h = moyGagnon1, col = "red", lty = 2)

# Convert frequency to "fires per millennium" for plotting
GagnonFF1$Freq <- GagnonFF1$Freq * 1000
GagnonFF1$Lo   <- GagnonFF1$Lo * 1000
GagnonFF1$Up   <- GagnonFF1$Up * 1000

### PLOT: Lake Gagnon fire frequency ###
# Dated fire events used for the reconstruction (from CharAnalysis output).
# NOTE: this is a hard-coded list of ages (cal. yr BP) — if the source data
# changes, this list needs to be regenerated/re-pasted manually. Consider
# reading these directly from the same CSV used above instead of retyping
# them, to avoid the two ever getting out of sync.

fire_dates2 <- data.frame(Age = c(
  117, 495, 603, 711, 1062, 1197, 1332, 1386, 1602, 1872, 1926, 2088, 2277,
  2358, 2466, 2682, 2817, 2925, 3060, 3276, 3384, 3627, 3897, 4113, 4329,
  4410, 4626, 5058, 5139, 5220, 5301, 5382, 5598, 5652, 5814, 6165, 6300,
  6354, 6462, 6516, 6678, 6786, 7083, 7191, 7299, 7461, 7569, 7812, 7893,
  8244, 8406, 8703, 8892, 9432
))

# Subset of the above that are considered "local" fire events for Gagnon
local_gagnon <- data.frame(Age = c(
  117, 711, 1062, 1197, 1332, 1872, 1926, 2088, 2817, 3897, 4113, 4329,
  4626, 5058, 5301, 5382, 5598, 5652, 5814, 6165, 6678, 6786, 7083, 7569,
  7812, 7893, 8406, 8892, 9432
))

pt1 <- ggplot(data = GagnonFF1) +
  geom_line(aes(x = Age, y = Freq), colour = "black", linewidth = 1.3) +
  geom_point(data = fire_dates2, aes(x = Age, y = 15), shape = 3, size = 1.5, color = "black", stroke = 1) +
  geom_point(data = local_gagnon, aes(x = Age, y = 14), shape = 3, size = 1.5, color = "red", stroke = 1) +
  geom_ribbon(aes(x = Age, ymin = Lo, ymax = Up), linetype = 2, alpha = 0.3, fill = "darkgreen") +
  scale_x_continuous(name = "cal. yr BP", trans = "reverse", breaks = c(10000, 8000, 6000, 4000, 2000, 0)) +
  scale_y_continuous(
    expression("Number of fires per millenium"),
    breaks = c(0, 5, 10, 15),
    limits = c(0, 15)
  ) +
  geom_hline(aes(yintercept = moyGagnon1 * 1000), color = "darkgreen", linetype = "dashed", lwd = 1.5) +
  theme_classic() +
  theme(
    panel.grid.minor.y = element_line(colour = "black"),
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.x = element_text(face = "bold", color = "black", size = 11, angle = 30),
    axis.text.y = element_text(face = "bold", color = "black", size = 11, angle = 30),
    axis.line = element_line(color = "black", size = 0.6, linetype = "solid"),
    axis.title.y = element_text(color = "black", size = 9, face = "bold.italic")
  )
pt1

########################################
#LAKE CASTOR
########################################

Lakes_fires_Castor <- read.csv(file.path(castor_dir, "Lake_fire_frequency.csv"), header = TRUE, sep = ";")
Castor <- Lakes_fires_Castor[, 1]
Castor <- Castor[!is.na(Castor)]
par(mar = rep(2, 4))

fevent  <- c(Castor)
ffCastor <- kdffreq(fevent, up = -72, lo = 10000, bandwidth = 600, interval = 10,
                    nbboot = 1000, alpha = 0.01, pseudo = TRUE)

mat_ffCastor <- matrix(cbind(ffCastor$age, ffCastor$ff), ncol = 2)

plot.kdffreq(ffCastor, ylim = c(0, 0.025), xlim = c(10000, 0), colour = "red",
             bty = "n", ylab = "Fire.millennium", main = "Fire occurrence Castor")

CastorFF <- as.data.frame(cbind(ffCastor$age, ffCastor$ff, ffCastor$lo, ffCastor$up))
colnames(CastorFF) <- c("Age", "Freq", "Lo", "Up")

moyCastor1 <- mean(ffCastor$ff)
abline(h = moyCastor1, col = "red", lty = 2)

CastorFF$Freq <- CastorFF$Freq * 1000
CastorFF$Lo   <- CastorFF$Lo * 1000
CastorFF$Up   <- CastorFF$Up * 1000

### PLOT: Lake Castor fire frequency ###
fire_dates <- data.frame(Age = c(
  93, 318, 438, 513, 978, 1023, 1098, 1308, 1353, 1413, 1698, 1743, 1773,
  1893, 2118, 2283, 2313, 2373, 2463, 2568, 2838, 3018, 3228, 3273, 3348,
  3468, 3603, 3648, 3693, 3753, 3963, 4143, 4488, 4533, 4638, 4773, 4923,
  4998, 5163, 5238, 5388, 5418, 5448, 5478, 5733, 5868, 6018, 6168, 6243,
  6363, 6438, 6543, 6603, 7008, 7578, 7728, 7818, 7848, 8058, 8088, 8328,
  8538, 8808, 9033, 9363
))

local_castor <- data.frame(Age = c(
  93, 1353, 1698, 1743, 1893, 2373, 2568, 2838, 3018, 3348, 3648, 3693,
  3963, 4488, 4533, 4638, 4998, 5163, 5238, 5388, 5478, 5733, 6168, 6603,
  7008, 7818, 7848, 8538
))

pt2 <- ggplot(data = CastorFF) +
  geom_line(aes(x = Age, y = Freq), colour = "black", linewidth = 1.3) +
  geom_ribbon(aes(x = Age, ymin = Lo, ymax = Up), linetype = 2, alpha = 0.3, fill = "darkblue") +
  geom_point(data = fire_dates, aes(x = Age, y = 15), shape = 3, size = 1.5, color = "black", stroke = 1) +
  geom_point(data = local_castor, aes(x = Age, y = 14), shape = 3, size = 1.5, color = "red", stroke = 1) +
  scale_x_continuous(name = "cal. yr BP", trans = "reverse", breaks = c(10000, 8000, 6000, 4000, 2000, 0)) +
  scale_y_continuous(
    expression("Number of fires per millenium"),
    breaks = c(0, 5, 10, 15),
    limits = c(0, 15)
  ) +
  geom_hline(aes(yintercept = moyCastor1 * 1000), color = "darkblue", linetype = "dashed", lwd = 1.5) +
  theme_classic() +
  theme(
    panel.grid.minor.y = element_line(colour = "black"),
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.x = element_text(face = "bold", color = "black", size = 11, angle = 30),
    axis.text.y = element_text(face = "bold", color = "black", size = 11, angle = 30),
    axis.line = element_line(color = "black", size = 0.6, linetype = "solid"),
    axis.title.y = element_text(color = "black", size = 9, face = "bold.italic")
  )
pt2



################################################################################
# 2 - Biomass burned (charcoal accumulation rate, "influx")
################################################################################

# --- Lake Gagnon ---
BBGagnon <- read.csv(file.path(gagnon_dir, "BBGagnon.csv"), header = TRUE, sep = ";", dec = ",")

pt3 <- ggplot(data = BBGagnon) +
  geom_line(aes(x = Cal_BP, y = Influx), colour = "grey", linewidth = 0.1) +
  geom_smooth(aes(x = Cal_BP, y = Influx), method = "loess", span = 0.3, colour = "darkgreen", fill = "darkgreen") +
  scale_x_continuous(trans = "reverse", breaks = c(10000, 8000, 6000, 4000, 2000, 0), name = "cal. yr BP") +
  scale_y_continuous(
    name = expression("Charcoal acc. rate (" * mm^2 * "\u00B7" * cm^{-2} * "\u00B7" * year^{-1} * ")"),
    breaks = c(0, 0.05, 0.1, 0.15),
    limits = c(0, 0.15)
  ) +
  theme_classic() +
  theme(
    panel.grid.major.y = element_line(colour = "black"),
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.x = element_text(face = "bold", color = "black", size = 11, angle = 30),
    axis.text.y = element_text(face = "bold", color = "black", size = 11, angle = 30),
    axis.line = element_line(color = "black", size = 0.6, linetype = "solid"),
    axis.title.x = element_text(color = "black", size = 11),
    axis.title.y = element_text(color = "black", size = 9, face = "bold.italic")
  ) +
  ggtitle("Lake Gagnon") +
  # NOTE: added na.rm = TRUE here for consistency with the Castor plot below
  # (original code omitted it, so this line would silently disappear if
  # BBGagnon$Influx contained any NA values)
  geom_hline(aes(yintercept = mean(Influx, na.rm = TRUE)), color = "darkgreen", linetype = "dashed", lwd = 1.5)
pt3


# --- Lake Castor ---
BBCastor <- read.csv(file.path(castor_dir, "BBCastor.csv"), header = TRUE, sep = ";", dec = ",")

pt4 <- ggplot(data = BBCastor) +
  geom_line(aes(x = Cal_BP, y = Influx), colour = "grey", linewidth = 0.1) +
  geom_smooth(aes(x = Cal_BP, y = Influx), method = "loess", span = 0.3, colour = "darkblue", fill = "darkblue") +
  scale_x_continuous(trans = "reverse", breaks = c(10000, 8000, 6000, 4000, 2000, 0), name = "cal. yr BP") +
  scale_y_continuous(
    name = expression("Charcoal acc. rate (" * mm^2 * "\u00B7" * cm^{-2} * "\u00B7" * year^{-1} * ")"),
    breaks = c(0, 0.05, 0.1, 0.15),
    limits = c(0, 0.15)
  ) +
  theme_classic() +
  theme(
    panel.grid.major.y = element_line(colour = "black"),
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.x = element_text(face = "bold", color = "black", size = 11, angle = 30),
    axis.text.y = element_text(face = "bold", color = "black", size = 11, angle = 30),
    axis.line = element_line(color = "black", size = 0.6, linetype = "solid"),
    axis.title.x = element_text(color = "black", size = 11),
    axis.title.y = element_text(color = "black", size = 9, face = "bold.italic")
  ) +
  ggtitle("Lake Castor") +
  geom_hline(aes(yintercept = mean(Influx, na.rm = TRUE)), color = "darkblue", linetype = "dashed", lwd = 1.5)
pt4

# Combined 2x2 figure: A = Castor biomass burned, B = Gagnon biomass burned,
# C = Castor fire frequency, D = Gagnon fire frequency
combined_fig <- plot_grid(pt4, pt3, pt2, pt1, nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))
combined_fig

################################################################################
# 3 - BINCOR correlations
################################################################################
# Compares Castor vs Gagnon time series (both biomass-burned and fire
# frequency) after binning them onto a common time axis, using the BINCOR
# package.

# --- Biomass-burned (influx) correlation ---
influxcastor <- cbind(BBCastor$Cal_BP[1:672], BBCastor$Influx[1:672])
influxgagnon <- cbind(BBGagnon$Cal_BP[1:374], BBGagnon$Influx[1:374])

test1_allperiod <- bin_cor(influxcastor, influxgagnon, FLAGTAU = 3, ofilename = "bincor_influx_log.txt")
binnedts <- test1_allperiod$Binned_time_series
bin_ts1  <- na.omit(test1_allperiod$Binned_time_series[, 1:2])
bin_ts2  <- na.omit(test1_allperiod$Binned_time_series[, c(1, 3)])
plot_ts(influxcastor, influxgagnon, bin_ts1, bin_ts2, "Castor", "Gagnon",
        colts1 = 1, colts2 = 2, colbints1 = 3, colbints2 = 4, device = "screen")
cor_ts(bin_ts1, bin_ts2, rmltrd = "n", KoCM = "spearman")
plot(bin_ts2[, 2], type = "l")

# --- Fire frequency correlation ---
ffcastor <- cbind(CastorFF$Age, CastorFF$Freq)
ffgagnon <- cbind(GagnonFF1$Age, GagnonFF1$Freq)

test2_allperiod <- bin_cor(ffcastor, ffgagnon, FLAGTAU = 3, ofilename = "bincor_freq_log.txt")
binnedts <- test2_allperiod$Binned_time_series
bin_ts1  <- na.omit(test2_allperiod$Binned_time_series[, 1:2])
bin_ts2  <- na.omit(test2_allperiod$Binned_time_series[, c(1, 3)])
plot_ts(ffcastor, ffgagnon, bin_ts1, bin_ts2, "Castor", "Gagnon",
        colts1 = 1, colts2 = 2, colbints1 = 3, colbints2 = 4, device = "screen")
cor_ts(bin_ts1, bin_ts2, rmltrd = "n", KoCM = "spearman")
plot(bin_ts2[, 2], type = "l")


