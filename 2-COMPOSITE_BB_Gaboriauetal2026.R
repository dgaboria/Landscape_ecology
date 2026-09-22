# Fire regime Composite (for Biomass Burned))
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
library(devtools)
library('ggpubr')
library('vegan')
library(locfit)
library(cowplot)
library(boot)
library(dplyr)
library(readr)
library(nlme)
library(mgcv)

####################
# File paths
####################

base_dir  <- "D:/PROJETS_POSTDOC/Aiguebelle/Datasets/Files"
setwd(base_dir)
############################################################
#Import data (Char acc. rate for lakes on hills and lowlands)
############################################################

# Hilltops
CAS <- read.csv("BBCastor.csv", header = TRUE, sep = ";", dec = ",")        # resolution = 15 years
GAG <- read.csv("BBGagnon.csv", header = TRUE, sep = ";", dec = ",")        # resolution = 27 years
DES <- read.csv("BBDesperiers.csv", header = TRUE, sep = ";", dec = ",")    # resolution = 20 years
PER <- read.csv("BBPerche.csv", header = TRUE, sep = ";", dec = ",")        # resolution = 14 years
LAB <- read.csv("BBLabelle.csv", header = TRUE, sep = ";", dec = ",")       # resolution = 12 years
#Mean resolution = 17,6 years

# Lowlands
FRA <- read.csv("BBFrancis.csv", header = TRUE, sep = ";", dec = ",")      # resolution = 25 years
PAS <- read.csv("BBPasdefondv2.csv", header = TRUE, sep = ";", dec = ",")  # resolution = 24 years
LIL <- read.csv("BBLili.csv", header = TRUE, sep = ";", dec = ",")         # resolution = 28 years
CLO <- read.csv("BBClo.csv", header = TRUE, sep = ";", dec = ",")         # resolution = 29 years
# Mean resolution = 26.5 years

# Plot raw data - hilltops
plot(x = LAB$Age, y = LAB$Influx_Labelle, type = "l", ylim = c(0, 0.5))
lines(x = CAS$Cal_BP, y = CAS$Influx, type = "l", col = "orange")
lines(x = GAG$Cal_BP, y = GAG$Influx, type = "l", col = "red")
lines(x = DES$Age_Cal_BP, y = DES$Influx, type = "l", col = "blue")
lines(x = PER$Cal_BP, y = PER$Influx, type = "l", col = "green")

# Plot raw data - lowlands
plot(x = FRA$Age, y = FRA$Influx_Francis, type = "l", ylim = c(0, 0.5))
lines(x = PAS$Cal_BP, y = PAS$Influx, type = "l", col = "orange")
lines(x = LIL$Cal_BP, y = LIL$Influx, type = "l", col = "blue")
lines(x = CLO$Cal_BP, y = CLO$Influx, type = "l", col = "green")

###HILLS###
### Data for CASTOR LAKE - file corresponding to INPUTS of CharAnalysis
files    <- c("Input_CHARANALYSIS/BBCastor.csv")
metadata <- c("Input_CHARANALYSIS/metadata_castor.csv")
mydata <- pfAddData(files = files, metadata = metadata, type = "CharAnalysis", yrInterp = 1, sep = ";", dec = ".")
TR1   <- pfTransform2(add = mydata, method = c("MinMax", "Box-Cox", "Z-Score"))
COMP2 <- pfCompositeLF2(TR1, hw = 300, nboot = 1000, tarAge = seq(-70, 10000, 10), conf = c(0.05, 0.95)) # hw = half window width
Castor <- as.data.frame(cbind(COMP2$Result$AGE, COMP2$Result$LocFit))
colnames(Castor) <- c("Age", "Locfit")

p1 <- ggplot(Castor, aes(x = Age, y = Locfit)) +
  geom_line(linewidth = 1.2, color = "brown") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.6) +
  scale_x_continuous(trans = "reverse", breaks = c(10000, 8000, 6000, 4000, 2000, 0), limits = c(10000, 0)) +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  theme_classic() +
  labs(
    y = "Charcoal Z-score",
    color = NULL,
    fill = NULL
  ) +
  theme(
    plot.margin = margin(10, 10, 10, 10),
    plot.title = element_text(hjust = 0, vjust = 1.5, face = "bold", size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(color = "black", size = 11, angle = 30),
    axis.text.y = element_text(size = 12, color = "black")
  ) +
  labs(title = "Castor")
p1

### Data for GAGNON LAKE - file corresponding to INPUTS of CharAnalysis
files    <- c("Input_CHARANALYSIS/BBGagnon.csv")
metadata <- c("Input_CHARANALYSIS/metadata_gagnon.csv")
mydata <- pfAddData(files = files, metadata = metadata, type = "CharAnalysis", yrInterp = 1, sep = ";", dec = ".")
TR1   <- pfTransform2(add = mydata, method = c("MinMax", "Box-Cox", "Z-Score"))
COMP2 <- pfCompositeLF2(TR1, hw = 300, nboot = 1000, tarAge = seq(-70, 10000, 10), conf = c(0.05, 0.95))
Gagnon <- as.data.frame(cbind(COMP2$Result$AGE, COMP2$Result$LocFit))
colnames(Gagnon) <- c("Age", "Locfit")

p2 <- ggplot(Gagnon, aes(x = Age, y = Locfit)) +
  geom_line(linewidth = 1.2, color = "brown") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.6) +
  scale_x_continuous(trans = "reverse", breaks = c(10000, 8000, 6000, 4000, 2000, 0), limits = c(10000, 0)) +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  theme_classic() +
  theme(
    plot.margin = margin(10, 10, 10, 10),
    plot.title = element_text(hjust = 0, vjust = 1.5, face = "bold", size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 11, angle = 30),
    axis.text.y = element_text(size = 12, color = "black")
  ) +
  labs(title = "Gagnon")
p2

### Data for PERCHE LAKE - file corresponding to INPUTS of CharAnalysis
files    <- c("Input_CHARANALYSIS/BBPerche.csv")
metadata <- c("Input_CHARANALYSIS/metadata_perche.csv")
mydata <- pfAddData(files = files, metadata = metadata, type = "CharAnalysis", yrInterp = 1, sep = ";", dec = ".")
TR1   <- pfTransform2(add = mydata, method = c("MinMax", "Box-Cox", "Z-Score"))
# NOTE: Perche's tarAge window (6500-9630) is much narrower and starts much
# later than the other hill lakes
COMP2 <- pfCompositeLF2(TR1, hw = 300, nboot = 1000, tarAge = seq(6500, 9630, 10), conf = c(0.05, 0.95))
Perche <- as.data.frame(cbind(COMP2$Result$AGE, COMP2$Result$LocFit))
colnames(Perche) <- c("Age", "Locfit")

p3 <- ggplot(Perche, aes(x = Age, y = Locfit)) +
  geom_line(linewidth = 1.2, color = "brown") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.6) +
  scale_x_continuous(trans = "reverse", breaks = c(10000, 8000, 6000, 4000, 2000, 0), limits = c(10000, 0)) +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  theme_classic() +
  theme(
    plot.margin = margin(10, 10, 10, 10),
    plot.title = element_text(hjust = 0, vjust = 1.5, face = "bold", size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 11, angle = 30),
    axis.text.y = element_text(size = 12, color = "black")
  ) +
  labs(title = "Perché")
p3

### Data for DESPERIERS LAKE - file corresponding to INPUTS of CharAnalysis
files    <- c("Input_CHARANALYSIS/BBDesperiers.csv")
metadata <- c("Input_CHARANALYSIS/metadata_desperiers.csv")
mydata <- pfAddData(files = files, metadata = metadata, type = "CharAnalysis", yrInterp = 1, sep = ";", dec = ".")
TR1   <- pfTransform2(add = mydata, method = c("MinMax", "Box-Cox", "Z-Score"))
COMP2 <- pfCompositeLF2(TR1, hw = 300, nboot = 1000, tarAge = seq(6640, 9920, 10), conf = c(0.05, 0.95))
Desperiers <- as.data.frame(cbind(COMP2$Result$AGE, COMP2$Result$LocFit))
colnames(Desperiers) <- c("Age", "Locfit")

p4 <- ggplot(Desperiers, aes(x = Age, y = Locfit)) +
  geom_line(linewidth = 1.2, color = "brown") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.6) +
  scale_x_continuous(trans = "reverse", breaks = c(10000, 8000, 6000, 4000, 2000, 0), limits = c(10000, 0)) +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  theme_classic() +
  theme(
    plot.margin = margin(10, 10, 10, 10),
    plot.title = element_text(hjust = 0, vjust = 1.5, face = "bold", size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 11, angle = 30),
    axis.text.y = element_text(size = 12, color = "black")
  ) +
  labs(title = "Despériers")
p4

### Data for LABELLE LAKE - file corresponding to INPUTS of CharAnalysis
files    <- c("Input_CHARANALYSIS/BBLabelle.csv")
metadata <- c("Input_CHARANALYSIS/metadata_labelle.csv")
mydata <- pfAddData(files = files, metadata = metadata, type = "CharAnalysis", yrInterp = 1, sep = ";", dec = ".")
TR1   <- pfTransform2(add = mydata, method = c("MinMax", "Box-Cox", "Z-Score"))
COMP2 <- pfCompositeLF2(TR1, hw = 300, nboot = 1000, tarAge = seq(-60, 10350, 10), conf = c(0.05, 0.95))
Labelle <- as.data.frame(cbind(COMP2$Result$AGE, COMP2$Result$LocFit))
colnames(Labelle) <- c("Age", "Locfit")

p5 <- ggplot(Labelle, aes(x = Age, y = Locfit)) +
  geom_line(linewidth = 1.2, color = "brown") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.6) +
  scale_x_continuous(trans = "reverse", breaks = c(10000, 8000, 6000, 4000, 2000, 0), limits = c(10000, 0)) +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  theme_classic() +
  theme(
    plot.margin = margin(10, 10, 10, 10),
    plot.title = element_text(hjust = 0, vjust = 1.5, face = "bold", size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 11, angle = 30),
    axis.text.y = element_text(size = 12, color = "black")
  ) +
  labs(title = "Labelle")
p5

###########################################################################
### LOWLANDS ###
###########################################################################
### Data for LILI LAKE - file corresponding to INPUTS of CharAnalysis
files    <- c("Input_CHARANALYSIS/BBLili.csv")
metadata <- c("Input_CHARANALYSIS/metadata_lili.csv")
mydata <- pfAddData(files = files, metadata = metadata, type = "CharAnalysis", yrInterp = 1, sep = ";", dec = ".")
TR1   <- pfTransform2(add = mydata, method = c("MinMax", "Box-Cox", "Z-Score"))
COMP2 <- pfCompositeLF2(TR1, hw = 300, nboot = 1000, tarAge = seq(-60, 8450, 10), conf = c(0.05, 0.95))
Lili <- as.data.frame(cbind(COMP2$Result$AGE, COMP2$Result$LocFit))
colnames(Lili) <- c("Age", "Locfit")

p6 <- ggplot(Lili, aes(x = Age, y = Locfit)) +
  geom_line(linewidth = 1.2, color = "steelblue") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.6) +
  scale_x_continuous(trans = "reverse", breaks = c(10000, 8000, 6000, 4000, 2000, 0), limits = c(10000, 0)) +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  theme_classic() +
  labs(
    y = "Charcoal Z-score",
    color = NULL,
    fill = NULL
  ) +
  theme(
    plot.margin = margin(10, 10, 10, 10),
    plot.title = element_text(hjust = 0, vjust = 1.5, face = "bold", size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(color = "black", size = 11, angle = 30),
    axis.text.y = element_text(size = 12, color = "black")
  ) +
  labs(title = "Lili")
p6

### Data for PAS-DE-FOND LAKE - file corresponding to INPUTS of CharAnalysis
files    <- c("Input_CHARANALYSIS/BBPasdefondv2.csv")
metadata <- c("Input_CHARANALYSIS/metadata_pasdefond.csv")
mydata <- pfAddData(files = files, metadata = metadata, type = "CharAnalysis", yrInterp = 1, sep = ";", dec = ".")
TR1   <- pfTransform2(add = mydata, method = c("MinMax", "Box-Cox", "Z-Score"))
COMP2 <- pfCompositeLF2(TR1, hw = 300, nboot = 1000, tarAge = seq(-40, 7700, 10), conf = c(0.05, 0.95))
Pasdefond <- as.data.frame(cbind(COMP2$Result$AGE, COMP2$Result$LocFit))
colnames(Pasdefond) <- c("Age", "Locfit")

p7 <- ggplot(Pasdefond, aes(x = Age, y = Locfit)) +
  geom_line(linewidth = 1.2, color = "steelblue") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.6) +
  scale_x_continuous(trans = "reverse", breaks = c(10000, 8000, 6000, 4000, 2000, 0), limits = c(10000, 0)) +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  theme_classic() +
  theme(
    plot.margin = margin(10, 10, 10, 10),
    plot.title = element_text(hjust = 0, vjust = 1.5, face = "bold", size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 11, angle = 30),
    axis.text.y = element_text(size = 12, color = "black")
  ) +
  labs(title = "Pas-de-Fond")
p7

### Data for FRANCIS LAKE - file corresponding to INPUTS of CharAnalysis
files    <- c("Input_CHARANALYSIS/BBFrancis.csv")
metadata <- c("Input_CHARANALYSIS/metadata_francis.csv")
mydata <- pfAddData(files = files, metadata = metadata, type = "CharAnalysis", yrInterp = 1, sep = ";", dec = ".")
TR1   <- pfTransform2(add = mydata, method = c("MinMax", "Box-Cox", "Z-Score"))
COMP2 <- pfCompositeLF2(TR1, hw = 300, nboot = 1000, tarAge = seq(0, 7100, 10), conf = c(0.05, 0.95))
Francis <- as.data.frame(cbind(COMP2$Result$AGE, COMP2$Result$LocFit))
colnames(Francis) <- c("Age", "Locfit")

p8 <- ggplot(Francis, aes(x = Age, y = Locfit)) +
  geom_line(linewidth = 1.2, color = "steelblue") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.6) +
  scale_x_continuous(trans = "reverse", breaks = c(10000, 8000, 6000, 4000, 2000, 0), limits = c(10000, 0)) +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  theme_classic() +
  theme(
    plot.margin = margin(10, 10, 10, 10),
    plot.title = element_text(hjust = 0, vjust = 1.5, face = "bold", size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 11, angle = 30),
    axis.text.y = element_text(size = 12, color = "black")
  ) +
  labs(title = "Francis")
p8

### Data for CLO LAKE - file corresponding to INPUTS of CharAnalysis
files    <- c("Input_CHARANALYSIS/BBClo.csv")
metadata <- c("Input_CHARANALYSIS/metadata_clo.csv")
mydata <- pfAddData(files = files, metadata = metadata, type = "CharAnalysis", yrInterp = 1, sep = ";", dec = ".")
TR1   <- pfTransform2(add = mydata, method = c("MinMax", "Box-Cox", "Z-Score"))
COMP2 <- pfCompositeLF2(TR1, hw = 300, nboot = 1000, tarAge = seq(-70, 8500, 10), conf = c(0.05, 0.95))
Clo <- as.data.frame(cbind(COMP2$Result$AGE, COMP2$Result$LocFit))
colnames(Clo) <- c("Age", "Locfit")

p9 <- ggplot(Clo, aes(x = Age, y = Locfit)) +
  geom_line(linewidth = 1.2, color = "steelblue") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.6) +
  scale_x_continuous(trans = "reverse", breaks = c(10000, 8000, 6000, 4000, 2000, 0), limits = c(10000, 0)) +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  theme_classic() +
  theme(
    plot.margin = margin(10, 10, 10, 10),
    plot.title = element_text(hjust = 0, vjust = 1.5, face = "bold", size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 11, angle = 30),
    axis.text.y = element_text(size = 12, color = "black")
  ) +
  labs(title = "Clo")
p9


#########
#Figure
#########
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9, nrow = 2, ncol = 5)

########################################################################
### COMPOSITE HILLS ###
########################################################################
# Merge the biomass-burned series of each hill lake
BBhills <- merge(Castor, Gagnon, by = "Age", all = TRUE)
BBhills <- merge(BBhills, Perche, by = "Age", all = TRUE)
BBhills <- merge(BBhills, Desperiers, by = "Age", all = TRUE)
BBhills <- merge(BBhills, Labelle, by = "Age", all = TRUE)
colnames(BBhills) <- c("Age", "Castor", "Gagnon", "Perche", "Desperiers", "Labelle")

# Compute the mean biomass-burned value across lakes to get a regional value (BBmean)
BBhills[, 7] <- apply(BBhills[, 2:ncol(BBhills)], 1, mean, na.rm = TRUE)
colnames(BBhills) <- c("Age", "Castor", "Gagnon", "Perche", "Desperiers", "Labelle", "BBmean")
head(BBhills)
# write.csv(BBhills, file = "D:/BB_merge_hills.csv", row.names = FALSE)

plot(BBhills$Age, BBhills$Castor, xlim = rev(range(BBhills$Age)), type = "l", ylim = c(-3, 3), lwd = 2, ylab = "", col = "blue", cex.axis = 1.3)
lines(BBhills$Age, BBhills$Gagnon, col = "forestgreen", lwd = 2)
lines(BBhills$Age, BBhills$Perche, col = "red", lwd = 2)
lines(BBhills$Age, BBhills$Desperiers, col = "black", lwd = 2)
lines(BBhills$Age, BBhills$Labelle, col = "orange", lwd = 2)

### REGIONAL BIOMASS BURNED (REGBB) ###

### HILLS ### WITH BOOTSTRAP

# File with biomass-burned values for each lake and their mean for each year
head(BBhills)
# Keep only the lake columns (not the mean) for the bootstrap
RegBBhills <- BBhills[, c("Castor", "Gagnon", "Perche", "Desperiers", "Labelle")]
# Keep only rows with data from at least two lakes
RegBBhills <- RegBBhills[1:1008, ]
# Confidence interval around RegBB via bootstrap
my.data <- as.data.frame(t(RegBBhills))  # transpose (format required by boot())
R <- 999  # number of bootstrap resamples
l <- length(my.data[1, ])
bootCIhills <- as.data.frame(matrix(ncol = 3, nrow = l))  # empty matrix to hold the output
boot.mean <- function(y, m) {   # statistic function: bootstrap mean
  z <- mean(y[m])
  z
}
for (i in 1:l) {
  y <- my.data[, i]
  y <- na.omit(y)
  boot.obj <- boot(y, boot.mean, R)
  IC.norm <- boot.ci(boot.obj, conf = 0.90, type = "norm")  # change the confidence level here if needed (currently 0.90)
  bootCIhills[i, ] <- IC.norm$normal
}

# Attach the CI to the data frame
# (bootCIhills columns, from boot.ci()$normal, are: confidence level, lower bound, upper bound)
BBBootCIhills <- cbind(BBhills$Age[1:1008], BBhills$BBmean[1:1008], bootCIhills)
colnames(BBBootCIhills) <- c("Age", "Mean", "CI", "Infe", "Supe")
BBmean_smoothhills <- smooth.spline(BBBootCIhills$Age, BBBootCIhills$Mean, spar = 0.4)  # smoothing

########################################################################
### COMPOSITE LOWLANDS ###
########################################################################
# Merge the biomass-burned series of each lowland lake
BBlow <- merge(Lili, Pasdefond, by = "Age", all = TRUE)
BBlow <- merge(BBlow, Francis, by = "Age", all = TRUE)
BBlow <- merge(BBlow, Clo, by = "Age", all = TRUE)
colnames(BBlow) <- c("Age", "Lili", "Pasdefond", "Francis", "Clo")

# Compute the mean biomass-burned value across lakes to get a regional value (BBmean)
BBlow[, 6] <- apply(BBlow[, 2:ncol(BBlow)], 1, mean, na.rm = TRUE)
colnames(BBlow) <- c("Age", "Lili", "Pasdefond", "Francis", "Clo", "BBmean")
head(BBlow)
# write.csv(BBlow, file = "D:/BB_merge_low.csv", row.names = FALSE)

plot(BBlow$Age, BBlow$Lili, xlim = rev(range(BBlow$Age)), type = "l", ylim = c(-3, 3), lwd = 2, ylab = "", col = "blue", cex.axis = 1.3)
lines(BBlow$Age, BBlow$Pasdefond, col = "forestgreen", lwd = 2)
lines(BBlow$Age, BBlow$Francis, col = "red", lwd = 2)
lines(BBlow$Age, BBlow$Clo, col = "black", lwd = 2)

### REGBB ###

### LOWLANDS ### WITH BOOTSTRAP

# File with biomass-burned values for each lake and their mean for each year
head(BBlow)
# Keep only the lake columns (not the mean) for the bootstrap
RegBBlow <- BBlow[, c("Lili", "Pasdefond", "Francis", "Clo")]
# Keep only rows with data from at least two lakes
RegBBlow <- RegBBlow[2:853, ]
# Confidence interval around RegBB via bootstrap
my.data <- as.data.frame(t(RegBBlow))
R <- 999
l <- length(my.data[1, ])
bootCIlow <- as.data.frame(matrix(ncol = 3, nrow = l))
boot.mean <- function(y, m) {
  z <- mean(y[m])
  z
}
for (i in 1:l) {
  y <- my.data[, i]
  y <- na.omit(y)
  boot.obj <- boot(y, boot.mean, R)
  IC.norm <- boot.ci(boot.obj, conf = 0.90, type = "norm")
  bootCIlow[i, ] <- IC.norm$normal
}

# Attach the CI to the data frame
BBBootCIlow <- cbind(BBlow$Age[2:853], BBlow$BBmean[2:853], bootCIlow)
colnames(BBBootCIlow) <- c("Age", "Mean", "CI", "Infe", "Supe")
BBmean_smoothlow <- smooth.spline(BBBootCIlow$Age, BBBootCIlow$Mean, spar = 0.4)  # smoothing

#########
# Figure
#########
graph1 <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 2, ncol = 5)
graph1


############ WILCOXON TEST ############
#########################################
# Holocene: BBhills vs BBlow (-60 to 8450)
#########################################

# --- DATASET 1 Holocene (8450; -60)---
df1 <- data.frame(
  valeur = c(BBhills$BBmean[5:4261],
             BBlow$BBmean[6:4262]),
  groupe = factor(c(
    rep("BBhills", length(5:4261)),
    rep("BBlow",   length(6:4262))
  ))
)

# Extract a single shared legend before removing it from the individual plots
legend_box <- get_legend(
  ggplot(df1, aes(x = groupe, y = valeur, fill = groupe)) +
    geom_boxplot(alpha = 0.4) +
    scale_fill_manual(values = c("BBhills" = "brown", "BBlow" = "steelblue")) +
    theme_classic() +
    theme(legend.position = "right")
)

# Boxplot 1 (keeps the Y axis label)
box1 <- ggplot(df1, aes(x = groupe, y = valeur, fill = groupe)) +
  geom_boxplot(alpha = 0.4) +
  scale_fill_manual(values = c("BBhills" = "brown", "BBlow" = "steelblue")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 9)
  ) +
  labs(y = "BBmean", title = "Holocene (8450 -60)")
box1
w1 <- wilcox.test(BBhills$BBmean[5:4261], BBlow$BBmean[6:4262], paired = FALSE)
w1

#########################################
# Late Holocene (-60 to 1000)
#########################################

df2 <- data.frame(
  valeur = c(BBhills$BBmean[5:537],
             BBlow$BBmean[6:538]),
  groupe = factor(c(
    rep("BBhills", length(5:537)),
    rep("BBlow",   length(6:538))
  ))
)

box2 <- ggplot(df2, aes(x = groupe, y = valeur, fill = groupe)) +
  geom_boxplot(alpha = 0.4) +
  scale_fill_manual(values = c("BBhills" = "brown", "BBlow" = "steelblue")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 9)
  ) +
  labs(title = "1000 -60")
box2
w2 <- wilcox.test(BBhills$BBmean[5:537], BBlow$BBmean[6:538], paired = FALSE)
w2

#########################################
# Holocene (8450 to 1000)
#########################################

df3 <- data.frame(
  valeur = c(BBhills$BBmean[537:4261],
             BBlow$BBmean[538:4262]),
  groupe = factor(c(
    rep("BBhills", length(537:4261)),
    rep("BBlow",   length(538:4262))
  ))
)

box3 <- ggplot(df3, aes(x = groupe, y = valeur, fill = groupe)) +
  geom_boxplot(alpha = 0.4) +
  scale_fill_manual(values = c("BBhills" = "brown", "BBlow" = "steelblue")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 9)
  ) +
  labs(title = "8450-1000")
box3
w3 <- wilcox.test(BBhills$BBmean[537:4261], BBlow$BBmean[538:4262], paired = FALSE)
w3

#########################################
# Mid-Holocene (8200 to 4200)
#########################################

df4 <- data.frame(
  valeur = c(BBhills$BBmean[428:828],
             BBlow$BBmean[428:828]),
  groupe = factor(c(
    rep("BBhills", length(428:828)),
    rep("BBlow",   length(428:828))
  ))
)
max4 <- max(df4$valeur, na.rm = TRUE)

box4 <- ggplot(df4, aes(x = groupe, y = valeur, fill = groupe)) +
  geom_boxplot(alpha = 0.4) +
  annotate("text", x = 1.5, y = 1, label = "p < 0.001", size = 4.2, fontface = "bold") +
  scale_fill_manual(values = c("BBhills" = "brown", "BBlow" = "steelblue")) +
  scale_y_continuous(limits = c(-1, 1)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 12)
  ) +
  labs(title = "Mid-Holocene (8,200-4,200 cal. yr BP)")
box4
w4 <- wilcox.test(BBhills$BBmean[428:828], BBlow$BBmean[428:828], paired = FALSE)
w4

#########################################
# Late Holocene (4200 to -60)
#########################################

df5 <- data.frame(
  valeur = c(BBhills$BBmean[2:428],
             BBlow$BBmean[2:428]),
  groupe = factor(c(
    rep("BBhills", length(2:428)),
    rep("BBlow",   length(2:428))
  ))
)

max5 <- max(df5$valeur, na.rm = TRUE)

box5 <- ggplot(df5, aes(x = groupe, y = valeur, fill = groupe)) +
  geom_boxplot(alpha = 0.4) +
  annotate("text", x = 1.5, y = 1, label = "p < 0.001", size = 4.2, fontface = "bold") +
  scale_fill_manual(values = c("BBhills" = "brown", "BBlow" = "steelblue")) +
  scale_y_continuous(limits = c(-1, 1)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 12)
  ) +
  labs(title = "Late Holocene (4,200-0 cal. yr BP)")
box5
w5 <- wilcox.test(BBhills$BBmean[2:428], BBlow$BBmean[2:428], paired = FALSE)
w5

##################################################################################
##################################################################################
# Plot RegBB and CI
BBmean_infhill <- smooth.spline(BBBootCIhills$Age, BBBootCIhills$Infe, spar = 0.4)
BBmean_suphill <- smooth.spline(BBBootCIhills$Age, BBBootCIhills$Supe, spar = 0.4)
plot(BBmean_infhill, type = "l", xlim = c(10000, 0), ylim = c(-2, 2), col = alpha("brown", .2), axes = FALSE, font.axis = 2, xlab = " ", ylab = " ", cex.axis = 0.8, lty = 1, bty = "n", lwd = 2, xaxt = "n")
lines(BBmean_suphill, col = alpha("brown", .2), lwd = 2, lty = 1)
polygon(c(BBmean_smoothhills$x, rev(BBmean_smoothhills$x)), c(BBmean_suphill$y, rev(BBmean_infhill$y)), col = alpha("brown", .2), border = NA)
linemeanBBhills <- mean(BBmean_smoothhills$y)
lines(BBmean_smoothhills$x, BBmean_smoothhills$y, col = "brown", lwd = 2)
title(ylab = "RegBB (no unit)")
axis(2, at = NULL, labels = TRUE, lty = 1, lwd = 2, font.axis = 2, cex.axis = 1.5, las = 1)
axis(3, at = c(0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000), lty = 1, lwd = 2, font.axis = 2, cex.axis = 1.5)

library(scales)  # for alpha()

# --- Smoothed data (recomputed for the ggplot version below) ---
BBmean_infhill     <- smooth.spline(BBBootCIhills$Age, BBBootCIhills$Infe, spar = 0.4)
BBmean_suphill     <- smooth.spline(BBBootCIhills$Age, BBBootCIhills$Supe, spar = 0.4)
BBmean_smoothhills <- smooth.spline(BBBootCIhills$Age, BBBootCIhills$Mean, spar = 0.4)
BBmean_inflow      <- smooth.spline(BBBootCIlow$Age, BBBootCIlow$Infe, spar = 0.4)
BBmean_suplow      <- smooth.spline(BBBootCIlow$Age, BBBootCIlow$Supe, spar = 0.4)
BBmean_smoothlow   <- smooth.spline(BBBootCIlow$Age, BBBootCIlow$Mean, spar = 0.4)

# --- Convert to data frames for ggplot ---
df_hill <- data.frame(
  Age = BBmean_smoothhills$x,
  Mean = BBmean_smoothhills$y,
  Lower = approx(BBmean_infhill$x, BBmean_infhill$y, xout = BBmean_smoothhills$x)$y,
  Upper = approx(BBmean_suphill$x, BBmean_suphill$y, xout = BBmean_smoothhills$x)$y
)

df_low <- data.frame(
  Age = BBmean_smoothlow$x,
  Mean = BBmean_smoothlow$y,
  Lower = approx(BBmean_inflow$x, BBmean_inflow$y, xout = BBmean_smoothlow$x)$y,
  Upper = approx(BBmean_suplow$x, BBmean_suplow$y, xout = BBmean_smoothlow$x)$y
)

# --- Add a series identifier ---
df_hill$zone <- "Hill"
df_low$zone  <- "Low"

# --- Combine ---
df_all <- rbind(df_hill, df_low)

# --- Colors ---
cols_line <- c(Hill = "brown",      Low = "steelblue")
cols_fill <- c(Hill = alpha("brown", 0.20),
               Low  = alpha("steelblue", 0.20))

# --- Plot ---
p <- ggplot(df_all, aes(x = Age, y = Mean, color = zone, fill = zone)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), color = NA) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_reverse(limits = c(10000, 0),
                  breaks = seq(0, 10000, by = 1000),
                  expand = c(0, 0)) +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  scale_color_manual(
    values = c(Hill = "brown", Low = "steelblue"),
    labels = c(Hill = "Hills", Low = "Lowlands")
  ) +
  scale_fill_manual(
    values = c(Hill = alpha("brown", 0.2), Low = alpha("steelblue", 0.2)),
    labels = c(Hill = "Hills", Low = "Lowlands")
  ) +
  labs(
    x = "Age (cal. yr BP)",
    y = "Charcoal Z-score",
    color = NULL,
    fill = NULL
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 13, face = "plain"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text    = element_text(size = 12, color = "black"),
    axis.ticks   = element_line(color = "black"),
    axis.line    = element_line(linewidth = 1),
    plot.margin  = margin(10, 10, 10, 10),
    legend.position = "top",
    legend.text = element_text(size = 12)
  )

p

# Figure
plot_grid(graph1, p, nrow = 2, ncol = 1)

# 1) Stack box4 and box5
right_boxes <- plot_grid(
  box4, box5,
  ncol = 1,
  rel_heights = c(1, 1),
  align = "v"
)

# 2) Bottom row: p + (box4 + box5)
bottom_row <- plot_grid(
  p, right_boxes,
  ncol = 2,
  rel_widths = c(3, 1),   # give p more width
  align = "h"
)

# 3) Final figure: graph1 spanning the full width
final_fig <- plot_grid(
  graph1,
  bottom_row,
  ncol = 1,
  rel_heights = c(1, 1.3)  # adjust to taste
)

final_fig


#################################################
# Correlations: GLS + AR1 & Spearman
#################################################

library(tidyverse)

#################################################
# 1) HILLS
#################################################

df <- BBhills
df_bins <- df %>%
  filter(Age <= 8500) %>%
  mutate(bin = cut(Age,
                   breaks = seq(0, 8500, by = 500),
                   include.lowest = TRUE,
                   right = FALSE))

table_bbmean <- df_bins %>%
  group_by(bin) %>%
  summarise(BBmean = mean(BBmean, na.rm = TRUE))

table_bbmean$bin <- seq(0, 8500, by = 500)

#################################################
# 2) LOWLANDS
#################################################

df2 <- BBlow

df_bins <- df2 %>%
  filter(Age <= 8500) %>%
  mutate(bin = cut(Age,
                   breaks = seq(0, 8500, by = 500),
                   include.lowest = TRUE,
                   right = FALSE))

table_bbmean2 <- df_bins %>%
  group_by(bin) %>%
  summarise(BBmean = mean(BBmean, na.rm = TRUE))

table_bbmean2$bin <- seq(0, 8500, by = 500)

#################################################
# 3) MERGE TABLES
#################################################

df_bi <- table_bbmean %>%
  rename(BBmean_Hills = BBmean) %>%
  inner_join(
    table_bbmean2 %>% rename(BBmean_Lowlands = BBmean),
    by = "bin"
  )

#################################################
# 4) FULL BIPLOT (0-8500)
#################################################

biplotbb1 <- ggplot(df_bi, aes(x = BBmean_Hills, y = BBmean_Lowlands)) +
  geom_point(size = 3, color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
  geom_text(aes(label = bin), nudge_x = 0.05, nudge_y = 0.05, size = 3) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(size = 9)) +
  labs(
    x = "BBmean Hills",
    y = "BBmean Lowlands",
    title = "BBmean (0-8,500 cal. yr BP)"
  )

biplotbb1

#################################################
# 1) COMPUTE CONFIDENCE INTERVALS FOR BBmean (HILLS + LOWLANDS)
#################################################

# --- HILLS ---
df_hills <- BBhills %>%   # <-- your BBhills table
  filter(Age <= 8500) %>%
  mutate(bin = cut(Age, breaks = seq(0, 8500, 500),
                   include.lowest = TRUE, right = FALSE))

IC_hills <- df_hills %>%
  group_by(bin) %>%
  summarise(
    mean  = mean(BBmean, na.rm = TRUE),
    lower = quantile(BBmean, 0.025, na.rm = TRUE),
    upper = quantile(BBmean, 0.975, na.rm = TRUE)
  )

# --- LOWLANDS ---
df_low <- BBlow %>%       # <-- your BBlow table
  filter(Age <= 8500) %>%
  mutate(bin = cut(Age, breaks = seq(0, 8500, 500),
                   include.lowest = TRUE, right = FALSE))

IC_low <- df_low %>%
  group_by(bin) %>%
  summarise(
    mean  = mean(BBmean, na.rm = TRUE),
    lower = quantile(BBmean, 0.025, na.rm = TRUE),
    upper = quantile(BBmean, 0.975, na.rm = TRUE)
  )

#################################################
# 2) MERGE + RENAME
#################################################

df_bi_IC <- IC_hills %>%
  rename(BBmean_Hills = mean,
         Hills_low = lower,
         Hills_up  = upper) %>%
  inner_join(
    IC_low %>%
      rename(BBmean_Lowlands = mean,
             Low_low = lower,
             Low_up  = upper),
    by = "bin"
  )

#################################################
# 3) ADD AGE CLASS (BLUE / RED)
#################################################

df_bi_IC <- df_bi_IC %>%
  mutate(
    bin_start = parse_number(as.character(bin)),
    AgeClass = case_when(
      bin_start < 4200 ~ "4200-0",       # blue
      TRUE ~ "8200-4200"                 # red
    )
  )

#################################################
# 4) BBmean BIPLOT WITH CI + COLORS
#################################################

biplotbb1 <- ggplot(df_bi_IC,
                    aes(x = BBmean_Hills,
                        y = BBmean_Lowlands,
                        color = AgeClass)) +
  
  # Horizontal CIs (Hills)
  geom_errorbarh(aes(xmin = Hills_low, xmax = Hills_up),
                 height = 0.05) +
  
  # Vertical CIs (Lowlands)
  geom_errorbar(aes(ymin = Low_low, ymax = Low_up),
                width = 0.05) +
  
  # Points
  geom_point(size = 3) +
  
  # Regression (optional)
  # geom_smooth(method = "lm", se = FALSE,
  #             color = "black", linewidth = 1) +
  
  # Blue / red colors
  scale_color_manual(
    values = c("8200-4200" = "red",
               "4200-0"   = "blue"),
    name = " "
  ) +
  
  theme_minimal(base_size = 14) +
  labs(
    x = "Hills",
    y = "Lowlands",
    title = "RegBB"
  )

biplotbb1

############################################
# 2) BIPLOT 0-4000
############################################
# Filter for 0-4000 cal BP only
df_bi_0_4000 <- df_bi %>% filter(bin <= 4000)

biplotbb2 <- ggplot(df_bi_0_4000, aes(x = BBmean_Hills, y = BBmean_Lowlands)) +
  geom_point(size = 3, color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
  geom_text(aes(label = bin), nudge_x = 0.05, nudge_y = 0.05, size = 3) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(size = 9)) +
  labs(
    x = "BBmean Hills",
    y = "BBmean Lowlands",
    title = "BBmean (0-4,000 cal. yr BP)"
  )

biplotbb2

############################################
# 3) BIPLOT 4000-8000
############################################
# Filter for 4000-8000 cal BP only
df_bi_4000_8000 <- df_bi %>% filter(bin >= 4000)

biplotbb3 <- ggplot(df_bi_4000_8000, aes(x = BBmean_Hills, y = BBmean_Lowlands)) +
  geom_point(size = 3, color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
  geom_text(aes(label = bin), nudge_x = 0.05, nudge_y = 0.05, size = 3) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(size = 9)) +
  labs(
    x = "BBmean Hills",
    y = "BBmean Lowlands",
    title = "BBmean (4000-8000 cal. yr BP)"
  )

biplotbb3

#############################################################
# Detrend Hills and Lowlands with a GAM
#############################################################

# Holocene: 8,450-0

# 1. Detrend
fit_hills <- gam(BBhills$BBmean ~ s(BBhills$Age, k = 50))
fit_low   <- gam(BBlow$BBmean ~ s(BBlow$Age, k = 50))

res_hills <- resid(fit_hills)
res_low   <- resid(fit_low)

# 2. Prepare the data

df <- data.frame(
  age = BBhills$Age[2:853],
  hills = res_hills[2:853],
  lowlands = res_low[2:853]
)

# 3. GLS + AR1
model_det <- gls(hills ~ lowlands, correlation = corAR1(form = ~ age), data = df)
summary(model_det)

# 4. Spearman correlation
cor.test(df$hills, df$lowlands, method = "spearman")

# 95% CI for Hills
pred_hills <- predict(fit_hills, se.fit = TRUE)
hills_df <- data.frame(
  age = BBhills$Age,
  fit = pred_hills$fit,
  se = pred_hills$se.fit,
  lower = pred_hills$fit - 1.96 * pred_hills$se.fit,
  upper = pred_hills$fit + 1.96 * pred_hills$se.fit
)

# 95% CI for Lowlands
pred_low <- predict(fit_low, se.fit = TRUE)
low_df <- data.frame(
  age = BBlow$Age,
  fit = pred_low$fit,
  se = pred_low$se.fit,
  lower = pred_low$fit - 1.96 * pred_low$se.fit,
  upper = pred_low$fit + 1.96 * pred_low$se.fit
)

df_plot <- data.frame(
  age = df$age,
  hills = df$hills,
  lowlands = df$lowlands,
  hills_se = hills_df$se[2:853],
  low_se = low_df$se[2:853]
)

A_plot <- ggplot(df_plot, aes(x = lowlands, y = hills)) +
  geom_errorbar(aes(ymin = hills - 1.96 * hills_se,
                    ymax = hills + 1.96 * hills_se),
                width = 0, alpha = 0.3) +
  geom_errorbarh(aes(xmin = lowlands - 1.96 * low_se,
                     xmax = lowlands + 1.96 * low_se),
                 height = 0, alpha = 0.3) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_minimal() +
  labs(
    x = "Lowlands (residuals from GAM)",
    y = "Hills (residuals from GAM)",
    title = "RegBB - 8,450 - 0 cal. BP"
  )

A_plot

#######################################
# 2. Prepare the data (window: 1,000-0 cal. BP)
df <- data.frame(
  age = BBhills$Age[2:108],
  hills = res_hills[2:108],
  lowlands = res_low[2:108]
)

# 3. GLS + AR1
model_det <- gls(hills ~ lowlands, correlation = corAR1(form = ~ age), data = df)
summary(model_det)

# 4. Spearman correlation
cor.test(df$hills, df$lowlands, method = "spearman")

# 95% CI for Hills
pred_hills <- predict(fit_hills, se.fit = TRUE)
hills_df <- data.frame(
  age = BBhills$Age,
  fit = pred_hills$fit,
  se = pred_hills$se.fit,
  lower = pred_hills$fit - 1.96 * pred_hills$se.fit,
  upper = pred_hills$fit + 1.96 * pred_hills$se.fit
)

# 95% CI for Lowlands
pred_low <- predict(fit_low, se.fit = TRUE)
low_df <- data.frame(
  age = BBlow$Age,
  fit = pred_low$fit,
  se = pred_low$se.fit,
  lower = pred_low$fit - 1.96 * pred_low$se.fit,
  upper = pred_low$fit + 1.96 * pred_low$se.fit
)

df_plot <- data.frame(
  age = df$age,
  hills = df$hills,
  lowlands = df$lowlands,
  hills_se = hills_df$se[2:108],
  low_se = low_df$se[2:108]
)

B_plot <- ggplot(df_plot, aes(x = lowlands, y = hills)) +
  geom_errorbar(aes(ymin = hills - 1.96 * hills_se,
                    ymax = hills + 1.96 * hills_se),
                width = 0, alpha = 0.3) +
  geom_errorbarh(aes(xmin = lowlands - 1.96 * low_se,
                     xmax = lowlands + 1.96 * low_se),
                 height = 0, alpha = 0.3) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_minimal() +
  labs(
    x = "Lowlands (residuals from GAM)",
    y = "Hills (residuals from GAM)",
    title = "RegBB - 1,000 - 0 cal. BP"
  )
B_plot
#######################################
# 3. Prepare the data (window: 8,450-1,000 cal. BP)
df <- data.frame(
  age = BBhills$Age[108:853],
  hills = res_hills[108:853],
  lowlands = res_low[108:853]
)

# 3. GLS + AR1
model_det <- gls(hills ~ lowlands, correlation = corAR1(form = ~ age), data = df)
summary(model_det)

# 4. Spearman correlation
cor.test(df$hills, df$lowlands, method = "spearman")

# 95% CI for Hills
pred_hills <- predict(fit_hills, se.fit = TRUE)
hills_df <- data.frame(
  age = BBhills$Age,
  fit = pred_hills$fit,
  se = pred_hills$se.fit,
  lower = pred_hills$fit - 1.96 * pred_hills$se.fit,
  upper = pred_hills$fit + 1.96 * pred_hills$se.fit
)

# 95% CI for Lowlands
pred_low <- predict(fit_low, se.fit = TRUE)
low_df <- data.frame(
  age = BBlow$Age,
  fit = pred_low$fit,
  se = pred_low$se.fit,
  lower = pred_low$fit - 1.96 * pred_low$se.fit,
  upper = pred_low$fit + 1.96 * pred_low$se.fit
)

df_plot <- data.frame(
  age = df$age,
  hills = df$hills,
  lowlands = df$lowlands,
  hills_se = hills_df$se[108:853],
  low_se = low_df$se[108:853]
)

C_plot <- ggplot(df_plot, aes(x = lowlands, y = hills)) +
  geom_errorbar(aes(ymin = hills - 1.96 * hills_se,
                    ymax = hills + 1.96 * hills_se),
                width = 0, alpha = 0.3) +
  geom_errorbarh(aes(xmin = lowlands - 1.96 * low_se,
                     xmax = lowlands + 1.96 * low_se),
                 height = 0, alpha = 0.3) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_minimal() +
  labs(
    x = "Lowlands (residuals from GAM)",
    y = "Hills (residuals from GAM)",
    title = "RegBB - 8,450 - 1,000 cal. BP"
  )

C_plot
#######################################
# 4. Prepare the data (window: 8,200-4,200 cal. BP)
df <- data.frame(
  age = BBhills$Age[428:828],
  hills = res_hills[428:828],
  lowlands = res_low[428:828]
)

# 3. GLS + AR1
model_det <- gls(hills ~ lowlands, correlation = corAR1(form = ~ age), data = df)
summary(model_det)

# 4. Spearman correlation
cor.test(df$hills, df$lowlands, method = "spearman")

# 95% CI for Hills
pred_hills <- predict(fit_hills, se.fit = TRUE)
hills_df <- data.frame(
  age = BBhills$Age,
  fit = pred_hills$fit,
  se = pred_hills$se.fit,
  lower = pred_hills$fit - 1.96 * pred_hills$se.fit,
  upper = pred_hills$fit + 1.96 * pred_hills$se.fit
)

# 95% CI for Lowlands
pred_low <- predict(fit_low, se.fit = TRUE)
low_df <- data.frame(
  age = BBlow$Age,
  fit = pred_low$fit,
  se = pred_low$se.fit,
  lower = pred_low$fit - 1.96 * pred_low$se.fit,
  upper = pred_low$fit + 1.96 * pred_low$se.fit
)

df_plot <- data.frame(
  age = df$age,
  hills = df$hills,
  lowlands = df$lowlands,
  hills_se = hills_df$se[428:828],
  low_se = low_df$se[428:828]
)

D_plot <- ggplot(df_plot, aes(x = lowlands, y = hills)) +
  geom_errorbar(aes(ymin = hills - 1.96 * hills_se,
                    ymax = hills + 1.96 * hills_se),
                width = 0, alpha = 0.3) +
  geom_errorbarh(aes(xmin = lowlands - 1.96 * low_se,
                     xmax = lowlands + 1.96 * low_se),
                 height = 0, alpha = 0.3) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_minimal() +
  labs(
    x = "Lowlands (residuals from GAM)",
    y = "Hills (residuals from GAM)",
    title = "RegBB - 8,200 - 4,200 cal. BP"
  )

D_plot
#######################################
# 5. Prepare the data (window: 4,200-0 cal. BP)
df <- data.frame(
  age = BBhills$Age[2:428],
  hills = res_hills[2:428],
  lowlands = res_low[2:428]
)

# 3. GLS + AR1
model_det <- gls(hills ~ lowlands, correlation = corAR1(form = ~ age), data = df)
summary(model_det)

# 4. Spearman correlation
cor.test(df$hills, df$lowlands, method = "spearman")

# 95% CI for Hills
pred_hills <- predict(fit_hills, se.fit = TRUE)
hills_df <- data.frame(
  age = BBhills$Age,
  fit = pred_hills$fit,
  se = pred_hills$se.fit,
  lower = pred_hills$fit - 1.96 * pred_hills$se.fit,
  upper = pred_hills$fit + 1.96 * pred_hills$se.fit
)

# 95% CI for Lowlands
pred_low <- predict(fit_low, se.fit = TRUE)
low_df <- data.frame(
  age = BBlow$Age,
  fit = pred_low$fit,
  se = pred_low$se.fit,
  lower = pred_low$fit - 1.96 * pred_low$se.fit,
  upper = pred_low$fit + 1.96 * pred_low$se.fit
)

df_plot <- data.frame(
  age = df$age,
  hills = df$hills,
  lowlands = df$lowlands,
  hills_se = hills_df$se[2:428],
  low_se = low_df$se[2:428]
)

E_plot <- ggplot(df_plot, aes(x = lowlands, y = hills)) +
  geom_errorbar(aes(ymin = hills - 1.96 * hills_se,
                    ymax = hills + 1.96 * hills_se),
                width = 0, alpha = 0.3) +
  geom_errorbarh(aes(xmin = lowlands - 1.96 * low_se,
                     xmax = lowlands + 1.96 * low_se),
                 height = 0, alpha = 0.3) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_minimal() +
  labs(
    x = "Lowlands (residuals from GAM)",
    y = "Hills (residuals from GAM)",
    title = "RegBB - 4,200 - 0 cal. BP"
  )

E_plot

# Optional: harmonize title sizes
A_plot <- A_plot + theme(plot.title = element_text(size = 10))
B_plot <- B_plot + theme(plot.title = element_text(size = 10))
C_plot <- C_plot + theme(plot.title = element_text(size = 10))
D_plot <- D_plot + theme(plot.title = element_text(size = 10))
E_plot <- E_plot + theme(plot.title = element_text(size = 10))

# Combine into a 2x3 grid
combined_plot <- plot_grid(
  A_plot, B_plot, C_plot,
  D_plot, E_plot, NULL,
  labels = c("A", "B", "C", "D", "E", ""),
  label_size = 12,
  ncol = 3,
  align = "hv"
)

combined_plot
