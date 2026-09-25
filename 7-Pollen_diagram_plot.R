# Pollen percentages, and influx stratigraphic diagrams
# UTF-8
# Author: Dorian Gaboriau - dorian.gaboriau@uqat.ca
# Last update : September 2026

# Deleting variables from the environment, and closing any open graphics device

rm(list = ls())
if (!is.null(grDevices::dev.list())) dev.off()

####################
# Import libraries
####################
library(rioja)
library(readxl)
library(tidyr)
library(dplyr)

#Set directory folder
getwd()
setwd('D:/PROJETS_POSTDOC/Aiguebelle/Datasets/Files/Pollen')

###############################################################################
# Pollen percentages
###############################################################################

# upload the raw data set of pollen percentages

data_original <- read_excel("Gagnon_pollen_2026.xlsx", sheet = "Percentage")

# organize the pollen data, reorder taxa and combine some taxa in new categories
# see Table S2

data_reorg <- data_original %>%
  mutate(
    # Combine tree species not listed individually
    Other_trees = rowSums(select(., c(
      "Acer rubrum", "Fraxinus nigra", "Juglans", "Carya ovata", "Tilia americana",
      "Ostrya", "Platanus occidentalis", "Castanea dentata"
    )), na.rm = T),
    
    # Combine shrub species not listed individually
    Other_shrubs = rowSums(select(., c(
     `Myrica gale`, `Salix type`, "Morus rubra","Corylus", "Viburnum"
    )), na.rm = TRUE),
    
    # Combine herb species not listed individually
    Other_herbs = rowSums(select(., c(
      "Rosaceae", "Epilobium", "Rumex", "Agoseris","Ilex", "Menyanthaceae", 
      "Sarcobatus vermiculatus", "Saxifraga", "Caryophyllaceae"
    )), na.rm = TRUE)
  ) %>%
  #  Keep only depths < 370 cm
  filter(Depth < 375) %>%
  # creating the new dataframe with all the species and categories represented
  select(
    Depth, age_best,
    `Picea spp. total`, "Betula spp.", "Pinus strobus", "Pinus cf. banksiana",
    "Pinus indiff. total", "Abies balsamea", "Populus cf. tremuloides",
    "Larix laricina", "Acer saccharum", "Quercus spp.", "Tsuga", "Ulmus",
    "Fraxinus type pennsylanica", "Fagus grandifolia", "Other_trees",
    "Alnus rugosa", "Alnus crispa", "Juniperus/Thuja", "Ericaceae",
    "Other_shrubs", "Poaceae", "Artemesia", "Aster", "Chenopodiaceae",
    "Ambrosia", "Other_herbs","Total terrestre", "Inconnu", "Concentration","Total influx"
  )


# confirming all values are numerical

data_reorg[] <- lapply(data_reorg, as.numeric)

# creating depth and age vectors 
depth <- data_reorg[,1]
age <- data_reorg[,2]

# creating a pollen vector containing the pollen percentages for the selected taxa
# and groups only 
pol <- data_reorg[3:28]

# Constrained cluster analysis to calculate pollen assemblage zones

diss <- dist(sqrt(pol/100)^2)
clust <- chclust(diss, method="coniss")
clust
bstick(clust) # broken stick model suggest 4 significant zones
#we added one additional based on visual observation

# creating vectors for curb and exaggeration curbs colors 
# one color for trees, there are 15 taxa represented

cce=c( rep("firebrick4",15), rep("chartreuse4",5),rep("goldenrod4",6),rep("black", 2))
cc1e=c(rep("firebrick2",15), rep("chartreuse2",5),rep("goldenrod2",6))

# creating a vector for taxa names, with italic format

taxa_names <- c(expression(italic(Picea)~"spp."), expression(italic(Betula)~"spp."),
                expression(italic(Pinus~strobus)),expression(italic(Pinus~banksiana)),
                expression(italic(Pinus)~"indiff."), expression(italic(Abies~balsamea)),
                expression(italic(Populus~tremuloides)), expression(italic(Larix~laricina)),
                expression(italic(Acer~saccharum)), expression(italic(Quercus)~"spp."),
                expression(italic(Tsuga~canadensis)), expression(italic(Ulmus)~"sp."),
                expression(italic(Fraxinus~pennsylvanica)), expression(italic(Fagus~grandifolia)),
                "Other trees", expression(italic(Alnus~rugosa)),
                expression(italic(Alnus~crispa)), "Cupressaceae",
                "Ericaceae","Other shrubs",expression(italic(Poaceae)),expression(italic(Artemisia~"sp.")),
                expression(italic(Aster)~"sp."),expression(italic(Chenopodiaceae)),
                expression(italic(Ambrosia)~"sp."),"Other herbs","Total grains analyzed",
                "Unknown and undetermined pollen (grains)")

# plotting total concentration and influx

# selecting and formation concentration and influx data
conc_influx <- data_reorg %>%
  mutate(Concentration= Concentration/10000)%>%
  select("Concentration","Total influx")

# creating the corresponding labels
label_conc_influx <- c(expression("Concentration\n (grains"~cm^-3~"\u00D7"~10000~")"),
                       expression("Influx (grains"~cm^-2~"\u00D7"~100~")"))

# plotting with age for the y-axis
conc <- strat.plot(conc_influx,yvar = age, scale.percent=TRUE,
           y.rev=TRUE,plot.bar= T, plot.poly = T, exag = T,
           col.poly = "black", col.bar = "grey", col.exag = "darkgrey", 
           xLeft = 0.06, xRight = 0.15,
           srt.xlabel = 45, x.names = label_conc_influx,
           cex.xlabel = 1, cex.ylabel = 1, cex.axis = 1, cex.yaxis = 1,
           x.pc.inc = 5, min.width = 5,
           ylabPos = 3, yTop = 0.8,
           y.tks = c(0,1000,2000,3000,4000,5000,
                     6000,7000,8000,9000,10000),
          ylabel="Age (cal. yr BP)")


# plotting the pollen diagramm next to concentration and influx
perc <- strat.plot(pol,yvar = age, scale.percent=TRUE,
           y.rev=TRUE,plot.bar= T, plot.poly = T, exag = T,
           col.poly = cce, col.bar = "black", col.exag = cc1e, 
           xLeft= 0.15, xRight = 0.8, y.axis = F,
           srt.xlabel = 45, x.names = taxa_names,
           cex.xlabel = 1, cex.ylabel = 2, cex.axis = 1, cex.yaxis = 2,
           x.pc.inc = 5,
           add = T)

# Number of terrestrial grains analyzed 
grains <- data_reorg[29]
# Number of unknown and undetermined grains
unknown <- data_reorg[30]

# add the graph for the number of unknown and undetermined grains analyzed
unk <- strat.plot(unknown, xLeft = 0.79, yvar = age, y.rev = T, plot.bar = T,
           plot.line = F, scale.percent=F,
           xRight = 0.83, y.axis = F, x.pc.inc = 5,
           col.bar="black",lwd.bar = 5, 
           srt.xlabel = 45, x.names = "Unknown and undetermined \npollen (grains)",
           cex.xlabel = 1, cex.ylabel = 1, cex.axis = 1, cex.yaxis = 2,
           add = T)


# add the graph for the number of terrestrial grains analyzed
tot <- strat.plot(grains, xLeft = 0.83, yvar = age, y.rev = T, plot.bar = T,
          lwd.line = 2, clust = clust, col.bar = "black",
           cex.xlabel = 1, cex.ylabel = 2, cex.axis = 1, cex.yaxis = 2,
           x.pc.inc = 20, col.poly =  "grey2",
           srt.xlabel = 45, x.names = "Total grains analyzed",
           xRight = 0.99, y.axis = F, add = T)


# adding the cluster zones on the plot computing by the previous cluster analysis
addClustZone(conc, clust, 5, col="grey25",lwd=1.5, lty = 2)
addClustZone(perc, clust, 5, col="grey25",lwd=1.5, lty = 2)
addClustZone(unk, clust, 5, col="grey25",lwd=1.5, lty = 2)
addClustZone(tot, clust, 5, col="grey25",lwd=1.5, lty = 2)

###############################################################################
# Pollen INFLUX  
###############################################################################
# same principle but with influx data

# upload the raw data set of pollen influx

data_original_influx <- read_excel("Gagnon_pollen_2026.xlsx", sheet = "Influx")

# organize the pollen data, reorder taxa and combine some taxa in new categories

data_reorg_influx <- data_original_influx %>%
  mutate(
    # Combine tree species not listed individually
    Other_trees = rowSums(select(., c(
      "Acer rubrum", "Fraxinus nigra", "Juglans", "Carya ovata", "Tilia americana",
      "Ostrya", "Platanus occidentalis", "Castanea dentata"
    )), na.rm = T),
    
    # Combine shrub species not listed individually
    Other_shrubs = rowSums(select(., c(
      `Myrica gale`, `Salix type`, "Morus rubra","Corylus", "Viburnum"
    )), na.rm = TRUE),
    
    # Combine herb species not listed individually
    Other_herbs = rowSums(select(., c(
      "Rosaceae", "Epilobium", "Rumex", "Agoseris","Ilex", "Menyanthaceae", 
      "Sarcobatus vermiculatus", "Saxifraga", "Caryophyllaceae"
    )), na.rm = TRUE)
  ) %>%
  #  Keep only depths < 370 cm
  filter(Depth < 375) %>%
  select(
    Depth,age_best,
    `Picea spp. total`, "Betula spp.", "Pinus strobus", "Pinus cf. banksiana",
    "Pinus indiff. total", "Abies balsamea", "Populus cf. tremuloides",
    "Larix laricina", "Acer saccharum", "Quercus spp.", "Tsuga", "Ulmus",
    "Fraxinus type pennsylanica", "Fagus grandifolia", "Other_trees",
    "Alnus rugosa", "Alnus crispa", "Juniperus/Thuja", "Ericaceae",
    "Other_shrubs", "Poaceae", "Artemesia", "Aster", "Chenopodiaceae",
    "Ambrosia", "Other_herbs","Total terrestre", "Inconnu", "Concentration","Total influx"
  )

data_reorg_influx[] <- lapply(data_reorg_influx, as.numeric)

age <- data_reorg[,2]

pol_influx <- data_reorg_influx[3:28]
pol_influx <- pol_influx %>%
  mutate(across(everything(), ~ .x * 10))

#Constrained cluster analysis to calculate pollen zones
# we used the same zones created previously

diss <- dist(sqrt(pol/100)^2)
clust <- chclust(diss, method="coniss")
clust
bstick(clust) # broken stick model suggest 4 significant zones

# creating vectors for curb and exageration curbs colours 

cce=c( rep("firebrick4",15), rep("chartreuse4",5),rep("goldenrod4",6),rep("black", 2))
cc1e=c(rep("firebrick2",15), rep("chartreuse2",5),rep("goldenrod2",6))


# creating a vector for taxa names

taxa_names <- c(expression(italic(Picea)~"spp."), expression(italic(Betula)~"spp."),
                expression(italic(Pinus~strobus)),expression(italic(Pinus~banksiana)),
                expression(italic(Pinus)~"indiff."), expression(italic(Abies~balsamea)),
                expression(italic(Populus~tremuloides)), expression(italic(Larix~laricina)),
                expression(italic(Acer~saccharum)), expression(italic(Quercus)~"spp."),
                expression(italic(Tsuga~canadensis)), expression(italic(Ulmus)~"sp."),
                expression(italic(Fraxinus~pennsylvanica)), expression(italic(Fagus~grandifolia)),
                "Other trees", expression(italic(Alnus~rugosa)),
                expression(italic(Alnus~crispa)), "Cupressaceae",
                "Ericaceae","Other shrubs",expression(italic(Poaceae)),expression(italic(Artemisia~"sp.")),
                expression(italic(Aster)~"sp."),expression(italic(Chenopodiaceae)),
                expression(italic(Ambrosia)~"sp."),"Other herbs","Total grains analyzed",
                "Unknown and undetermined pollen (grains)")


# plotting total concentration and influx

conc_influx_influx <- data_reorg_influx %>%
  mutate(Concentration= Concentration/10000)%>%
  select("Concentration","Total influx")

label_conc_influx_influx <- c(expression("Concentration\n (grains"~cm^-3~"\u00D7"~10000~")"),
                       expression("Influx (grains"~cm^-2~"\u00D7"~100~")"))

conc_influx <- strat.plot(conc_influx_influx,yvar = age, scale.percent=TRUE,
                   y.rev=TRUE,plot.bar= T, plot.poly = T, exag = T,
                   col.poly = "black", col.bar = "grey", col.exag = "darkgrey", 
                   xLeft = 0.06, xRight = 0.15,
                   srt.xlabel = 45, x.names = label_conc_influx_influx,
                   cex.xlabel = 1, cex.ylabel = 1, cex.axis = 1, cex.yaxis = 1,
                   x.pc.inc = 5, min.width = 5, add = F,
                   ylabPos = 3, yTop = 0.8, 
                   y.tks = c(0,1000,2000,3000,4000,5000,
                             6000,7000,8000,9000,10000),
                   ylabel="Age (cal. yr BP)")


# plotting the pollen diagramm based on the influx
infl <- strat.plot(pol_influx,yvar = age, scale.percent=TRUE,
                   y.rev=TRUE,plot.bar= T, plot.poly = T, exag = T,
                   col.poly = cce, col.bar = "black", col.exag = cc1e, 
                   xLeft= 0.15, xRight = 0.8, y.axis = F,
                   srt.xlabel = 45, x.names = taxa_names,
                   cex.xlabel = 1, cex.ylabel = 2, cex.axis = 1, cex.yaxis = 2,
                   x.pc.inc = 5,
                   add = T)

# Number of terrestrial grains analyzed 
grains <- data_reorg[29]
# Number of unknown and undetermined grains
unknown <- data_reorg[30]

# add the graph for the number of unknown and undetermined grains analyzed
unk <- strat.plot(unknown, xLeft = 0.79, yvar = age, y.rev = T, plot.bar = T,
                  plot.line = F, scale.percent=F,
                  xRight = 0.83, y.axis = F, x.pc.inc = 5,
                  col.bar="black",lwd.bar = 5, 
                  srt.xlabel = 45, x.names = "Unknown and undetermined \npollen (grains)",
                  cex.xlabel = 1, cex.ylabel = 1, cex.axis = 1, cex.yaxis = 2,
                  add = T)


# add the graph for the number of terrestrial grains analyzed
tot <- strat.plot(grains, xLeft = 0.83, yvar = age, y.rev = T, plot.bar = T,
                  lwd.line = 2, clust = clust, col.bar = "black",
                  cex.xlabel = 1, cex.ylabel = 2, cex.axis = 1, cex.yaxis = 2,
                  x.pc.inc = 20, col.poly =  "grey2",
                  srt.xlabel = 45, x.names = "Total grains analyzed",
                  xRight = 0.99, y.axis = F, add = T)


# adding the cluster zones on the plot 
addClustZone(conc, clust, 5, col="grey25",lwd=1.5, lty = 2)
addClustZone(perc, clust, 5, col="grey25",lwd=1.5, lty = 2)
addClustZone(unk, clust, 5, col="grey25",lwd=1.5, lty = 2)
addClustZone(tot, clust, 5, col="grey25",lwd=1.5, lty = 2)
