### Mort Modelling###
## Load Libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggpubr)


#### Load Data ####
load("Data/TotalData.RData") #load Data
### Load Functions
source("01_Functions.R")

#### --- Create Geary's C Plot ################################################
ggpubr::ggarrange(GearyCPlot(TotalData = TotalData, 
                             sex="female", LastYearObs = 2017),
                  GearyCPlot(TotalData = TotalData, 
                             sex="male", LastYearObs = 2017),
                  nrow = 2)
ggpubr::annotate_figure(GearyC, top = text_grob("Mean Estimates of Geary's C", 
                                                face = "bold", size = 12))