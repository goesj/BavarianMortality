### Mort Modelling###
## Load Libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(rstan)
# if youre having trouble installing rstan see
#https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

#### Load Data ####
load("Data/TotalData.RData") #load Data
### Load Functions
source("01_Functions.R")

####---- Estimate Stan Model ##################################################
## Example with APC_BYM2 Model

#1. Get Data for stan function
# For other model change Parameter in Function accordingly
DataAPC_BYM2_Stan <- StanData(Data=TotalData, 
                              LastYearObs = 2014, #last Year of InSample Data 
                              sex="female", #Sex 
                              ModelType = "APC", #APC or RH
                              TFor = 3) #Forecast Horzion 

#2.) run Stan Code
#See Appendix for Parameters of HMC Modelling For each Model
##### !!!!!ATTENTION: TAKES A LONG TIME (3-5 hours) !!!!!!!!
options(mc.cores = parallel::detectCores())
StanAPC_BYM2_F <- rstan::stan(file.path(getwd(),"StanCode/APC_BYM2.stan"),
                              data=DataAPC_BYM2_Stan,
                              chains = 4, iter=4000,warmup=2000, 
                              save_warmup=FALSE, thin=4,
                              control = list(adapt_delta = 0.81))
#cannot be uploaded into Git, file too big (1.8 GB)
save(StanAPC_BYM2_F_1517, 
     file="../Results/StanAPC_BYM2_F_1517.RData")
