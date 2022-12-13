### Mort Modelling###
## Load Libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(rstan)

#### Load Data ####
load("Data/TotalData.RData") #load Data
### Load Functions
source("01_Functions.R")

####---- Estimate Stan Model ##################################################
## Example with APC_BYM2 Model

#1. Get Data for stan function
# For other model change Parameter in Function accordingly
DataAPC_BYM2_Stan <- StanData(Data=TotalData, 
                              LastYearObs = 2016, #last Year of InSample Data 
                              AdjMatType = 3, #3 being binary indicators (1 & 2 currently not in use)
                              sex="female", #Sex (weiblich = female, männlich = male)
                              ModelType = "APC", #APC or RH
                              RegionType = "BYM2", #BYM2 or SAR
                              Cohort=TRUE, #Cohort Index True or False 
                              #e.g. for APC Family (TRUE = APC, FALSE = AP)
                              TFor = 1) #Forecast Horzion (how many years to forecast)

#2.) run Stan Code
#See Appendix for Parameters of HMC Modelling For each Model
##### !!!!!ATTENTION: TAKES A LONG TIME (3-5 hours) !!!!!!!!
options(mc.cores = parallel::detectCores())
StanAPC_BYM2_F <- rstan::stan(file.path(getwd(),"StanCode/APC_BYM2.stan"),
                              data=DataAPC_BYM2_Stan,
                              chains = 4, iter=4000,warmup=2000, 
                              save_warmup=FALSE, thin=4,
                              control = list(adapt_delta = 0.81))
save(StanAPC_BYM2_F, 
     file="../Results/StanAPC_BYM2_F.RData")