#### Load Data ####
load("Data/TotalData.RData") #load Data
### Load Functions
source("01_Functions.R")

##########---- Posterior Predictive Checks -------------########################
##Attention, no such file in Git (too big to upload). 
#Run 04_RunModel first and save result
load(file="../Results/StanAPC_BYM2_F.RData")


## Extract replicated rates from Stan Object FC##
StanFit_APC_BYM2 <- rstan::extract(StanAPC_BYM2_F, permuted=TRUE, 
                                   pars="lambdahat") #rates for PPC's

#Create replicate Data  
InSampleData <- TrainingDataFun(TotalData, 
                                sex="female", 
                                LastYearObs=2014) 

InSampObs <-nrow(InSampleData$Data)

#create 1000 replicas for each observation
DrawsAPC_BYM2 <- matrix(0, nrow = 1000, ncol = InSampObs)

#Draw from poisson distribution with estimated Lambda from Model
for (i in 1:InSampObs) {
  DrawsAPC_BYM2[,i] <- rpois(n=1000, StanFit_APC_BYM2_AllD$lambdahat[,i])
}

#1 Compare replicated with observed Data
#taking only the first 100 for visual inspection (a lot faster)
YRepDensity(Draws = DrawsAPC_BYM2[1:100,], 
            Deaths=InSampleData$Data$Deaths,
            sex="female")


#2.) Compare proportion of zeros
prop_zero <- function(x) mean(x == 0)
prop_zero(InSampleData$Data$Deaths) # check proportion of zeros in y

bayesplot::ppc_stat_grouped(InSampleData$Data$Deaths, DrawsAPC_BYM2, 
                            stat = "prop_zero", 
                            group = InSampleData$Data$AgeID)+
  scale_x_continuous(n.breaks = 3)


