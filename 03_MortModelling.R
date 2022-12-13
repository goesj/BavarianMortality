### Mort Modelling###
## Load Libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(rstan, ggpubr)


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
# ATTENTION: TAKES A LONG TIME (3-5 hours)
options(mc.cores = parallel::detectCores())
StanAPC_BYM2_F <- rstan::stan(file.path(getwd(),"StanCode/APC_BYM2.stan"),
                            data=DataAPC_BYM2_Stan,
                            chains = 4, iter=4000,warmup=2000, 
                            save_warmup=FALSE, thin=4,
                            control = list(adapt_delta = 0.81))
save(StanAPC_BYM2_F, 
     file="../Results/StanAPC_BYM2_F.RData")

##########---- Posterior Predictive Checks -------------########################
load(file="../Results/StanAPC_BYM2_F.RData")
## Extract Stan FC##
StanFit_APC_BYM2 <- rstan::extract(StanAPC_BYM2_F, permuted=TRUE, 
                                   pars=c("mufor","lambdahat","MHat"))

#Create replicate Data  
InSampleData <- TestDataFun(TotalData, sex="female", 
                            LastYearObs=2016, 
                            AdjMatType = 3)
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


####### -----------   MODEL EVALUATION ------------------- #####################  
#Get Out of Sample Data
DataOOS <- OutOfSampleData(Data = TotalData,sex="female", 
                           LastYearObs = 2016, h=1)


FCMatStanAPC_BYM2 <- StanFit_APC_BYM2$mufor[1:1000,]

#PIT of Model
FY <- EmpCDFFun(ObservedCount = DataOOS$D, 
                FCMat=FCMatStanAPC_BYM2, Exposure = DataOOS$ExposureFC)
FY1 <- EmpCDFFun(ObservedCount = (DataOOS$D-1), 
                 FCMat=FCMatStanAPC_BYM2, Exposure = DataOOS$ExposureFC)

#plot PIT Value (see Appendix)
par(cex.axis=1.2, cex.lab=1.5, cex.main=1.5, cex.sub=1)
plot(pit(J=20,x=DataOOS$D, Px=FY,Pxm1=FY1),
     ylab="Relative frequency", ylim=c(0,2.5),
     main="APC_BYM2 Model")
abline(h=1, lty=2)



### Model Evaluation ###
set.seed(420)
PQuant <- PIlevel(80) #80% Quantil for coverage

#Get Out of Sample Data
DataOOS <- OutOfSampleData(Data = TotalData, Region="Bayern",
                           sex="female", LastYearObs = 2016, h=1)

#Create list with Models to be evaluated
ModelList <- list(APC_BYM2 = FCMatStanAPC_BYM2)

#Create a data Frame in Long Format for Model Evaluation 
TotFCDatFrame <- ModelList %>% 
  lapply(., function(x) FCDataFrame(x, Exposure = DataOOS$ExposureFC,
                                    ObservedCount = DataOOS$D, PQuant=PQuant)) %>% 
  bind_rows(., .id="ModelID") %>%
  mutate(ModelName = transformer(ModelID,names(ModelList)))

# Add Information for Out of Sample Evaluation
TotFCDatFrame <- TotFCDatFrame %>% 
  mutate("Region"=rep(DataOOS$FCSubset$RegionName,length(unique(ModelName))), 
         "RegionNum"=rep(DataOOS$FCSubset$RegionNumber,length(unique(ModelName))), 
         "Age"=rep(DataOOS$FCSubse$AgeID,length(unique(ModelName))),#Age
         "RegID"=rep(DataOOS$FCSubset$KreisID,length(unique(ModelName))),#RegionID,
         "D"=rep(DataOOS$FCSubset$Deaths,length(unique(ModelName))), #Deaths
         "Jahr"=rep(DataOOS$FCSubset$Year, length(unique(ModelName)))) #Year

#Model Evaluation
TotFCDatFrame %>% group_by(ModelName) %>% #group_by 
  summarise("MeanLog"=mean(logScore),
            "MeanDSS"=mean(DSS),
            "MeanRPS"=mean(CRPSEmp),
            "MAE"=abs(D-MeanVal) %>% mean(),
            "RMSE"=sqrt(mean((D-MeanVal)^2)), #Root Mean Squared Error
            "Coverage"=pi_accuracy(PIL=PIL_mid,
                                   PIU=PIU_mid, yobs=D)$Success)





### ----- STACKING EXAMPLE -----------------------------------------------------

#Example for Males
#Calculation of Stacking Weights InSample 2001-2010, FC= 2011-2014
DataAPC_BYM2_Stan_1114<- StanData(Data=TotalData, 
                              LastYearObs = 2011, #last Year of InSample Data 
                              AdjMatType = 3,
                              sex="male", #Sex (female = female, männlich = male)
                              ModelType = "APC", #APC or RH
                              RegionType = "BYM2", #BYM2 or SAR
                              Cohort=TRUE, #Cohort Index True or False 
                              TFor = 3) #Forecast Horzion (how many years to forecast)

DataRH_BYM2_Stan_1114 <- StanData(Data=TotalData, 
                                  LastYearObs = 2011, #last Year of InSample Data 
                                  AdjMatType = 3,
                                  sex="male", #Sex (weiblich = female, männlich = male)
                                  ModelType = "RH", #APC or RH
                                  RegionType = "BYM2", #BYM2 or SAR
                                  Cohort=TRUE, #Cohort Index True or False 
                                  TFor = 3)

#See Appendix for Parameters of HMC Modelling For each Model
library(rstan)
options(mc.cores = parallel::detectCores())
StanAPC_BYM2_M_1114 <- rstan::stan(file.path(getwd(),"StanCode/APC_BYM2.stan"),
                              data=DataAPC_BYM2_Stan_1114,
                              chains = 4, iter=5000,warmup=2500, 
                              save_warmup=FALSE, thin=5,
                              control = list(adapt_delta = 0.81))


#See Appendix for Parameters of HMC Modelling For each Model
library(rstan)
options(mc.cores = parallel::detectCores())
StanRH_BYM2_M_1114 <- rstan::stan(file.path(getwd(),"StanCode/RH_BYM2.stan"),
                                   data=DataRH_BYM2_Stan_1114,
                                   chains = 4, iter=5000,warmup=2500, 
                                  save_warmup=FALSE, thin=5,
                                   control = list(adapt_delta = 0.81))


#Get Model Parameters and Calculation of Stacking Weights
FCMatStanAPC_BYM2_1114<- rstan::extract(StanAPC_BYM2_M_1114, permuted=TRUE, pars="mufor")$mufor
FCMatStanRH_BYM2_1114<- rstan::extract(StanRH_BYM2_M_1114, permuted=TRUE, pars="mufor")$mufor

#Evaluation of Stacking Matrix
ModelListStacking <- list(APC_BYM2=FCMatStanAPC_BYM2_1114,
                          RH_BYM2=FCMatStanRH_BYM2_1114)

DataOOS <- OutOfSampleData(Data = TotalData, Region="Bayern",
                           sex="male", LastYearObs = 2011, h=3)

SWeights1114 <- StackingWeights(ObservedCount = DataOOS$D,
                          FCMatList = ModelListStacking,
                          Exposure = DataOOS$ExposureFC)

#For Calculation of Matrix use StackingMat function with the according weights
#here In-Sample Matrix
InSampleStackingMat <- list(rstan::extract(FCMatStanAPC_BYM2_1114, permuted=TRUE, pars="MHat")$Mhat,
                            rstan::extract(FCMatStanRH_BYM2_1114, permuted=TRUE, pars="MHat")$Mhat) %>% 
                        StackingMat(Weights = SWeights1114, FCMatList = .)
