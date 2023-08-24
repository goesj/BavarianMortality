#### Load Data ####
load("Data/TotalData.RData") #load Data
### Load Functions
source("01_Functions.R")
pacman::p_load(rstan)

######---- Create Data and run models----#######################################
#Example of Stacking For males
#Calculation of Stacking Weights InSample 2001-2010, FC= 2011-2014
DataAPC_BYM2_Stan_1114 <- StanData(Data=TotalData, 
                                  LastYearObs = 2011, #last Year of InSample Data 
                                  AdjMatType = 3,
                                  sex="male", #Sex 
                                  ModelType = "APC", #APC or RH
                                  TFor = 3) #Forecast Horzion 

DataRH_BYM2_Stan_1114 <- StanData(Data=TotalData, 
                                  LastYearObs = 2011, #last Year of InSample Data 
                                  AdjMatType = 3,
                                  sex="male", #Sex 
                                  ModelType = "RH", #APC or RH
                                  TFor = 3)

#See Appendix for Parameters of HMC Modelling For each Model
######--- RUN APC MODEL---######################################################
options(mc.cores = parallel::detectCores())
StanAPC_BYM2_M_1114 <- rstan::stan(file.path(getwd(),"StanCode/APC_BYM2.stan"),
                                   data=DataAPC_BYM2_Stan_1114,
                                   chains = 4, iter=5000,warmup=2500, 
                                   save_warmup=FALSE, thin=5,
                                   control = list(adapt_delta = 0.81))
save(StanAPC_BYM2_M_1114, 
     file="../Results/StanAPC_BYM2_M_1114.RData")
#####---- RUN RH MODEL -----####################################################
#See Appendix for Parameters of HMC Modelling For each Model
options(mc.cores = parallel::detectCores())
StanRH_BYM2_M_1114 <- rstan::stan(file.path(getwd(),"StanCode/RH_BYM2.stan"),
                                  data=DataRH_BYM2_Stan_1114,
                                  chains = 4, iter=5000,warmup=2500, 
                                  save_warmup=FALSE, thin=5,
                                  control = list(adapt_delta = 0.81))

save(StanRH_BYM2_M_1114, 
     file="../Results/StanRH_BYM2_M_1114.RData")


######---- LOAD DATA and calculate Stacking Weights#############################
load(file="../Results/StanAPC_BYM2_M_1114.RData")
load(file="../Results/StanRH_BYM2_M_1114.RData")

#Get Forecastes Mortality Rates
FCMatStanAPC_BYM2_1114<- rstan::extract(StanAPC_BYM2_M_1114, 
                                        permuted=TRUE, pars="mufor")$mufor
FCMatStanRH_BYM2_1114<- rstan::extract(StanRH_BYM2_M_1114, 
                                       permuted=TRUE, pars="mufor")$mufor

### Calculate Stacking Weights ##
#Create Model List
ModelListStacking <- list(APC_BYM2=FCMatStanAPC_BYM2_1114,
                          RH_BYM2=FCMatStanRH_BYM2_1114)

#Get Out Of Sample Data (Test Data)
DataOOS <- OutOfSampleData(Data = TotalData, Region="Bayern",
                           sex="male", LastYearObs = 2011, h=3)

#Calculate Stacking Weights 
SWeights1114 <- StackingWeights(ObservedCount = DataOOS$D,
                                FCMatList = ModelListStacking,
                                Exposure = DataOOS$ExposureFC)

#Calculate Stacking Matrix
# First Create List of Stacking Forecasts to combine (here IN-SAMPLE!!!!)
InSampleStackingMat <- list(rstan::extract(FCMatStanAPC_BYM2_1114, 
                                           permuted=TRUE, pars="MHat")$Mhat,
                            rstan::extract(FCMatStanRH_BYM2_1114, 
                                           permuted=TRUE, pars="MHat")$Mhat) %>% 
  StackingMat(Weights = SWeights1114, FCMatList = .) #combine to one Matrix



# Calculation of OOS Stacking Matrix #
OOSStackingMat <- list(FCMatStanAPC_BYM2_1114,
                       FCMatStanRH_BYM2_1114) %>% 
  StackingMat(SWeights1114, FCMatList=.)

###### --- Calculate Life Expectancy Quantiles for Stacking Matrix #############
LifeExpFrameMaleStacking11114 <- 
  LifeExpSampleFunction2(FCMatIn = InSampleStackingMat, 
                         FCMatOut = OOSStackingMat, 
                         PI=c(80,50), sex="male", 
                         YearsIn = 2001:2011,
                         YearsOut = 2012:2014,
                         LastYearIn = 2011)





