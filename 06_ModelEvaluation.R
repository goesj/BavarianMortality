#### Load Data ####
load("Data/TotalData.RData") #load Data
### Load Functions
source("01_Functions.R")


####### -----------   MODEL EVALUATION ------------------- #####################  
load(file="../Results/StanAPC_BYM2_F.RData")
#Get Out of Sample Data (Test Data)
DataOOS <- OutOfSampleData(Data = TotalData,sex="female", 
                           LastYearObs = 2014, h=3)

#Get Samples from Posterior
StanFit_APC_BYM2 <- rstan::extract(StanAPC_BYM2_F, inc_warmup=FALSE)

#Get 1000 Draws of predicted mortality rates
FCMatStanAPC_BYM2 <- StanFit_APC_BYM2$mufor[1:1000,]

##Create PIT plot (see Appendix)
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

#Create list with Models to be evaluated (may be extended with multiple models)
ModelList <- list(APC_BYM2 = FCMatStanAPC_BYM2)

#Estimate Scores and Measures with FCDataFrame function
TotFCDatFrame <- ModelList %>% 
  lapply(., function(x) FCDataFrame(x, Exposure = DataOOS$ExposureFC,
                                    ObservedCount = DataOOS$D, PQuant=PQuant)) %>% 
  bind_rows(., .id="ModelID") %>% #bind sublists together
  mutate(ModelName = transformer(ModelID,names(ModelList))) #create Name ID column

# Add Information for Out of Sample Evaluation
TotFCDatFrame <- TotFCDatFrame %>% 
  mutate("Region"=rep(DataOOS$FCSubset$RegionName,length(unique(ModelName))), 
         "RegionNum"=rep(DataOOS$FCSubset$RegionNumber,length(unique(ModelName))), 
         "Age"=rep(DataOOS$FCSubse$AgeID,length(unique(ModelName))),#Age
         "RegID"=rep(DataOOS$FCSubset$KreisID,length(unique(ModelName))),#RegionID,
         "D"=rep(DataOOS$FCSubset$Deaths,length(unique(ModelName))), #Deaths
         "Jahr"=rep(DataOOS$FCSubset$Year, length(unique(ModelName)))) #Year

#Model Evaluation
TotFCDatFrame %>% group_by(ModelName) %>% #group_by Model
  summarise("MeanLog"=mean(logScore), #calculate mean scores
            "MeanDSS"=mean(DSS),
            "MeanRPS"=mean(RPSEmp),
            "MAE"=mean(abs(D-MeanVal)),
            "RMSE"=sqrt(mean((D-MeanVal)^2)), #Root Mean Squared Error
            "Coverage"=pi_accuracy(PIL=PICoherent.1,
                                   PIU=PICoherent.2, yobs=D)$Success)




