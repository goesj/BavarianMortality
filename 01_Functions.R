## Functions needed in Main Script ##

## For Data Generation ##
# From Matrix into Neighborhood Structure for Stan
CARData4Stan <- function(NeighborhoodMatrix){ #Region in Stan Matrix
  N <- nrow(NeighborhoodMatrix) #Amount of Regions
  N_edges <- sum(NeighborhoodMatrix) #Amount of Edges
  node1 <- vector(mode="numeric", length=N_edges)
  node2 <- vector(mode="numeric", length=N_edges)
  iEdge <- 1 #Helpvector
  for (j in 1:nrow(NeighborhoodMatrix)) {
    NumEdge <- sum(NeighborhoodMatrix[j,]) #Sum Neighbors
    TypeEdge <- which(NeighborhoodMatrix[j,]!=0) #Overview which Neighbor 
    node1[iEdge:(iEdge+NumEdge-1)] <- rep(j,NumEdge) 
    node2[iEdge:(iEdge+NumEdge-1)] <- TypeEdge 
    iEdge <- iEdge+NumEdge #Update Help Vector
  }
  
  node1Names <- rownames(NeighborhoodMatrix)[node1] 
  node2Names <- rownames(NeighborhoodMatrix)[node2] 
  
  return (list("N"=N,"N_edges"=N_edges,"node1"=node1,"node2"=node2,
               "node1Names"=node1Names, node2Names=node2Names))
} 

#Function for Creating of Adjacency Matrix (for SAR MODEL)
AdjMatFun <- function(AdjMat, type=2, TypeVek){
  AdjMatNew <- as.matrix(AdjMat) #Create new Matrix
  WMat <- matrix(0, nrow = nrow(AdjMatNew), ncol=ncol(AdjMatNew))
  if(type==1){ #type 1: Maximum counts
    for (i in 1:nrow(AdjMat)) {
      for (j in 1:ncol(AdjMat)) {
        if(i==j){ #0 on diagonal 
          next
        }else{
          AdjMatNew[i,j] <- AdjMat[i,j]* #Check if neighbor
            max(TypeVek[i],TypeVek[j]) #multiply with highest values 
        }
      }
    }
  } else if(type==2) { #type 2: Distance counts 
    for (i in 1:nrow(AdjMat)) {
      for (j in 1:ncol(AdjMat)) {
        if(i==j){  #0 on diagonal 
          next
        }else{
          AdjMatNew[i,j] <- AdjMat[i,j]* #Check if neighbor
            (abs(TypeVek[i]-TypeVek[j])+1) #multiply with highest values 
        }
      }
    }
  } else { #type 3: Normal Neighborhood structure
    AdjMatNew <- AdjMat
  }
  #Row Standardization
  for(i in 1:nrow(AdjMatNew)){
    WMat[i,] <- AdjMatNew[i,]/sum(AdjMatNew[i,])
  }
  return(list("WeightMat"=AdjMatNew,
              "RowStandMat"=WMat))
}

#Compute Scaling Factor for BYM2 Model
ScalingFacBYM2 <- function(Nodes, AdjMat){
  #Add a small jitter to the diagonal for numerical stability (optional but recommended)
  Q <- Diagonal(Nodes$N, rowSums(AdjMat)) - AdjMat
  # Compute the diagonal elements of the covariance matrix subject to the 
  # constraint that the entries of the ICAR sum to zero.
  #See the inla.qinv function help for further details.
  Q_pert <-  Q + Diagonal(Nodes$N) * max(diag(Q)) * sqrt(.Machine$double.eps)
  
  # Compute the diagonal elements of the covariance matrix subject to the 
  # constraint that the entries of the ICAR sum to zero.
  #See the inla.qinv function help for further details.
  Q_inv <-  inla.qinv(Q_pert, constr=list(A = matrix(1,1,Nodes$N),e=0))
  
  #Compute the geometric mean of the variances, which are on the diagonal of Q.inv
  scaling_factor <-  exp(mean(log(diag(Q_inv))))
  
  return(scaling_factor)
}



#Calculation of cohort Index
CohortIndex <- function(agegroup, time, maxAge, M){
  #reshuffle age groups (the lowest two age groups are joined together)
  agegroup <- ifelse(agegroup==1,agegroup,agegroup-1) 
  return(k=M*(maxAge-agegroup)+time)
}

## Function for Creation of Test Data
TestDataFun <- function(Data,Sex,LastYearObs,Region="Bayern", AdjMatType=2){
  if(Region!="Bayern"){ #OberFranken
    Data <- Data %>% filter(Kreis.Nummer %in% OberFranken$ARS & 
                              Jahr.R > 2000)
  } else { #Bayern
    Data <- Data %>% filter(Jahr.R > 2000)
  }
  #Filter Year and Sex
  Data <- Data %>% filter(., Geschlecht==Sex & Jahr.R <=LastYearObs)
  
  #create IDS
  Data$KreisID <- match(Data$Kreis.Nummer,unique(Data$Kreis.Nummer))
  Data$YearID <- match(Data$Jahr.R , unique(Data$Jahr.R))
  Data$AgeID <- match(Data$AgeGroup, unique(Data$AgeGroup)) #only lowest age group
  Data$REID <- 1:nrow(Data) #Random effect ID
  Data$CohortID <- CohortIndex(agegroup = Data$AgeID, 
                               time = Data$YearID, maxAge = max(Data$AgeID)-1,
                               M=5)
  
  
  Data$KreisYearID <- as.numeric(interaction(Data$KreisID,Data$YearID, drop=TRUE))
  
  
  
  #Insert Nodes for BYM2 Model and SAR Model
  if(Region!="Bayern"){ #Oberfranken
    nm.adjOF <- poly2nb(OberFranken)
    AdjMatOberFranken <- as(nb2mat(nm.adjOF, style = "B"), "Matrix")
    Nodes <-  CARData4Stan(NeighborhoodMatrix = AdjMatOberFranken) 
    scalingFac <- ScalingFacBYM2(Nodes = Nodes, AdjMat = AdjMatOberFranken)
    
    AdjMatSAR <- AdjMatFun(AdjMat = AdjMatOberFranken, type=AdjMatType, TypeVek = OberFranken$SN_KTYP4)
    
  } else {
    nm.adj <- poly2nb(Bayern)
    AdjMatBayern <- as(nb2mat(nm.adj, style = "B"), "Matrix")
    Nodes <- CARData4Stan(NeighborhoodMatrix = AdjMatBayern) 
    scalingFac <- ScalingFacBYM2(Nodes = Nodes, AdjMat = AdjMatBayern)
    
    AdjMatSAR <- AdjMatFun(AdjMat=AdjMatBayern, type=AdjMatType, TypeVek= Bayern$SN_KTYP4)
    
  }
  
  return(list(Data=Data, #Data
              Nodes=Nodes, #BYM2 /ICAR
              scalingFac=scalingFac, #BYM2
              AdjMat=AdjMatSAR)) #SAR
}

#Function for Gerneration of Data as Input for Stan Models
StanData <- function(Data, LastYearObs=2016,Sex, Region="Bayern",AdjMatType=3, 
                     ModelType, RegionType="BYM2", Cohort=TRUE, TFor=1){
  
  DataInput <- TestDataFun(Data=Data, Region=Region, AdjMatType=AdjMatType, LastYearObs=LastYearObs, Sex=Sex)
  
  #Basics
  DataOut <- list("T"=max(DataInput$Data$YearID),"A"=max(DataInput$Data$AgeID),
                "R"=max(DataInput$Data$KreisID),"y"=DataInput$Data$Deaths,
                "E"=DataInput$Data$Exposure,
                "TFor"=TFor)
  
  #Check for Model Type (APC or RH)
  if(ModelType=="APC"){
    DataOut <- append(DataOut,
                      list("TInd"=DataInput$Data$YearID, 
                           "AInd"=DataInput$Data$AgeID,
                           "RInd"=DataInput$Data$KreisID,
                           "ReInd"=1:nrow(DataInput$Data)))
    if(Cohort==TRUE){ #Cohort Index APC?
      DataOut <- append(DataOut,
                        list("CInd"=DataInput$Data$CohortID,"M"=5))
    } #else keep as it is
  } else if(Cohort == TRUE){ #Check for Cohort Index RH Family
    DataOut <- append(DataOut,
                      list("M"=5))
  }
  #Check for Regional Type (BYM2, SAR)
  if(RegionType=="BYM2"){
    DataOut <- append(DataOut,
                      list("N_edges"=DataInput$Nodes$N_edges,
                           "node1"=DataInput$Nodes$node1, 
                           "node2"=DataInput$Nodes$node2,
                           "scaling_factor"=DataInput$scalingFac))
  } else if(RegionType=="SAR"){
    DataOut <- append(DataOut,
                      list("W"=DataInput$AdjMat$RowStandMat,
                           "eigenWsar"=eigen(DataInput$AdjMat$RowStandMat)$values))
  }
  return(DataOut)
}


### Functions for Evaluation of Forecasts ###
OutOfSampleData <- function(Data, Region="Bayern", Sex="weiblich", LastYearObs=2016, h=1){
  LastYearTest <- LastYearObs+h
  if(Region!="Bayern"){
    FCSubset <- Data %>% filter(Data$Jahr.R>LastYearObs & Data$Jahr.R<=LastYearTest & Data$Geschlecht==Sex & 
                                  Data$Kreis.Nummer %in% OberFranken$ARS)
  } else {
    FCSubset <- Data %>% filter(Data$Jahr.R>LastYearObs & Data$Jahr.R<=LastYearTest & Data$Geschlecht==Sex)
  }
  
  FCSubset$KreisID <- match(FCSubset$Kreis.Nummer,unique(FCSubset$Kreis.Nummer))
  FCSubset$YearID <- match(FCSubset$Jahr.R , unique(FCSubset$Jahr.R))
  FCSubset$AgeID <- match(FCSubset$AgeGroup, unique(FCSubset$AgeGroup))
  
  D <- FCSubset %>% dplyr::select(Deaths) %>% pull()
  ExposureFC <- FCSubset %>%  dplyr::select(Exposure) %>% pull()
  return(list("FCSubset"=FCSubset,
              "D"=D,
              "ExposureFC"=ExposureFC))
}

#Get PI Level
PIlevel <- function(Percent){
  Percent <- ifelse(Percent>1,Percent/100,Percent)
  UpperBound <- 0.5+(Percent/2)
  LowerBound <- 0.5-(Percent/2)
  return(list(Up=UpperBound,
              Lo=LowerBound))
}

PointMeasures <- function(ObservedCount, FCMat, Exposure){
  YHat <- matrix(0, nrow = length(ObservedCount), ncol=nrow(FCMat)) #create empty matrix
  MAEVek <- numeric(length(ObservedCount))
  RMSEVek <- numeric(length(ObservedCount))
  for (r in 1:nrow(YHat)) { #sample over all forecasted values
    YHat[r,] <- rpois(n=nrow(FCMat), lambda = exp(FCMat[,r])*Exposure[r]) #draw deaths for each mu_i
    MAEVek[r] <- abs(YHat[r,]-ObservedCount[r]) %>% mean() #MAE über alle Iterationen
    RMSEVek[r] <- (YHat[r,]-ObservedCount[r])^2 %>% mean() %>% sqrt() #RMSE über alle Iterationen
    
  }
  return(list("MAe"=MAEVek,
              "RMSE"=RMSEVek))
}
#Calculation of LogScore
logScore <- function(ObservedCount,FCMat, Exposure){
  logS <- numeric(length(ObservedCount)) #create empty vector
  for (i in 1:length(ObservedCount)) {
    lambdaVec <- exp(FCMat[,i])*Exposure[i] #calculate lambda of first observation (all posterior draws)
    #alternative -dpois(y,lambda,log=TRUE) (scoring Rules)
    logS[i] <- -log(mean(dpois(ObservedCount[i], lambdaVec))) #get mean log score of that observation
  }
  return(logS)
}
#DSS Score
DssScore <- function(ObservedCount, FCMat, Exposure){
  MeanVal <- apply(FCMat, 2, function(x) mean(exp(x)))*Exposure #Get Mean Value (law of iterated Expectation)
  VarLambda <- apply(FCMat, 2, function(x) var(exp(x))) #Get Variance of lambda
  SdVal <- sqrt(MeanVal + (Exposure)^2*VarLambda) # Law of total Variance
  DSS <- ((ObservedCount - MeanVal) / SdVal)^2 + 2*log(SdVal) # see Czado. et al (2009)
  return(DSS)
}

DssEmp <- function(ObservedCount, FCMat, Exposure){
  YHat <- matrix(0, nrow = length(ObservedCount), ncol=nrow(FCMat))
  for (r in 1:nrow(YHat)) {
    YHat[r,] <- rpois(n=nrow(FCMat), lambda = exp(FCMat[,r])*Exposure[r])
  }
  DSSEmp <- dss_sample(y = ObservedCount, dat= YHat)
  return(DSSEmp)
}

#ranked probability score
CRPSEmp <- function(ObservedCount, FCMat, Exposure){
  YHat <- matrix(0, nrow = length(ObservedCount), ncol=nrow(FCMat)) #create empty matrix
  for (r in 1:nrow(YHat)) { #sample over all forecasted values
    YHat[r,] <- rpois(n=nrow(FCMat), lambda = exp(FCMat[,r])*Exposure[r]) #draw deathds for each mu_i
  }
  DSSEmp <- scoringRules::crps_sample(y = ObservedCount, dat= YHat, method="edf") #Use Function of scoring rules
  return(DSSEmp)
}

#Empirical CDF
EmpCDFFun <- function(ObservedCount, FCMat, Exposure){
  N <- length(ObservedCount)
  ECDF <- numeric(N)
  for (i in 1:N) {
    ECDF[i] <- mean(ppois(ObservedCount[i],exp(FCMat[,i])*Exposure[i]))
  }
  return(ECDF)
}

#ranked probability score (Czado et. al (2009)_p.1257)
rps <- function(ObservedCount, FCMat, Exposure, kmax = 1500){
  RPS <- numeric(length(ObservedCount))
  for (i in 1:length(ObservedCount)) {
    kmax <- ObservedCount[i]+100
    #Calculate ECDF for all values of k
    ECDFVek <- apply(matrix(0:kmax, ncol = 1),1,  
                     function(x){
                       mean(ppois(x,exp(FCMat[,i])*Exposure[i])) 
                     })
    
    #Calulate value of ranked probability score
    RPS[i] <-  (ECDFVek - ifelse(ObservedCount[i] <= 0:kmax,1,0))^2 %>% 
                      sum()
    
  }
  
  return(RPS)
}

#Pit Histogram##
## Taken from Riebler INLA Group
# https://groups.google.com/g/r-inla-discussion-group/c/eSwZJ5Iegcc
## non-randomized version of the PIT histogram
##
## Params: 
## u - real number 0<=u<=1 
## Px - P(X <= x)
## Pxm1 - P(X <= x-1)
pit.one <- function(u,x,Px,Pxm1){
  F_u <- ifelse(u <= Pxm1 , 0, pmin(1,(u-Pxm1)/(Px-Pxm1) ) )
  F_u[x==0] <- pmin(1,u/Px)[x==0] 
  if(u == 1){
    F_u <- 1
  }
  if(u == 0){
    F_u <- 0
  }
  #print(F_u)
  return(mean(F_u))
}
## Params: 
## J - number of bins
pit <- function(J=10, x, Px, Pxm1){
  F_u.bar <- sapply((0:J)/J,pit.one, x=x, Px=Px, Pxm1=Pxm1)
  f_j <- J*diff(F_u.bar)
  
  erg <- list(breaks=(0:J)/J,counts=f_j, density=f_j,mids=(0:(J-1))/J+diff((0:J)/J)/2,
              xname="Probability Integral Transform",equidist=TRUE)
  class(erg) <- "histogram"
  return(erg)
}


#Function for Calculation of Prediction Interval via Monte Carlo
PIFunctionO1 <- function(FCMat,ExposureFC,prob=0.025){
  Drawvec <-  numeric(ncol(FCMat))
  for (r in 1:ncol(FCMat)) {
    Drawvec[r] <- rpois(n=nrow(FCMat), lambda = exp(FCMat[,r])*ExposureFC[r]) %>% quantile(probs=prob)
  }
  return(Drawvec)
}

# Function for Calculation of Mid Prediction Interval
PIFunctionO1_mid <- function(FCMat,ExposureFC,prob=0.025){
  Drawvec <- numeric(ncol(FCMat))
  for (r in 1:ncol(FCMat)) {
    Drawvec[r] <- Qtools::midquantile(rpois(n=nrow(FCMat), lambda = exp(FCMat[,r])*ExposureFC[r]),probs=prob)$y #Calculation of midquantile
  }
  return(Drawvec)
}

#Function for Creation of PI Interval (own Method)
PIFunction02 <- function(FCMat,ExposureFC,prob=0.025){
  QuantileVec <-  numeric(ncol(FCMat))
  for (r in 1:ncol(FCMat)) {
    Draws <- rpois(n=nrow(FCMat), lambda = exp(FCMat[,r])*ExposureFC[r]) #get draws
    ECDF <- Draws %>% ecdf(.) #compute the empirical cummulative distribution function
    QuantileVec[r] <-  Draws %>% quantile(.,prob=prob,type=1) #calculate the wanted quantile
    #check if quantile is zero
    if(QuantileVec[r]==0){
      next #do nothing since quantile cannot be lower than zero 
    } else if(ECDF(QuantileVec[r])>prob){ #check if quantile is bigger than probability ?
        Qobs <- ECDF(QuantileVec[r]) #calculate specific value of quantile
        
        #check if difference
        #knots(ECDF) # gives unique values of ECDF
        QValLower <- knots(ECDF)[which(knots(ECDF)==QuantileVec[r])- #find the value of the Quantile
                                   1] #go 1 position lower in knot vector (could mean going multiple values lower)
        Qobs1 <- ECDF(QValLower) 
        
        #calulate parameter of binomial draw
        ptilde <- (prob-Qobs1)/(Qobs-Qobs1)
        omega <- rbinom(n = 1,size=1,prob = ptilde)
        
        QuantileVec[r] <- omega*QuantileVec[r]+(1-omega)*QValLower
      }
  }
  return(QuantileVec)
}
#Calculation of Prediction Interval Method Anne (Insample)
PIFunction02_InSample <- function(Draws,prob=0.1){
  QuantileVec <-  numeric(ncol(Draws))
  for (r in 1:ncol(Draws)) {
    ECDF <- Draws[,r] %>% ecdf(.) #compute the empirical cumulative distribution function
    QuantileVec[r] <-  Draws[,r] %>% quantile(.,prob=prob,type=1) #calculate the wanted quantile
    #check if quantile is zero
    if(QuantileVec[r]==0){
      next #do nothing since quantile cannot be lower than zero 
    } else if(ECDF(QuantileVec[r])>prob){ #check if quantile is bigger than probability ?
      Qobs <- ECDF(QuantileVec[r]) #calculate specific value of quantile
      
      #check if difference
      #knots(ECDF) # gives unique values of ECDF
      QValLower <- knots(ECDF)[which(knots(ECDF)==QuantileVec[r])- #find the value of the Quantile
                                 1] #go 1 position lower in knot vector (could mean going multiple values lower)
      Qobs1 <- ECDF(QValLower) 
      
      #calulate parameter of binomial draw
      ptilde <- (prob-Qobs1)/(Qobs-Qobs1)
      omega <- rbinom(n = 1,size=1,prob = ptilde)
      
      QuantileVec[r] <- omega*QuantileVec[r]+(1-omega)*QValLower
    }
  }
  return(QuantileVec)
}

#Coverage 
pi_accuracy <- function(PIL,PIU, yobs){# checks the success of prediction intervals of an object of class 
  if(length(yobs) != length(PIL)){
    stop("yobs needs to be the same length as the forecast period.")
  }
  n <- length(yobs) #Length of Observations
  In <- (yobs >= PIL & yobs <= PIU) #Is estimate in Quantile
  ZeroPIU <- sum(PIL == 0 & yobs == 0)
  In <- ifelse(PIL == 0 & yobs == 0,TRUE,In) #Check if lower PI is zero
  Success <- mean(In)
  return(list(In = In, Success = Success, n = n,
              ZeroPIU=ZeroPIU))
}


##Create Data Frame of Forecasts##
FCDataFrame <- function(FCMat,Exposure,PQuant,ObservedCount){ #Function for the forecast dataFrame
  DataFrame <- data.frame("MeanVal"=(apply(FCMat,2, function(x) mean (exp(x)))*Exposure),
                          "logScore"=logScore(ObservedCount,FCMat,Exposure),
                          "DSS"=DssScore(ObservedCount, FCMat, Exposure),
                          "CRPSEmp"=CRPSEmp(ObservedCount,FCMat, Exposure),
                          "PIL"=PIFunction02(FCMat, ExposureFC = Exposure, prob=PQuant$Lo),
                          "PIL_mid"=PIFunctionO1_mid(FCMat, ExposureFC = Exposure, prob=PQuant$Lo),
                          "PIU"=PIFunction02(FCMat, ExposureFC = Exposure, prob=PQuant$Up),
                          "PIU_mid"=PIFunctionO1_mid(FCMat, ExposureFC = Exposure, prob=PQuant$Up),
                          "Obs"=1:ncol(FCMat))
  return(DataFrame)
}



### Life Expectancy Functions###
#Life Expectancy FC Function##
LifeExpFunPI <- function(FCMat,PILevel, LastYear=2017, Year=Jahr, sex="weiblich", OutOfSample=TRUE,hMax=15){
  LifeExpFrame <- data.frame("RegNumber"=unique(InSampleData$Data$Kreis.Nummer))
  rownames(LifeExpFrame) <- unique(InSampleData$Data$Kreise.Name)
  
  sexEnglish <- ifelse(sex=="weiblich","female","male") #translate sex 
  
  if(OutOfSample==TRUE){
    IndexHelper <- expand.grid("Age"=unique(InSampleData$Data$AgeID), #number of age effects
                               "Region"=unique(InSampleData$Data$Kreise.Name), #number region effect
                               "Year"=(LastYear+1):(LastYear+hMax)) #year effects
    
  } else { #InSample Fit
    IndexHelper <- expand.grid("Age"=unique(InSampleData$Data$AgeID), #number of age effects
                               "Region"=unique(InSampleData$Data$Kreise.Name), #number region effect
                               "Year"=unique(InSampleData$Data$Jahr.R))
    
  }
  
  for (i in 1:nrow(LifeExpFrame)) { #for all regions
    #1. Find Position in Matrix
    SubSetPosition <- which(IndexHelper$Year==Year & 
                              IndexHelper$Region==rownames(LifeExpFrame)[i])
    
    #Mean and Quantile Bestimmen Rate
    QuantileMat <- apply(FCMat[,SubSetPosition], 2 , function(x) c(mean(exp(x)),
                                                                   quantile(exp(x),
                                                                            prob=c(PILevel$Lo,PILevel$Up))))
    
    LifeExpFrame$Mean[i] <- #Mean life expectancy 
      MortCast::life.table(QuantileMat[1,], #mean Value
                           sex = sexEnglish, abridged = TRUE, open.age = 100)[1,"ex"] #Lifetable 
    
    LifeExpFrame$PIUp[i] <- #Upper Quantile life expectancy 
      MortCast::life.table(QuantileMat[2, 
      ],
      sex = sexEnglish, abridged = TRUE, open.age = 100)[1,"ex"] #Lifetable 
    
    LifeExpFrame$PILo[i] <- #Lower Quantile life expectancy 
      MortCast::life.table(QuantileMat[3, 
      ],
      sex = sexEnglish, abridged = TRUE, open.age = 100)[1,"ex"] #Lifetable 
    
  } 
  return(LifeExpFrame)
  
}

### with mcmc intervals ##
LifeExpFunIt <- function(FCMat, LastYear=2017, Year=Jahr, sex="weiblich", OutOfSample=TRUE,hMax=15){
  
  LifeExpMat <- matrix(0, nrow = 1000, #iterations
                       ncol=96, 
                       dimnames = list("It"=1:1000,"Reg"=unique(InSampleData$Data$Kreise.Name)))
  if(OutOfSample==TRUE){
    IndexHelper <- expand.grid("Age"=unique(InSampleData$Data$AgeID), #number of age effects
                               "Region"=unique(InSampleData$Data$Kreise.Name), #number region effect
                               "Year"=(LastYear+1):(LastYear+hMax)) #year effects
    FCMat <- exp(FCMat) #exponate out of sample FC
  } else {
    IndexHelper <- expand.grid("Age"=unique(InSampleData$Data$AgeID), #number of age effects
                               "Region"=unique(InSampleData$Data$Kreise.Name), #number region effect
                               "Year"=unique(InSampleData$Data$Jahr.R))
  }
  sexEnglish <- ifelse(sex=="weiblich","female","male") #translate sex 
  
  #Calculate LE for each iteration for all regions
  for (i in 1:ncol(LifeExpMat)) { #for all regions
    #1. Find Position in Matrix
    SubSetPosition <- which(IndexHelper$Year==Year & #Jahr 2017
                              IndexHelper$Region==colnames(LifeExpMat)[i])
  
    
    LifeExpMat[,i] <- 
      apply(FCMat[,SubSetPosition],1, FUN = function(x) 
                                   MortCast::life.table(x,sexEnglish, abridged = TRUE, open.age = 100)[1,"ex"])
    
    

  }
  return(LifeExpMat)
  
}

LifeExpFunIt2 <- function(FCMat, LastYear=2017, Year=Jahr, sex="weiblich", OutOfSample=TRUE,hMax=15){
  if(OutOfSample==TRUE){
    IndexHelper <- expand.grid("Age"=unique(InSampleData$Data$AgeID), #number of age effects
                               "Region"=unique(InSampleData$Data$Kreise.Name), #number region effect
                               "Year"=(LastYear+1):(LastYear+hMax)) #year effects
  FCMat <- exp(FCMat) #exponate Out Of Sample List
  } else {
    IndexHelper <- expand.grid("Age"=unique(InSampleData$Data$AgeID), #number of age effects
                               "Region"=unique(InSampleData$Data$Kreise.Name), #number region effect
                               "Year"=unique(InSampleData$Data$Jahr.R))
  }
  
  sexEnglish <- ifelse(sex=="weiblich","female","male") #translate sex
  RegNames <- unique(InSampleData$Data$Kreise.Name)
  LifeExpList <- 
    lapply(RegNames, function(u){ #get subsets via lapply
      FCMat[1:1000,which(IndexHelper$Year==Year & IndexHelper$Region==u)]}) %>% 
    future.apply::future_lapply(X=., function(y){ #apply function to each list using future apply
      apply(y,1, function(z){ #calculate Life Table for each iteration
        MortCast::life.table(z, sex = "male", abridged = TRUE, open.age = 100)[1,"ex"]
      })
    })
  names(LifeExpList) <- unique(InSampleData$Data$Kreise.Name)
  return(LifeExpList)

}

#Calculation of Stacking Weights##
#Taken from Barigou Stan MoMo
log_sum_exp <- function(u) {
  max_u <- max(u);
  a <- 0;
  for (n in 1:length(u)) {
    a <- a + exp(u[n] - max_u);
  }
  return(max_u + log(a));
}

#Function to get log Likelihood of test Values
LoglikeMat <- function(ObservedCount, FCMat, Exposure){
  loglikeMat <- matrix(0, nrow = nrow(FCMat), ncol=length(ObservedCount)) #create empty vector
  for (i in 1:length(ObservedCount)) {
    lambdaVec <- exp(FCMat[,i])*Exposure[i] #calculate lambda of first observation (all posterior draws)
    loglikeMat[,i] <- log(dpois(ObservedCount[i], lambdaVec)) #get mean log score of that observation
  }
  return(loglikeMat)
}
#Calulation of Stacking Weights
StackingWeights <- function(ObservedCount, FCMatList, Exposure){
  LoglikeList <- lapply(FCMatList, function(x) {LoglikeMat(ObservedCount = ObservedCount,
                                                              FCMat = x, 
                                                              Exposure = Exposure)})
  #Log sum exp more stable?
  #Calculation of Mean of Log Like Values (See https://mc-stan.org/docs/2_29/stan-users-guide/log-sum-of-exponentials.html 17.4.3)
  lpd <- lapply(LoglikeList,function(x){apply(x,2,log_sum_exp)-log(nrow(x))}) #Calculation of Mean
  lpd_point <- simplify2array(lpd) #Mean lpd Point values mean(log(p(y|y,theta)))
  Stacking <- StackFun(lpd_point = lpd_point)
  return(Stacking)
}

#Function from loo package. Changes so that results gives Names of Models
StackFun <- function (lpd_point, optim_method = "BFGS", optim_control = list()) {
  stopifnot(is.matrix(lpd_point))
  N <- nrow(lpd_point)
  K <- ncol(lpd_point)
  if (K < 2) {
    stop("At least two models are required for stacking weights.")
  }
  exp_lpd_point <- exp(lpd_point)
  negative_log_score_loo <- function(w) {
    stopifnot(length(w) == K - 1)
    w_full <- c(w, 1 - sum(w))
    sum <- 0
    for (i in 1:N) {
      sum <- sum + log(exp(lpd_point[i, ]) %*% w_full)
    }
    return(-as.numeric(sum))
  }
  gradient <- function(w) {
    stopifnot(length(w) == K - 1)
    w_full <- c(w, 1 - sum(w))
    grad <- rep(0, K - 1)
    for (k in 1:(K - 1)) {
      for (i in 1:N) {
        grad[k] <- grad[k] + 
          (exp_lpd_point[i, k] - exp_lpd_point[i, K])/(exp_lpd_point[i, ] %*% w_full)
      }
    }
    return(-grad)
  }
  ui <- rbind(rep(-1, K - 1), diag(K - 1))
  ci <- c(-1, rep(0, K - 1))
  w <- constrOptim(theta = rep(1/K, K - 1), f = negative_log_score_loo, 
                   grad = gradient, ui = ui, ci = ci, method = optim_method, 
                   control = optim_control)$par

  wts <- structure(c(w, 1 - sum(w)), 
                   names = colnames(lpd_point), 
                   class = c("stacking_weights"))
  return(wts)
}

StackingMat <- function(Weights, FCMatList){
  SMat <- Weights[1]*FCMatList[[1]] #copy Matrix
  for (i in 2:length(Weights)) {
    Add <- Weights[i]*FCMatList[[i]]
    SMat <- SMat+Add
  }
  return(SMat)
}


#Function for Plotting Geary C
#Function For the Creation of GeryC Plot
GearyCPlot <- function(TotalData, Sex, LastYearObs){
  #Get InSample Data
  InSampleData <- TestDataFun(TotalData, Sex=Sex, LastYearObs=LastYearObs,Region="Bayern", AdjMatType = 1)
  
  InSampleData$Data <- InSampleData$Data %>% 
    mutate("NormD" = (Deaths/Exposure)*1000)
  
  #Get Neighborhood Structure
  nm.adj <- spdep::poly2nb(Bayern) #Neighborhood List Bayern
  ListBY <- spdep::nb2listw(nm.adj, style = "W") #Adjacency matrix in List Form
  
  #Create Result data Frame
  ResDat <- data.frame("Val"=numeric(max(InSampleData$Data$AgeID)*max(InSampleData$Data$YearID)),
                       "Sd"=numeric(max(InSampleData$Data$AgeID)*max(InSampleData$Data$YearID)))
  
  i <- 1
  for(t in unique(InSampleData$Data$Jahr.R)){
    for ( a in unique(InSampleData$Data$AgeID)){
      
      Ind <- which(InSampleData$Data$Jahr.R == t & 
                     InSampleData$Data$AgeID == a)
      Res <- spdep::geary.test(InSampleData$Data$NormD[Ind], listw = ListBY) #do geary test and check for spatial association
      ResDat$Val[i] <- Res$estimate[1]
      ResDat$Sd[i] <- Res$estimate[3]
      ResDat$pVal[i] <- Res$p.value
      i <-  i+1
    }
  }
  
  #Store Result Data
  PlotData <- expand.grid("A"=1:max(InSampleData$Data$AgeID),
                          "t"=unique(InSampleData$Data$Jahr.R)) %>% 
    mutate("GearyC_Mean"=ResDat$Val,
           "GearyC_Sd"=ResDat$Sd)
  
  #Colors for Plot
  Col <- if(Sex=="männlich") c("#b3cde0","#011f4b") else c("#bf7fbf","#400040") #select colors
  Main <- ifelse(Sex=="männlich","Males","Females")
  
  #Plot Result
  P <- ggplot(data=PlotData, aes(x = t, y=GearyC_Mean, group=t))+
    geom_boxplot(fill=Col[1], color=Col[2])+
    xlab("Year")+
    theme_bw()+
    ggtitle(Main)
  
  return(P)
}



LifeExpSampleFunction2 <- function(FCMatIn, FCMatOut, PI, sex="weiblich",
                                   YearsIn, YearsOut){
  
  InSampleList <- future.apply::future_lapply(YearsIn, function(x) 
    LifeExpFunIt(FCMat = FCMatIn,LastYear=2017, Year=x, 
                 sex=sex, OutOfSample = FALSE,hMax=15)) 
  
  OutOfSampleList <- future.apply::future_lapply(YearsOut, function(x) 
    LifeExpFunIt(FCMat = FCMatOut,LastYear=2017, Year=x, 
                 sex=sex, OutOfSample = TRUE,hMax=15))
  
  #Bind Sublists by Columsn
  InSampleList <- lapply(InSampleList, bind_cols)
  OutOfSampleList <- lapply(OutOfSampleList, bind_cols)
  
  #Give List Names
  names(InSampleList) <- YearsIn
  
  
  #Append Last Insample Year to Out of Sample Year
  OutOfSampleList <- append(list(InSampleList[[length(InSampleList)]]),
                            OutOfSampleList)
  
  names(OutOfSampleList) <- c(tail(YearsIn,1),YearsOut)

  
  PQuant1 <- PIlevel(PI[1])
  PQuant2 <- PIlevel(PI[2])
  
  #Calculate Quantiles of Both Lists
  QInList1 <- 
    lapply(InSampleList, function(x) 
      apply(x, 2, function(y) #calculate both mean and quantiles by Column
        c(mean(y), quantile(y, probs=c(PQuant1$Lo,PQuant1$Up))))) %>% 
    lapply(., function(z) as.data.frame(t(z))) %>% #transpose result
    lapply(., function(u) { #change columns names
      dimnames(u)[[2]] <- c("Mean","PiLo", "PiUp")
      u <- mutate(u, "Width"=PI[1]/100, #add with
                  "RegNumber"=unique(InSampleData$Data$Kreis.Nummer))
      return(u)
    })
  
  QInList2 <- 
    lapply(InSampleList, function(x) 
      apply(x, 2, function(y) #calculate both mean and quantiles by Column
        c(mean(y), quantile(y, probs=c(PQuant2$Lo,PQuant2$Up))))) %>% 
    lapply(., function(z) as.data.frame(t(z))) %>%  #transpose result
    lapply(., function(u) { #change columns names
      dimnames(u)[[2]] <- c("Mean","PiLo", "PiUp")
      u <- mutate(u, "Width"=PI[2]/100, #add with
                  "RegNumber"=unique(InSampleData$Data$Kreis.Nummer))
      return(u)
    })
  
 #OUT OF SAMPLE
  QOutList1 <- 
    lapply(OutOfSampleList, function(x) 
      apply(x, 2, function(y) #calculate both mean and quantiles by Column
        c(mean(y), quantile(y, probs=c(PQuant1$Lo,PQuant1$Up))))) %>% 
    lapply(., function(z) as.data.frame(t(z))) %>% #transpose result
    lapply(., function(u) { #change columns names
      dimnames(u)[[2]] <- c("Mean","PiLo", "PiUp")
      u <- mutate(u, "Width"=PI[1]/100, #add with
                  "RegNumber"=unique(InSampleData$Data$Kreis.Nummer))
      return(u)
    })
  
  QOutList2 <- 
    lapply(OutOfSampleList, function(x) 
      apply(x, 2, function(y) #calculate both mean and quantiles by Column
        c(mean(y), quantile(y, probs=c(PQuant2$Lo,PQuant2$Up))))) %>% 
    lapply(., function(z) as.data.frame(t(z))) %>%  #transpose result
    lapply(., function(u) { #change columns names
      dimnames(u)[[2]] <- c("Mean","PiLo", "PiUp")
      u <- mutate(u, "Width"=PI[2]/100, #add with
                  "RegNumber"=unique(InSampleData$Data$Kreis.Nummer))
      return(u)
    })
  
  
  #add all to singe Data Frame 
  Out <- 
    rbind( #InSample
      rbind( #Bind InSample rows Together
        bind_rows(QInList1,.id="Year"),
        bind_rows(QInList2,.id="Year")
      ) %>% mutate("Type"="InSample"), #add Type
      #Out of Sample
      rbind( #Bind Out of Sample rows Together
        bind_rows(QOutList1,.id="Year"),
        bind_rows(QOutList2,.id="Year")
    ) %>% mutate("Type"="OutOfSample") #add Tyoe
    ) %>% mutate("Sex"=sex) #add Sex
  return(Out)
  
}

# LifeExpSampleFunction <- function(FCMatIn, FCMatOut, PI, sex="weiblich"){
#   PQuant1 <- PIlevel(PI[1])
#   InSampleList1 <- future.apply::future_lapply(2001:2017, function(x) 
#     LifeExpFunPI(FCMat = FCMatIn,PILevel = PQuant1,LastYear=2017, Year=x, 
#                  sex="weiblich", OutOfSample = FALSE,hMax=15)) 
#   
#   OutOfSampleList1 <- future.apply::future_lapply(2018:2032, function(x) 
#     LifeExpFunPI(FCMat = FCMatOut,PILevel = PQuant1,LastYear=2017, Year=x, 
#                  sex="weiblich", OutOfSample = TRUE,hMax=15)) 
#   
#   #hard coded for length of 2
#   PQuant2 <- PIlevel(PI[2])
#   InSampleList2 <- future.apply::future_lapply(2001:2017, function(x) 
#     LifeExpFunPI(FCMat = FCMatIn,PILevel = PQuant2,LastYear=2017, Year=x, 
#                  sex="weiblich", OutOfSample = FALSE,hMax=15)) 
#   
#   OutOfSampleList2 <- future.apply::future_lapply(2018:2032, function(x) 
#     LifeExpFunPI(FCMat = FCMatOut,PILevel = PQuant2,LastYear=2017, Year=x, 
#                  sex="weiblich", OutOfSample = TRUE,hMax=15)) 
#   
#   #Append Last Insample Year to Out of Sample Year
#   OutOfSampleList1 <- append(list(InSampleList1[[length(InSampleList1)]]),
#                              OutOfSampleList1)
#   
#   OutOfSampleList2 <- append(list(InSampleList2[[length(InSampleList2)]]),
#                              OutOfSampleList2)
#   browser()
#   
#   #add all to singe Data Frame 
#   Out <- 
#     rbind( #InSample
#       rbind(do.call("rbind", InSampleList1),
#             do.call("rbind", InSampleList2)) %>% 
#         mutate("width"=c(rep(PI[1]/100, 96*17), #width for geom_lineribbon
#                          rep(PI[2]/100, 96*17)),
#                "Type"="InSample"), #add Type
#       #Out of Sample
#       rbind(do.call("rbind", OutOfSampleList1),
#             do.call("rbind", OutOfSampleList2)) %>% 
#         mutate("width"=c(rep(PI[1]/100, 96*16), #width for geom_lineribbon
#                          rep(PI[2]/100, 96*16)), #16 years
#                "Type"= "OutOfSample") #add Tyoe
#     ) %>% 
#     mutate("Year"=c(rep(rep(2001:2017,each=96), length(PI)),
#                     rep(rep(2017:2032,each=96), length(PI))),
#            "Sex"=sex)
#   return(Out)
#   
# }