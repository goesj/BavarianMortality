## Libraries needed for Functions###
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,scoringRules,Qtools,MortCast,
               future.apply,spdep,ggpubr,INLA,bayesplot)

################ 1.1 Functions for Generation of Data ##########################

# Transform Neighborhood Matrix into List for Stan
CARData4Stan <- function(NeighborhoodMatrix){ #Region in Stan Matrix
  N <- nrow(NeighborhoodMatrix) #Amount of Regions
  N_edges <- sum(NeighborhoodMatrix) #Amount of Edges
  node1 <- vector(mode="numeric", length=N_edges)
  node2 <- vector(mode="numeric", length=N_edges)
  iEdge <- 1 #Helpvector
  for (j in 1:nrow(NeighborhoodMatrix)) {
    NumEdge <- sum(NeighborhoodMatrix[j,]) #Sum Neighbors
    TypeEdge <- which(NeighborhoodMatrix[j,]!=0) #Check which Region is Neighbor 
    node1[iEdge:(iEdge+NumEdge-1)] <- rep(j,NumEdge) 
    node2[iEdge:(iEdge+NumEdge-1)] <- TypeEdge 
    iEdge <- iEdge+NumEdge #Update Help Vector
  }
  
  node1Names <- rownames(NeighborhoodMatrix)[node1] 
  node2Names <- rownames(NeighborhoodMatrix)[node2] 
  
  return (list("N"=N,"N_edges"=N_edges,"node1"=node1,"node2"=node2,
               "node1Names"=node1Names, node2Names=node2Names))
} 

#Currently not in use !!!!
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
#Code taken from Morris (2019) 
#https://mc-stan.org/users/documentation/case-studies/icar_stan.html
ScalingFacBYM2 <- function(Nodes, AdjMat){
  #Add a small jitter to the diagonal for numerical stability 
  Q <- Diagonal(Nodes$N, rowSums(AdjMat)) - AdjMat
  # Compute the diagonal elements of the covariance matrix subject to the 
  # constraint that the entries of the ICAR sum to zero.
  #See the inla.qinv function help for further details.
  Q_pert <-  Q + Diagonal(Nodes$N) * max(diag(Q)) * sqrt(.Machine$double.eps)
  
  # Compute the diagonal elements of the covariance matrix subject to the 
  # constraint that the entries of the ICAR sum to zero.
  #See the inla.qinv function help for further details.
  Q_inv <-  INLA::inla.qinv(Q_pert, constr=list(A = matrix(1,1,Nodes$N),e=0))
  
  #Compute the geometric mean of the variances, which are on the diagonal of Q.inv
  scaling_factor <-  exp(mean(log(diag(Q_inv))))
  
  return(scaling_factor)
}


#Calculation of Cohort Index
CohortIndex <- function(agegroup, time, maxAge, M){
  #reshuffle age groups (the lowest two age groups are joined together)
  agegroup <- ifelse(agegroup==1,agegroup,agegroup-1) 
  return(k=M*(maxAge-agegroup)+time)
}

#Generation of Training Data
TrainingDataFun <- function(Data, sex, LastYearObs){
  
  #Filter Year >2000, since 2000 does not have values for exposure
  Data <- Data %>% filter(Year > 2000)
  
  #Filter Year and Sex
  Data <- Data %>% filter(., Sex=={{sex}} & Year <=LastYearObs)
  
  #create IDS
  Data$KreisID <- match(Data$RegionNumber,unique(Data$RegionNumber))
  Data$YearID <- match(Data$Year , unique(Data$Year))
  Data$AgeID <- match(Data$AgeGroup, unique(Data$AgeGroup)) #only lowest age 
  Data$REID <- 1:nrow(Data) #Random effect ID
  Data$CohortID <- CohortIndex(agegroup = Data$AgeID, 
                               time = Data$YearID, maxAge = max(Data$AgeID)-1,
                               M=5)
  
  
  Data$KreisYearID <- 
    as.numeric(interaction(Data$KreisID,Data$YearID, drop=TRUE))
  
  
  #Insert Nodes for BYM2 Model and SAR Model
  nm.adj <- poly2nb(Bayern)
  #Create Adjacency Matrix
  AdjMatBayern <- as(nb2mat(nm.adj, style = "B"), "Matrix") 
  #Get Nodes for Stan
  Nodes <- CARData4Stan(NeighborhoodMatrix = AdjMatBayern) 
  #scaling Factor BYM2 Model
  scalingFac <- ScalingFacBYM2(Nodes = Nodes, AdjMat = AdjMatBayern) 

  
  return(list(Data=Data, #Data
              Nodes=Nodes, #BYM2 /ICAR
              scalingFac=scalingFac #BYM2
              )) 
}

#Function for Gerneration of Data as Input for Stan Models
StanData <- function(Data, LastYearObs=2016, sex, 
                     ModelType, RegionType="BYM2", Cohort=TRUE, TFor=1){
  #Get actual Data
  DataInput <- TrainingDataFun(Data=Data, 
                           AdjMatType=AdjMatType, 
                           LastYearObs=LastYearObs, sex=sex)
  
  #Basics
  DataOut <- list("T"=max(DataInput$Data$YearID),"A"=max(DataInput$Data$AgeID),
                "R"=max(DataInput$Data$KreisID), "y"=DataInput$Data$Deaths,
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
  #Check for Regional Type (BYM2)
  if(RegionType=="BYM2"){
    DataOut <- append(DataOut,
                      list("N_edges"=DataInput$Nodes$N_edges,
                           "node1"=DataInput$Nodes$node1, 
                           "node2"=DataInput$Nodes$node2,
                           "scaling_factor"=DataInput$scalingFac))
  } 
  return(DataOut)
}


### Creation of Out of Sample Data
OutOfSampleData <- function(Data, sex="female", LastYearObs=2016, TFor=1){
  
  LastYearTest <- LastYearObs+TFor #Last Year of Test Data
  
  #Filter OOS Data
  FCSubset <- Data %>% filter(Year>LastYearObs & 
                              Year<=LastYearTest & 
                              Sex=={{sex}})
  
  
  FCSubset$KreisID <- match(FCSubset$RegionNumber,unique(FCSubset$RegionNumber))
  FCSubset$YearID <- match(FCSubset$Year , unique(FCSubset$Year))
  FCSubset$AgeID <- match(FCSubset$AgeGroup, unique(FCSubset$AgeGroup))
  
  #Get Deaths and Exposure
  D <- FCSubset$Deaths
  ExposureFC <- FCSubset$Exposure
  
  return(list("FCSubset"=FCSubset,
              "D"=D,
              "ExposureFC"=ExposureFC))
}

####### 1.2 Functions for Model Evaluation #####################################
#Get Levels of Prediction Intervals
PIlevel <- function(Percent){
  Percent <- ifelse(Percent>1,Percent/100,Percent)
  UpperBound <- 0.5+(Percent/2)
  LowerBound <- 0.5-(Percent/2)
  return(list(Up=UpperBound,
              Lo=LowerBound))
}

#Calculation of LogScore
logScore <- function(ObservedCount,FCMat, Exposure){
  logS <- numeric(length(ObservedCount)) #create empty vector
  for (i in 1:length(ObservedCount)) {
    #calculate lambda of all all posterior draws for each observation 
    lambdaVec <- exp(FCMat[,i])*Exposure[i] 
    #get mean log score of that observation
    logS[i] <- -log(mean(dpois(ObservedCount[i], lambdaVec))) 
    #alternative -dpois(y,lambda,log=TRUE) 
  }
  return(logS)
}
#DSS Score
DssScore <- function(ObservedCount, FCMat, Exposure){
  #Get Mean Value (law of iterated Expectation)
  MeanVal <- apply(FCMat, 2, function(x) mean(exp(x)))*Exposure 
  #Get Variance of lambda
  VarLambda <- apply(FCMat, 2, function(x) var(exp(x))) 
  # Law of total Variance
  SdVal <- sqrt(MeanVal + (Exposure)^2*VarLambda) 
  #see Czado. et al (2009)
  DSS <- ((ObservedCount - MeanVal) / SdVal)^2 + 2*log(SdVal) 
  return(DSS)
}

#ranked probability score
CRPSEmp <- function(ObservedCount, FCMat, Exposure){
  #create empty matrix
  YHat <- matrix(0, nrow = length(ObservedCount), ncol=nrow(FCMat)) 
  #draw deaths for each mu_i
  for (r in 1:nrow(YHat)) { #sample over all forecasted values
    YHat[r,] <- rpois(n=nrow(FCMat), lambda = exp(FCMat[,r])*Exposure[r]) 
  }
  DSSEmp <- scoringRules::crps_sample(y = ObservedCount, dat= YHat, method="edf") 
  return(DSSEmp)
}

#Empirical CDF (for PIT histogram)
EmpCDFFun <- function(ObservedCount, FCMat, Exposure){
  N <- length(ObservedCount)
  ECDF <- numeric(N)
  for (i in 1:N) {
    ECDF[i] <- mean(ppois(ObservedCount[i],exp(FCMat[,i])*Exposure[i]))
  }
  return(ECDF)
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
  
  erg <- list(breaks=(0:J)/J,counts=f_j, density=f_j,
              mids=(0:(J-1))/J+diff((0:J)/J)/2,
              xname="Probability Integral Transform",equidist=TRUE)
  class(erg) <- "histogram"
  return(erg)
}


#Function for Calculation of Prediction Interval via Monte Carlo
PIFunctionO1 <- function(FCMat,ExposureFC,prob=0.025){
  Drawvec <-  numeric(ncol(FCMat))
  for (r in 1:ncol(FCMat)) {
    Drawvec[r] <- rpois(n=nrow(FCMat), lambda = exp(FCMat[,r])*ExposureFC[r]) %>% 
                  quantile(probs=prob)
  }
  return(Drawvec)
}

PIFunctionCoherent <- function(FCMat, ExposureFC, PI){
  
  #Create Matrix for final PI's
  PIMatFinal <- matrix(0, nrow = ncol(FCMat),
                       ncol = 2)
  
  for(r in 1:ncol(FCMat)){ #loop over all forecasts
    
    #coherent PI's as described by Homburg et. al (2021)
    #Get Poisson Draws
    PoiDraws <- rpois(n=nrow(FCMat),
                      lambda = exp(FCMat[,r])*ExposureFC[r])
    
    #1.) First Compute the larges Integer L in N s.t. calculate P(X<L)>=1-Cov
    UniqueValues <- 0:max(PoiDraws) #Get entire Sample Range 
    
    #calculate P(X<L)
    ProbVals <- sapply(UniqueValues,function(x) mean(PoiDraws < x)) #Calculate probability
    
    #Find integer L (here L is integer not position in vector)
    L <- UniqueValues[sum(ProbVals<=(1-PI))] #largest integer for which this holds
    
    #2.) Calculate L+1 PI's 
    PIMat <- matrix(data = 0,nrow = L+1,ncol = 4,
                    dimnames = list(0:L,
                                    c("Lower","Upper","Coverage","Width")))
    
    #for l = 0,..,L compute the smallest u s.t. P(l <= X <= u)>= Cov
    ProbValsUpper <- sapply(UniqueValues, function(z) mean(PoiDraws <= z))
    
    #3. for each l find a minumum u
    for(i in 1:(L+1)){ #L+1 PI's 
      
      l <- i-1 #starting with l = 0, i is helper index
      
      PIMat[i,1] <- l #set l as lower bound 
      
      uVec <- L:max(PoiDraws) #possible values of u
      #find upper bound of u in L,...,max(Y)
      #P(X in [Xu,Xl]) = P(X <= Xu)- P(X < Xl)
      upperBound <- which(
        sapply(uVec, #Find u in L,...,max(Y) 
               function(y) (ProbValsUpper[y+1]-ProbVals[i]))>=
          PI) %>% min()
      
      PIMat[i,2] <- uVec[upperBound] #upper Bound (Integer)
      
      #Coverage P(X in Xu,Xl) = P(X <= Xu)- P(X < Xl)
      #plus 1 since upper bound is integer and not position in vector
      PIMat[i,3] <- ProbValsUpper[PIMat[i,2]+1]-  ProbVals[i]
      
      PIMat[i,4] <- PIMat[i,2]-PIMat[i,1] #Width
    }
    
    #4) chose PI of smallest Width, then choose the one with greatest Coverage
    PIMatFinal[r,] <- PIMat %>% data.frame() %>% 
      arrange(Width, desc(Coverage)) %>% #smallest Width and greatest Coverage
      slice(1) %>% select(1:2) %>% as.numeric() #Get values in final Matrix
  }
  return(PIMatFinal)
}


#Coverage 
pi_accuracy <- function(PIL,PIU, yobs){
  if(length(yobs) != length(PIL)){
    stop("yobs needs to be the same length as the forecast period.")
  }
  n <- length(yobs) #Length of Observations
  In <- (yobs >= PIL & yobs <= PIU) #Is estimate in Interval?
  ZeroPIU <- sum(PIL == 0 & yobs == 0) #Amount of zeros estimated
  In <- ifelse(PIL == 0 & yobs == 0,TRUE,In) #Check if lower PI is zero
  Success <- mean(In)
  return(list(In = In, Success = Success, n = n,
              ZeroPIU=ZeroPIU))
}


##Create Data Frame with Model Evalution Score for each Forecast.
# For each Forecast (only on test Data) multiple score are calculated  
FCDataFrame <- function(FCMat,Exposure,PQuant,ObservedCount){ #Function for the forecast dataFrame
  DataFrame <- data.frame("MeanVal"=(apply(FCMat,2, #Mean Forecasts
                                           function(x) mean (exp(x)))*Exposure),
                          "logScore"=logScore(ObservedCount,FCMat,Exposure), #Log Score 
                          "DSS"=DssScore(ObservedCount, FCMat, Exposure), #DSS Score
                          "RPSEmp"=CRPSEmp(ObservedCount,FCMat, Exposure), #RPS Score
                          "PICoherent" = PIFunctionCoherent(FCMat,         #PI's
                                                            ExposureFC = Exposure,
                                                            PI = 0.8),
                          "Obs"=1:ncol(FCMat)) #Observation Index
  return(DataFrame)
}

#Helper function to get Model Name
transformer <- function(x, values) {
  if(length(unique(x)) != length(values)){
    stop("Must be of same length")
  }
  return(values[match(x, unique(x))])
}


############## 1.3. Functions for calcualtion of Stacking Weights##############
#Helper Function
#(See https://mc-stan.org/docs/2_29/stan-users-guide/log-sum-of-exponentials.html 17.4.3)
log_sum_exp <- function(u) {
  max_u <- max(u);
  a <- 0;
  for (n in 1:length(u)) {
    a <- a + exp(u[n] - max_u);
  }
  return(max_u + log(a));
}

#Function to get log Likelihood of test Data
LoglikeMat <- function(ObservedCount, FCMat, Exposure){
  loglikeMat <- matrix(0, nrow = nrow(FCMat), ncol=length(ObservedCount)) #create empty vector
  for (i in 1:length(ObservedCount)) {
    lambdaVec <- exp(FCMat[,i])*Exposure[i] #calculate lambda 
    loglikeMat[,i] <- log(dpois(ObservedCount[i], lambdaVec)) #get mean log score
  }
  return(loglikeMat)
}

#Calculation of Stacking Weights##
#Taken from Barigou StanMoMo Package (https://github.com/kabarigou/StanMoMo/)
StackingWeights <- function(ObservedCount, FCMatList, Exposure){
  LoglikeList <- lapply(FCMatList, function(x) {
    LoglikeMat(ObservedCount = ObservedCount,FCMat = x,Exposure = Exposure)}
    )
  #Calculation of Mean of Log Like Values 
  lpd <- lapply(LoglikeList,
                function(x){apply(x,2,log_sum_exp)-log(nrow(x))}
                ) 
  lpd_point <- simplify2array(lpd) #Mean lpd Point values mean(log(p(y|y,theta)))
  Stacking <- StackFun(lpd_point = lpd_point)
  return(Stacking)
}

#Function from loo package stacking_weights. 
#slight Change so that results show Names of Models
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
                   names = colnames(lpd_point),  #add model names
                   class = c("stacking_weights"))
  return(wts)
}

#Function for creation of Stacking Matrix
StackingMat <- function(Weights, FCMatList){
  SMat <- Weights[1]*FCMatList[[1]] #copy Matrix
  for (i in 2:length(Weights)) {
    Add <- Weights[i]*FCMatList[[i]]
    SMat <- SMat+Add
  }
  return(SMat)
}



############## 1.4. Functions for Evaluation of Life Expecancy##################
#Get Life Expectancy for each posterior predictive draw
LifeExpFunIt <- function(FCMat, LastYear=2017, Year=Jahr, 
                         sex="female", OutOfSample=TRUE, hMax=15){
  
  #Creation of Empty Life Expectancy Matrix
  LifeExpMat <- matrix(0, nrow = 1000, #iterations
                       ncol= 96, #Regions
                       dimnames = list("It"=1:1000,
                                       "Reg"=unique(InSampleData$Data$RegionName)))
  if(OutOfSample==TRUE){
    #Index Helper to find Rows in Matrix
    IndexHelper <- expand.grid("Age"=unique(InSampleData$Data$AgeID),
                               "Region"=unique(InSampleData$Data$RegionName), 
                               "Year"=(LastYear+1):(LastYear+hMax)) 
    
    FCMat <- exp(FCMat) #exp out of Sample FC to get actual rates
  } else { #InSample Fit
    #Index Helper to find Rows in Matrix 
    IndexHelper <- expand.grid("Age"=unique(InSampleData$Data$AgeID), 
                               "Region"=unique(InSampleData$Data$RegionName),
                               "Year"=unique(InSampleData$Data$Year))
  }

  #Calculate LE for each iteration for all regions
  for (i in 1:ncol(LifeExpMat)) { #for all regions
    #1. Find Position in Matrix
    SubSetPosition <- which(IndexHelper$Year==Year & #Jahr 2017
                            IndexHelper$Region==colnames(LifeExpMat)[i])
  
    #Calculate Life Expecancy for all samples
    LifeExpMat[,i] <- 
      apply(FCMat[,SubSetPosition],1, 
            FUN = function(x) 
            MortCast::life.table(x,sex, abridged = TRUE, open.age = 100)[1,"ex"]) 
    
  }
  return(LifeExpMat)
  
}

## Function for calculation of Life Expectancy Quantiles
LifeExpSampleFunction2 <- function(FCMatIn, FCMatOut, PI, sex="female",
                                   YearsIn, YearsOut){
  #Get List with InSample LE
  InSampleList <- future.apply::future_lapply(YearsIn, function(x) 
    LifeExpFunIt(FCMat = FCMatIn,LastYear=2017, Year=x, 
                 sex=sex, OutOfSample = FALSE,hMax=15)) 
  

  #Get List with OutofSample LE
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
  QInList1 <- #Calcualte PI for smaller Inverval (50%-PI)
    lapply(InSampleList, function(x) 
      apply(x, 2, function(y) #calculate both mean and quantiles by Column
        c(mean(y), quantile(y, probs=c(PQuant1$Lo,PQuant1$Up))))) %>% 
    lapply(., function(z) as.data.frame(t(z))) %>% #transpose result
    lapply(., function(u) { #change columns names
      dimnames(u)[[2]] <- c("Mean","PiLo", "PiUp")
      u <- mutate(u, "Width"=PI[1]/100, #add width
                  "RegNumber"=unique(InSampleData$Data$RegionNumber))
      return(u)
    })
  
  QInList2 <- #Calcualte PI for wider Inverval (80%-PI)
    lapply(InSampleList, function(x) 
      apply(x, 2, function(y) #calculate both mean and quantiles by Column
        c(mean(y), quantile(y, probs=c(PQuant2$Lo,PQuant2$Up))))) %>% 
    lapply(., function(z) as.data.frame(t(z))) %>%  #transpose result
    lapply(., function(u) { #change columns names
      dimnames(u)[[2]] <- c("Mean","PiLo", "PiUp")
      u <- mutate(u, "Width"=PI[2]/100, #add width
                  "RegNumber"=unique(InSampleData$Data$RegionNumber))
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
                  "RegNumber"=unique(InSampleData$Data$RegionNumber))
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
                  "RegNumber"=unique(InSampleData$Data$RegionNumber))
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

############## 1.5. Function for Graphics ######################################
#Function for Plotting Geary C
GearyCPlot <- function(TotalData, LastYearObs,sex){
  #Get InSample Data
  InSampleData <- TrainingDataFun(TotalData, sex=sex, 
                              LastYearObs=LastYearObs, 
                              AdjMatType = 1)
  
  ConstantRate <- InSampleData$Data %>% 
    group_by(Jahr.R, AgeID) %>% 
    reframe("CRate"=sum(Deaths)/sum(Exposure)) #ungrouped dataframe
  
  InSampleData$Data <- 
    InSampleData$Data %>% 
    arrange(Jahr.R,AgeID) %>% #Arrange by Year, Age to past constant rate
    mutate("CR"=rep(ConstantRate$CRate, each = 96)) %>%  #paste each element R times
    mutate("ExCounts"=CR*Exposure) %>% #Expected Death Counts
    mutate("ZCounts"=(Deaths - ExCounts) / sqrt(ExCounts)) %>% 
    arrange(Jahr.R, KreisID,AgeID) #original layout
  
  #Get Neighborhood Structure
  nm.adj <- spdep::poly2nb(Bayern) #Neighborhood List Bayern
  ListBY <- spdep::nb2listw(nm.adj, style = "B") #Adjacency matrix in List Form
  
  #Create Result data Frame
  ResDat <- data.frame("Val"=numeric(max(InSampleData$Data$AgeID)*
                                       max(InSampleData$Data$YearID)),
                       "Sd"=numeric(max(InSampleData$Data$AgeID)*
                                      max(InSampleData$Data$YearID)))
  
  i <- 1
  for(t in unique(InSampleData$Data$Year)){
    for ( a in unique(InSampleData$Data$AgeID)){
      
      Ind <- which(InSampleData$Data$Year == t & 
                     InSampleData$Data$AgeID == a)
      
      #run geary test and check for spatial association
      Res <- spdep::geary.test(InSampleData$Data$ZCounts[Ind], listw = ListBY) 
      ResDat$Val[i] <- Res$estimate[1]
      ResDat$Sd[i] <- Res$estimate[3]
      ResDat$pVal[i] <- Res$p.value
      i <-  i+1
    }
  }
  
  #Store Result Data
  PlotData <- expand.grid("A"=1:max(InSampleData$Data$AgeID),
                          "t"=unique(InSampleData$Data$Year)) %>% 
    mutate("GearyC_Mean"=ResDat$Val,
           "GearyC_Sd"=ResDat$Sd)
  
  #Colors for Plot
  Col <- if(sex=="male") {c("#b3cde0","#011f4b")} else {c("#bf7fbf","#400040")} #select colors
  Main <- ifelse(sex=="male","Males","Females")
  
  #Plot Result
  P <- ggplot(data=PlotData, aes(x = t, y=GearyC_Mean, group=t))+
    geom_boxplot(fill=Col[1], color=Col[2])+
    xlab("Year")+
    theme_bw()+
    ggtitle(Main)
  
  return(P)
}

#Function for Creation of Yrep Density 
YRepDensity <- function(Draws, Deaths, sex){
  if(sex=="female"){
    bayesplot::color_scheme_set("purple")
  } else {
    bayesplot::color_scheme_set("blue")
  }
  g1 <- bayesplot::ppc_dens_overlay(y = Deaths, yrep =Draws)+xlim(0,10)
  g2 <- bayesplot::ppc_dens_overlay(y = Deaths, yrep =Draws)+ xlim(10, 20)
  g3 <- bayesplot::ppc_dens_overlay(y = Deaths, yrep =Draws)+ xlim(20, 100)
  g4 <- bayesplot::ppc_dens_overlay(y = Deaths, yrep =Draws)+ xlim(100, 300)
  
  return(
    ggpubr::ggarrange(g1,g2,
                      g3,g4,
                      ncol=2,
                      nrow=2)
  )
}

######### 1.6. PI Functions for OOS ############################################
#Adjusted PI Functions, where "external" Draws (i.e. Deaths) are used as input
# to calculate Prediction intervals
PIFunction01_Ext <- function(PoiDraws,prob=0.025){
  Drawvec <-  numeric(ncol(PoiDraws)) 
  for (r in 1:ncol(PoiDraws)) {
    Drawvec[r] <- PoiDraws[,r] %>% 
      quantile(probs=prob)
  }
  return(Drawvec)
}

#Function for calculation of Mid-quantiles
PIFunction_Mid_Ext <- function(PoiDraws,prob=0.025){
  Drawvec <- numeric(ncol(PoiDraws))
  for (r in 1:ncol(PoiDraws)) {
    Drawvec[r] <- Qtools::midquantile(PoiDraws[,r],
                                      probs=prob)$y #Calculation of midquantile
  }
  return(Drawvec)
}

#Adjusted PI Function that runs faster. Takes less time for simulation study
PIFunctionCoherent_Ext <- function(PoiDraws, PI){
  
  #Create Matrix for final PI's
  PIMatFinal <- matrix(0, nrow = ncol(PoiDraws),
                       ncol = 2)
  
  for(r in 1:ncol(PoiDraws)){ #loop over all forecasts
    
    #coherent PI's as described by Homburg et. al (2021)
    #Get Poisson Draws
    
    #1.) First Compute the larges Integer L in N s.t. calculate P(X<L)>=1-Cov
    UniqueValues <- 0:max(PoiDraws[,r]) #Get entire Sample Range 
    
    #calculate P(X<L) 
    ProbVals <- sapply(UniqueValues,function(x) mean(PoiDraws[,r] < x))
    
    #Find integer L (here L is integer not position in vector)
    L <- UniqueValues[sum(ProbVals<=(1-PI))] #largest integer for which this holds
    
    #2.) Calculate L+1 PI's 
    PIMat <- matrix(data = 0,nrow = L+1,ncol = 4,
                    dimnames = list(0:L,
                                    c("Lower","Upper","Coverage","Width")))
    
    #for l = 0,..,L compute the smallest u s.t. P(l <= X <= u)>= Cov
    ProbValsUpper <- sapply(UniqueValues, function(z) mean(PoiDraws[,r] <= z))
    
    #Check if multiple low values (0,..,) have probability mass of 0
    #saves compuation time since entire loop does not have to be evaluated
    if(!all(ProbValsUpper>0)){ 
      #If lambda is large, multiple values from 0 onward will have a mass of zero
      LastValZero <- max(which(ProbValsUpper<=0)) #Last Values of Zero in Vector
      
      PIMat[1:(LastValZero),1] <- 0:(LastValZero-1) #includes the zero, hence -1
      
      uVec <- L:max(PoiDraws[,r]) #possible values of u
      
      upperBound <- which(
        sapply(uVec, #Find u in L,...,max(Y) 
               function(y) (ProbValsUpper[y+1]-ProbVals[1]))>=
          PI) %>% min()
      
      PIMat[1:(LastValZero),2] <- uVec[upperBound]
      PIMat[1:(LastValZero),3] <- ProbValsUpper[uVec[upperBound] + 1] -
        ProbVals[1]
      
      for(i in (LastValZero+1):(L+1)){ 
        # over all PI's 
        
        l <- i-1 #helper Index for Bound
        
        PIMat[i,1] <- l #set l as lower bound 
        
        uVec <- L:max(PoiDraws[,r]) #possible values of u
        #find upper bound of u in L,...,max(Y)
        #P(X in [Xu,Xl]) = P(X <= Xu)- P(X < Xl)
        upperBound <- which(
          sapply(uVec, #Find u in L,...,max(Y) 
                 function(y) (ProbValsUpper[y+1]-ProbVals[i]))>=
            PI) %>% min()
        
        PIMat[i,2] <- uVec[upperBound] #upper Bound (Integer)
        
        #Coverage P(X in Xu,Xl) = P(X <= Xu)- P(X < Xl)
        #plus 1 since upper bound is integer and not position in vector
        PIMat[i,3] <- ProbValsUpper[PIMat[i,2]+1]-  ProbVals[i]
        
        #PIMat[i,4] <- PIMat[i,2]-PIMat[i,1] #Width
      }
    } else { #do normal loop
      #3. for each l find a minumum u
      for(i in 1:(L+1)){ 
        # over all PI's 
        
        l <- i-1 #starting with l = 0, i is helper index
        
        PIMat[i,1] <- l #set l as lower bound 
        
        uVec <- L:max(PoiDraws[,r]) #possible values of u
        #find upper bound of u in L,...,max(Y)
        #P(X in [Xu,Xl]) = P(X <= Xu)- P(X < Xl)
        upperBound <- which(
          sapply(uVec, #Find u in L,...,max(Y) 
                 function(y) (ProbValsUpper[y+1]-ProbVals[i]))>=
            PI) %>% min()
        
        PIMat[i,2] <- uVec[upperBound] #upper Bound (Integer)
        
        #Coverage P(X in Xu,Xl) = P(X <= Xu)- P(X < Xl)
        #plus 1 since upper bound is integer and not position in vector
        PIMat[i,3] <- ProbValsUpper[PIMat[i,2]+1]-  ProbVals[i]
        
        #PIMat[i,4] <- PIMat[i,2]-PIMat[i,1] #Width
      } 
    } 
    PIMat[,4] <- PIMat[,2]-PIMat[,1]
    
    #4) chose PI of smallest Width, then choose the one with greatest Coverage
    PIMatFinal[r,] <- PIMat %>% data.frame() %>% 
      arrange(Width, desc(Coverage)) %>% #smallest Width and greatest Coverage
      slice(1) %>% select(1:2) %>% as.numeric() #Get values in final Matrix
  }
  return(PIMatFinal)
}

