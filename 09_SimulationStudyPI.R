if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,Qtools,future.apply)

### Load Functions
source("01_Functions.R")


##### Simulation Study #########################################################
#Function that runs the simulation study, to be repeated S times
SimFunction <- function(M){
  
  #True, "Observed" Data
  PoissonMat <- sapply(1:M, function(y) rpois(n = 2000,lambda = y))
  
  
  #### 2) Calculation of PI's ################
  PIArray <- array(data=0, dim=c(M, 2, 3),
                   dimnames = list(1:M,
                                   c("PIL","PIU"),
                                   c("CoherentHomburg","CoherentStandard",
                                     "PIMidStandard")))
  
  Coverage <- 0.8
  PIQuants <- PIlevel(0.8)
  
  #Draws on which to calculate the PI's
  PoiDraws <- sapply(1:M, function(y) rpois(n = 2000,lambda = y))
  
  #PI Coherent Homburg
  PIArray[,,1] <- PIFunctionCoherent_Ext(PoiDraws, PI = Coverage)
  #PI Own
  
  #PI Standard (with Quantiles)
  PIArray[,,2] <- cbind(PIFunction01_Ext(PoiDraws, prob = PIQuants$Lo),
                        PIFunction01_Ext(PoiDraws, prob = PIQuants$Up))
  
  #PI Mid 
  PIArray[,,3] <- cbind(PIFunction_Mid_Ext(PoiDraws, prob = PIQuants$Lo),
                        PIFunction_Mid_Ext(PoiDraws, prob = PIQuants$Up))
  
  #Multiple Performance Measures 
  # 1.) average excedance
  # 2.) mean
  # 3.) sd
  ResultMatrix <- matrix(0, nrow = 3, ncol = 3)
  for(i in 1:3){ #over all models
    CovReal <- sapply(1:M,function(z)
      mean(PIArray[z,1,i]<= PoissonMat[,z] & 
             PoissonMat[,z]<= PIArray[z,2,i])
    )
    ResultMatrix[i,1] <- mean(CovReal[which(CovReal>=Coverage)]-Coverage) #Avg. Exceedance
    ResultMatrix[i,2] <- mean(CovReal)
    ResultMatrix[i,3] <- sd(CovReal) #sample standard deviation
  }
  return(ResultMatrix) 
}

########## Run Simulation Study in parallel ####################################

### Attention, may take a while ...
S <- 200
set.seed(42)
plan(multisession)
ResList <- future.apply::future_lapply(1:S, function(x) SimFunction(500), 
                                       future.seed = TRUE)



######## Plot Results ##########################################################
do.call("rbind",ResList) %>% 
  data.frame() %>% 
  setNames(c("Average Exceedance", "Mean","Standard Deviation")) %>% 
  mutate(Type=
           rep(c("Coherent PI","Classical PI",
                 "Mid PI"),S),
         .before=1) %>% 
  mutate(Iteration = rep(1:S, each = 3),.after=1) %>% 
  pivot_longer(., cols = 3:ncol(.), names_to = "Measure", values_to = "Value") %>% 
  ggplot(., aes(x = Type, y = Value, fill = Type))+
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  facet_wrap(~Measure, scales = "free")+
  scale_fill_brewer(palette="Dark2")+
  theme_bw()+
  theme(axis.text = element_text(size = 25), #change size of axis text
        axis.text.x = element_text(size = 25, angle = 20),
        axis.title = element_text(size=25 , face="bold"), #change size of axis title
        panel.grid = element_line(colour = "grey92"),
        plot.title = element_text(hjust=0.5, size=25, face="bold"),
        panel.border = element_blank(),
        strip.text.x=element_text(face="bold", size = 22),
        strip.background = element_blank(),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.key.size = unit(2,"line"),
        legend.position = "bottom")

