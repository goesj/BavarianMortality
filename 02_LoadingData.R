## Loading Data ##
library(sf);library(openxlsx); library(tidyverse)

### Regional Data####
GermanyData <- read_sf(file.path(getwd(),"Data/Maps/vg250_01-01.tm32.shape.ebenen/vg250_ebenen_0101/VG250_KRS.shp"),as_tibble = FALSE)

#Regionaltype 
RegionalType <- read.xlsx(file.path(getwd(),"Data/Maps/ge1000.tm32.shape/ge1000/ktyp4_1000/KTYP4_1000_Tabelle.xlsx"),
                        sheet = 1,colNames = TRUE)

BY <- which(GermanyData$SN_L=="09") #Only Bavarias Region
Bayern <- GermanyData[BY,] #Dataset with only Bavaria

#Add RegionalType to German Data
Bayern$SN_KTYP4 <- RegionalType[match(Bayern$ARS, RegionalType$RS),"SN_KTYP4"]
Bayern$KTYP4 <- RegionalType[match(Bayern$ARS, RegionalType$RS),"KTYP4"]

#Change direction( 4 highest 1 lowest)
Bayern$SN_KTYP4 <- ifelse(Bayern$SN_KTYP4==1,4,
                               ifelse(Bayern$SN_KTYP4==2,3,
                                      ifelse(Bayern$SN_KTYP4==3,2,1)))

## Population Data
PopulationM <- read.xlsx(xlsxFile = file.path(getwd(),"Data/BY_12411(Population)_2000-17.xlsx"),
                          sheet = "BevölkerungMännlich", colNames = TRUE)

PopulationW <- read.xlsx(xlsxFile = file.path(getwd(),"Data/BY_12411(Population)_2000-17.xlsx"),
                          sheet = "BevölkerungWeiblich", colNames = TRUE)


DeathM <- read.xlsx(xlsxFile = file.path(getwd(),"Data/BY_12613(Deaths)_2000-17.xlsx"),
                    sheet = "SterbefälleMännlich", colNames = TRUE)

DeathW <- read.xlsx(xlsxFile = file.path(getwd(),"Data/BY_12613(Deaths)_2000-17.xlsx"),
                    sheet = "SterbefälleWeiblich", colNames = TRUE)

AgeGroups <- length(colnames(PopulationM[-c(1:5)])) #AgeGroups
NamesAge <- colnames(PopulationM[-c(1:5)]) #Names
RowsData <- nrow(PopulationM) 

#Transform Data Into Long Format
PM <- PopulationM %>% pivot_longer(., cols = -(1:5),names_to = "AgeGroup", values_to="Population")
DM <- DeathM %>% pivot_longer(., cols = -(1:5),names_to = "AgeGroup", values_to="Death")

MaleData <- PM %>% mutate("Deaths"=pull(DM[,ncol(DM)])) #add last column to data.frame (Vector of Deaths)


#Same for Female Data
PF <- PopulationW %>% pivot_longer(., cols = -(1:5),names_to = "AgeGroup", values_to="Population")
DF <- DeathW %>% pivot_longer(., cols = -(1:5),names_to = "AgeGroup", values_to="Death")

FemaleData <- PF %>% mutate("Deaths"=pull(DF[,ncol(DF)]))

#Calculate Exposure
#Exposure= (N(t)-N(t-1)/log(N(t)/N(t-1)))
ExposureDfMale <- (MaleData$Population[which(MaleData$Jahr.R>2000)]- #N(T)
                     MaleData$Population[which(MaleData$Jahr.R<2017)])/ #N(t-1)
  log(MaleData$Population[which(MaleData$Jahr.R>2000)]/
        MaleData$Population[which(MaleData$Jahr.R<2017)])

ExposureDfFemale <- (FemaleData$Population[which(FemaleData$Jahr.R>2000)]- #N(T)
                       FemaleData$Population[which(FemaleData$Jahr.R<2017)])/ #N(t-1)
  log(FemaleData$Population[which(FemaleData$Jahr.R>2000)]/
        FemaleData$Population[which(FemaleData$Jahr.R<2017)])

#Attention if NaN, then Population N(t)=N(t-1), in that case, Exposure = N(t)
#Finished Female Dataframe
FemaleData <- FemaleData %>% mutate(Exposure=c(rep(NA,times=which(FemaleData$Jahr.R>2000)[1]-1), #NA in 2000 (times = first observation that is of year 2001)
                                               ifelse(is.nan(ExposureDfFemale),Population, ExposureDfFemale))) %>% #Input Exposure
  mutate(MxRate= Deaths/Exposure) #Berechnung Death Rate


#Finished Male Dataframe
MaleData <- MaleData %>% mutate(Exposure=c(rep(NA,times=which(MaleData$Jahr.R>2000)[1]-1), #NA in 2000 (times = (first observation that is of year 2001))
                                           ifelse(is.nan(ExposureDfMale),Population, ExposureDfMale))) %>% 
  mutate(MxRate= Deaths/Exposure)



#Combine the Datasets
TotalData <- rbind(FemaleData,MaleData)


### LIFE EXPECTANCY Bevölkerungsvorausberechnung Landesamt Statistik###
LE_RBVB <- read.xlsx(xlsxFile = file.path(getwd(),"Data/Lebenserwartung_regBVB_2020-2040.xlsx"),
                      startRow = 4)

#Create Helper Data for Extraction of Names
HelpData <- FemaleData %>% filter(Jahr.R==2001, AgeGroup == "unter.1")

#Add Proper Region Number and Region Name
LE_RBVB <- LE_RBVB %>% 
  mutate("Geschlecht"=ifelse(Geschlechtsgruppe==1, "männlich","weiblich"),
         "Kreis.Nummer" = HelpData %>% 
           slice(match(AGS.kurz,substr(HelpData$Kreis.Nummer,3,5))) %>%  #Match Part of Regional Number 
           select(Kreis.Nummer) %>%  #Select Column
           pull(),
         "Kreise.Name" = HelpData %>% 
           slice(match(AGS.kurz,substr(HelpData$Kreis.Nummer,3,5))) %>% 
           select(Kreise.Name) %>% 
           pull())


save(TotalData,
     Bayern,
     OberFranken,
     LE_RBVB, file = file.path(getwd(),"TotalData.RData"))
