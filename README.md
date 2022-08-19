# BavarianMortality
Estimation of regional mortality rates using Bayesian hierarchical models

## Models 
There are a collection of bayesian hierarchical models found in `StanCode` that were used for the analysis. 
Deaths are modelled using a poisson likelihood where the death rate is scaled to person-years of exposure. 

The models are fitted using **Stan** 

## Data 
Within the `Data` folder the following can be found: Deaths and Population data by Age, Sex, Time and Region which was downloaded from <https://www.statistikdaten.bayern.de/genesis/online/logon.> the database of the Bavarian statistical institute. 
Also within the subfolder `Maps` are shapefiles which are needed for the creation of Maps

In addition, the results of the final forecasts may be found as an RData type called *LifeExpFrameStacking.RData*. These are stacked forecasts of multiple models for both males and females (see Text for details) and is needed for the visualization of Life Expecancy predictions and Shiny App.  

## R Code
There are multiple R files each requiring different packages to run.  

* `01_Functions.R` Includes all self written functions needed for both visualization and evaluation of the results.
* `02_LoadingData.R` Script whichs loads the data from excel, calculates person-years of exposure and transforms into long data format. Finished Data for Analysis is then saved as *TotalData.RData* which can also be found in the `Data` folder 
* `03_MortModelling.R` Examplary fitting and Evaluation of a Model for Females
* `04_LifeExpectancyVisualization.R` Script for all Visulaizationsplots found in the PDF. 

The parameteric Model Code requires the following packages: 
`rstan`, `tidyverse`, `sf`, `INLA`

The analysis part needs the following packages (in addition to the above ones)
`spdep`, `Qtools`, `scoringRules`, `MortCast`, `future.apply`

The visualizations part needs 
`RColorBrewer`, `cowplot`, `ggdist`, `ggridges`; `ggpubr`
