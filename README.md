# BavarianMortality
Estimation of regional mortality rates using Bayesian hierarchical models.
Regional Mortality is forecasted using Bayesian implementations of the APC and RH model. They are extended with the BYM2 prior to account for spatial dependencies. Both models are then used within a stacking approach to generate more robust forecasts. 


## Models 
There are a collection of Bayesian hierarchical models found in `StanCode` that were used for the analysis. 
Deaths are modelled using a poisson likelihood where the death rate is scaled to person-years of exposure. 

The models are fitted using **Stan** 

## Data 
Within the `Data` folder the following can be found: Deaths and Population data by Age, Sex, Time and Region which was downloaded from <https://www.statistikdaten.bayern.de/genesis/online/logon.> the database of the Bavarian statistical institute. 
Also within the subfolder `Maps` are shapefiles which are needed for the creation of Maps

In addition, the results of the final forecasts may be found as an RData type called *LifeExpFrameStacking.RData*. These are stacked forecasts of multiple models for both males and females (see Text for details) and is needed for the visualization of Life Expecancy predictions and Shiny App.  

## R Code
There are multiple R files each requiring different packages to run.  

* `01_Functions.R` Includes all self written functions needed for both visualization and evaluation of the results.
* `02_LoadingData.R` Script whichs loads the data from excel, calculates person-years of exposure and transforms into long data format. Finished Data for Analysis is then saved as *TotalData.RData* which can also be found in the `Data` folder.
* `03_Geary'sCPlot.R` Code for the creation of Geary's C plot as shown in the Appendix. 
* `04_RunModel.R` Exemplary Code to run and save one of the Models found in the `StanCode` subfolder.
* `05_PPC.R` Code to run all posterior predictive checks as described in the text. 
* `06_ModelEvaluation` Code for the forecast evaluation of one or multiple models as described in the text. 
* `07_Stacking` Exemplary code for the creation of stacked forecasts using leave-future-out validation.  
* `08_LifeExpectancyVisualization.R` Script for all Visulaizationsplots found in the PDF. 

The parameteric model code requires the following packages: 
`rstan`, `tidyverse`, `sf`, `INLA`

The analysis part needs the following packages (in addition to the above ones)
`spdep`, `Qtools`, `scoringRules`, `MortCast`, `future.apply`

The visualizations part needs 
`RColorBrewer`, `cowplot`, `ggdist`, `ggridges`, `ggpubr`,`ggrepel`

## Shiny App
A shiny app for better visualization of all results can be found in the `ShinyApp` folder. 
To run the shiny app, simply open the `RunShiny.R` file and run the entire code. Note, that the app needs access to the `Data`folder. 

