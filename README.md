# BavarianMortality
Estimation of regional mortality rates using Bayesian hierarchical models

## Models 
There are a collection of bayesian hierarchical models found in `StanCode` that were used for the analysis. 
Deaths are modelled using a poisson likelihood where the death rate is scaled using exposure. 

The models are fitted using **Stan** 

## Data 
Within the `Data` folder the following can be found: Deaths and Population data by Age, Sex, Time and Region which was downloaded from <https://www.statistikdaten.bayern.de/genesis/online/logon.> the database of the Bavarian statistical institute. 
Also within the subfolder `Maps` are shapefiles which are needed for the creation of Maps

In addition, the results of the final forecasts may be found as an RData type called *LifeExpFrameStacking.RData*. These are stacked forecasts of multiple models for both males and females (see Text for details) and is needed for the visualization of Life Expecancy predictions and Shiny App.  
