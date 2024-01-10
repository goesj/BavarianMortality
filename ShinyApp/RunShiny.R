if (!require("pacman")) install.packages("pacman")
pacman::p_load(shiny)
runApp(appDir = file.path(getwd(),"ShinyApp/app.R"))



