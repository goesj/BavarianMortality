if (!require("pacman")) install.packages("pacman")
pacman::pload(shiny)
runApp(appDir = file.path(getwd(),"ShinyApp/app.R"))



