if (!require("pacman")) install.packages("pacman")
pacman::p_load(shiny,RColorBrewer,cowplot,
               ggplot2, tidyverse, ggdist)

## First load the data and functions needed
source("../01_Functions.R")
#### Load Data ####
load(file = "../Data/TotalData.RData") #load Data
load(file = "../Data/LifeExpFrameStacking.RData") #Load Results

RegionList <- TotalData$RegionNumber %>% unique() %>% as.list()

NamesEnglish <- gsub("(Lkr)", "Region", unique(TotalData$RegionName)) #Change Lkr to Region
NamesEnglish <- gsub("(Krfr.St)", "City", NamesEnglish) #Change Krfr. St to City

names(RegionList) <- NamesEnglish  

Mypal <- RColorBrewer::brewer.pal(n = 11, name = "RdBu")#Create own Color Palette

# Define UI ----
ui <- fluidPage(
  titlePanel("Bavarian Mortality by Region"),
  
  
  fluidRow(

    column(4,
           selectInput(inputId = "Region",
                       label = strong("Region"),
                       choices = RegionList[sort(names(RegionList))], #sort by names
                       selected = "Aichach-Friedberg (Lkr)")),
    column(7, offset = 1,
           sliderInput(inputId = "Year",
                       label=strong("Year"),
                       min = 2001, max = 2030,
                       value=2017,  #start year 2017
                       step = 1, sep=""))
    ),
  fluidRow(
    column(8, offset = 3,
           h4("Life Expecancy at Birth", align="left"),
           tableOutput("LifeTable"))
  ),
  plotOutput("RegionPlot"),
  plotOutput("YearPlot")
)

# Define server logic required to draw plots
server <- function(input, output) {
  # 1. Output Map Plot
  output$RegionPlot <- renderPlot({ 
    Bayern <- Bayern %>% 
      mutate("LifeExpStack_F" = LifeExpFrameFemaleStacking %>% 
               filter(Year == input$Year & Width == 0.8) %>% 
               slice(1:96) %>% select(Mean) %>% pull(),
             "LifeExpStack_M" =  LifeExpFrameMaleStacking %>% 
               filter(Year == input$Year & Width == 0.8) %>% 
               slice(1:96) %>% select(Mean) %>% pull(),
             "Rank_F" = rank(-LifeExpStack_F, ties.method = "random"),
             "Rank_M" = rank(-LifeExpStack_F, ties.method = "random"),
             "SelectedRegion"=ifelse(RS == input$Region,"Yes","No"),
             "alphaRegion" = ifelse(RS==input$Region, 1,0.2))
    
    #Plot for Females
    g1F <- ggplot(data = Bayern) +
      geom_sf(aes(fill = LifeExpStack_F, color = SelectedRegion,alpha=alphaRegion), size = 0) +
      scale_color_manual(values = c("No"=NA,"Yes"="red"), guide="none")+
      scale_alpha_continuous(range=c(0.2,1), guide="none")+
      scale_fill_gradientn(name="LifeExpectany",colors=Mypal, #name is for legend
                           labels = function(x) {round(x,0)})+ #rounds values in legend
      ggtitle("Female")+
      theme_void() + 
      theme(plot.title = element_text(size = 15, face = "bold"),
            legend.position = "bottom") 
    
    #Plot for Males
    g1M <- ggplot(data = Bayern) +
      geom_sf(aes(fill = LifeExpStack_M, color = SelectedRegion,alpha=alphaRegion), size = 0) +
      scale_color_manual(values = c("No"=NA,"Yes"="red"), guide="none")+
      scale_alpha_continuous(range=c(0.2,1), guide="none")+
      scale_fill_gradientn(name="LifeExpectany",colors=Mypal, 
                           labels = function(x) {round(x,0)})+
      ggtitle("Male")+
      theme_void() + 
      theme(plot.title = element_text(size = 15, face = "bold"),
            legend.position = "bottom") 
    
    #adding both
    cowplot::plot_grid(g1F, g1M)
    
    
  })
  
  #2. Table Output
  output$LifeTable <- renderTable({
    Table <- bind_rows(#Female
      LifeExpFrameFemaleStacking %>% filter(RegNumber==input$Region & Year == input$Year) %>%
        slice(1:2) %>% #only first two rows
        select(1:5) %>% #only first 5 columns
        pivot_wider(id_cols = 1:2, names_from = 5, values_from = 3:4) %>% #pivot wide
        select(c(1:4,6,5)) %>% #reorder
        mutate("Width 50%"=PiUp_0.5-PiLo_0.5,
               "Width 80%"=PiUp_0.8-PiLo_0.8), 
      #Male
      LifeExpFrameMaleStacking %>% filter(RegNumber== input$Region & Year == input$Year) %>%
        slice(1:2) %>% #only first two rows
        select(1:5) %>% #only first 5 columns
        pivot_wider(id_cols = 1:2, names_from = 5, values_from = 3:4) %>% #pivot wide
        select(c(1:4,6,5)) %>% #reorder
        mutate("Width 50%"=PiUp_0.5-PiLo_0.5,
               "Width 80%"=PiUp_0.8-PiLo_0.8)
    ) %>% mutate("Sex"=c("Female","Male"),.before=1) %>%
      rename(., "10%" = PiLo_0.8,
             "25%"=PiLo_0.5,
             "75%"= PiUp_0.5,
             "90%" = PiUp_0.8) %>% 
      mutate_if(is.numeric, ~round(., 2)) # round numeric columns (later transformed into character)
    
    #Transpose Table and save as tible
    # as_tibble(cbind(names(Table), t(Table))) %>% #transpose
    #   rename(.," " = V1,                     #rename columns
    #         "Female"=V2,
    #          "Male" =V3) %>%
    #   slice(-1)
  })
  
  #3. Output Line Plot of Life Expectancy
  output$YearPlot <- renderPlot({
    
    #Plot for Females
    PlotF <- LifeExpFrameFemaleStacking %>% 
      filter(RegNumber==input$Region) %>% #Differentiate Mean Color of men and woman
      mutate("WUnique"=paste0(Width, substring(Sex,1,1))) %>%  
      ggplot(data=.)+
      ggdist::geom_lineribbon(aes(x = Year, y = Mean,lty=Type, fill=WUnique, ymin=PiLo, ymax=PiUp),
                      col = "#400040",alpha=0.5,linewidth=0.8,group=1)+
      scale_fill_manual(values = c("0.8w"="#e5cce5","0.5w"="#bf7fbf",
                                   "0.8m"="#d1e1ec", "0.5m"="#b3cde0"))+ #values of ribbon (second value is 50%PI)
      geom_vline(xintercept = as.character(input$Year), lty = 4)+
      theme_bw()+
      ylab("Life Expectancy")+xlab("Year")+
      scale_x_discrete(breaks=c(2001,2010,2020,2030))+
      ylim(c(75,90))+
      theme(legend.position = "none",
            axis.text = element_text(size = 13), #change size of axis text
            axis.title = element_text(size=13), #change size of axis title
            panel.grid = element_line(colour = "grey92"))
    
    #Plot for Males
    PlotM <- LifeExpFrameMaleStacking %>% 
      filter(RegNumber==input$Region) %>% #Differentiate Mean Color of men and woman
      mutate("WUnique"=paste0(Width, substring(Sex,1,1))) %>%  
      ggplot(data=.)+
      ggdist::geom_lineribbon(aes(x = Year, y = Mean,lty=Type, fill=WUnique, ymin=PiLo, ymax=PiUp),
                      col = "#011f4b",alpha=0.5,linewidth=0.8,group=1)+
      scale_fill_manual(values = c("0.8w"="#e5cce5","0.5w"="#bf7fbf",
                                   "0.8m"="#d1e1ec", "0.5m"="#b3cde0"))+ #values of ribbon (second value is 50%PI)
      geom_vline(xintercept = as.character(input$Year), lty = 4)+
      theme_bw()+
      ylab("Life Expectancy")+xlab("Year")+
      scale_x_discrete(breaks=c(2001,2010,2020,2030))+
      ylim(c(75,90))+
      theme(legend.position = "none",
            axis.text = element_text(size = 13), #change size of axis text
            axis.title = element_text(size=13), #change size of axis title
            panel.grid = element_line(colour = "grey92"))
    
    #adding Both
    cowplot::plot_grid(PlotF, PlotM)
  })
  
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)
