##Load Libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2,RColorBrewer, cowplot, ggdist)


### Load Functions
source("01_Functions.R")
#### Load Data ####
load("Data/TotalData.RData") #load Data

## Load Life Expectancy Stacking Frame Results
load("Data/LifeExpFrameStacking.RData")

## Map Plot Year 17 (next to each other) ##
InSampleData <- TestDataFun(TotalData, Sex="male",LastYearObs = 2017)

#Save Results in Bavarian Data
Bayern$LifeExpStack_Female<- LifeExpFrameFemaleStacking %>% 
                                          filter(Year == 2017 & 
                                                 Type == "InSample" & 
                                                 Width == 0.8) %>%
                                          select(Mean) %>% pull()
                                                              
Bayern$RankLifeExpF <- rank(-Bayern$LifeExpStack_Female, ties.method = "random") #Rank with highest value being 1


#Save Results in Bavarian Data
Bayern$LifeExpStack_Male<- LifeExpFrameMaleStacking %>% 
                                        filter(Year == 2017 & 
                                                 Type == "InSample" & 
                                                 Width == 0.8) %>%
                                        select(Mean) %>% pull()

Bayern$RankLifeExpM <- rank(-Bayern$LifeExpStack_Male, ties.method = "random") #Rank with highest value being 1


Mypal <- brewer.pal(n = 11, name = "RdBu") #Create own Color Palette## For Females ###
#Creation of Plot in Text
g1F <- ggplot(data = Bayern) +
  geom_sf(aes(fill = LifeExpStack_Female), alpha = 0.9, colour = "transparent", size = 0) +
  geom_sf_text(aes(label=RankLifeExpF, color=LifeExpStack_Female <83 | LifeExpStack_Female > 84.5 ),fontface=2, size=4.5)+
  scale_color_manual(guide = FALSE, values = c("black", "white"))+ #plot color depending on value see (https://stackoverflow.com/questions/47281365/text-color-based-on-contrast-against-background)
  scale_fill_gradientn(colors=Mypal)+
  ggtitle("Female")+
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_text(size = 30, face = "bold")) 


#Create Dataframe for ggplot object
BayernData <- data.frame("Mean"=Bayern$LifeExpStack_Female,
                         "Group"=1)

g2F <- ggplot(BayernData, aes(x = Mean, y = Group, fill = stat(x))) +
  ggridges::geom_density_ridges_gradient(scale = 3, rel_min_height = 0.001) +
  scale_fill_gradientn(colours = Mypal)+
  theme_void()+
  theme(legend.position = "bottom",
        legend.title=element_blank(),
        legend.key.width = unit(0.8, 'cm'),
        legend.text=element_text(face = "bold", size=15))

g3F <- get_legend(g2F) #extract the legend 
g2FDash <- g2F+ theme(legend.position = "none") #remove legend
g3Fplot <- as_ggplot(g3F) #legend as own plot
  
x11()
#for single picture, change font size text =3, theme text =15 in g1, legend text = 10 in g2
gFemale <- 
  ggdraw() +
  draw_plot(g1F)+
  draw_plot(g2FDash, x = 0.71, y = .8, width = .29, height = .15)+
  draw_plot(g3Fplot, x = 0.80, y = .715, width = .1, height = .1)


### FOR MALES ###
g1M <- ggplot(data = Bayern) +
    geom_sf(aes(fill = LifeExpStack_Male), alpha = 0.9, colour = "transparent", size = 0) +
    geom_sf_text(aes(label=RankLifeExpM, color=LifeExpStack_Male <78.1 | LifeExpStack_Male > 80.5 ),fontface=2, size=4.5)+
    scale_color_manual(guide = FALSE, values = c("black", "white"))+ #plot color depending on value see (https://stackoverflow.com/questions/47281365/text-color-based-on-contrast-against-background)
    scale_fill_gradientn(colors=Mypal)+
    ggtitle("Male")+
    theme_void() + 
    theme(legend.position = "none",
          plot.title = element_text(size = 30, face = "bold")) 
  
  
#Create Dataframe for ggplot object
BayernDataM <- data.frame("Mean"=Bayern$LifeExpStack_Male,
                           "Group"=1)
  
g2M <- ggplot(BayernDataM, aes(x = Mean, y = Group, fill = stat(x))) +
    geom_density_ridges_gradient(scale = 3, rel_min_height = 0.001) +
    scale_fill_gradientn(colours = Mypal)+
    theme_void()+
    theme(legend.position = "bottom",
          legend.title=element_blank(),
          legend.key.width = unit(0.8, 'cm'),
          legend.text=element_text(face = "bold", size=15))
  
g3M <- get_legend(g2M) #extract the legend 
g2MDash <- g2M+ theme(legend.position = "none") #remove legend
g3Mplot <- as_ggplot(g3M) #legend as own plot
  
x11()
#for single picture, change font size text =3, theme text =15 in g1, legend text = 10 in g2
#pdf(file=file.path(getwd(),"Code/Bilder/LifeExpMale_17_Stack.pdf"))
gMale <- ggdraw() +
  draw_plot(g1M)+
  draw_plot(g2MDash, x = 0.71, y = .8, width = .29, height = .15)+
  draw_plot(g3Mplot, x = 0.80, y = .715, width = .1, height = .1) 


x11()
plot_grid(gFemale, gMale) #save/look at in 16 to 9 for proper format




## Difference in LE for Females and Males between 2001  and 2032 ##
DiffExpFemale <- (LifeExpFrameFemaleStacking %>% 
                  filter(Year == 2030 & 
                         Width == 0.8) %>% #get only one width as they are double 
                  select(Mean) %>% pull()) - 
                (LifeExpFrameFemaleStacking %>% 
                   filter(Year == 2001 & 
                          Width == 0.8) %>% #get only one width as they are double 
                select(Mean) %>% pull())

DiffExpMale <- (LifeExpFrameMaleStacking %>% 
                  filter(Year == 2030 & 
                         Width == 0.8) %>% #get only one width  
                  select(Mean) %>% pull()) - 
            (LifeExpFrameMaleStacking %>% 
                filter(Year == 2001 & 
                      Width == 0.8) %>% 
                select(Mean) %>% pull())

##Calculate Quantiles for Breaks (both Males and Femals)
QB <- quantile(c(DiffExpFemale, DiffExpMale),probs=seq(0,1,1/9)) %>% round(.,2)
#Get Labels
labels <- paste0(QB[1:9],"-",
                 QB[2:10])

#plot for Females
gFDiff <- Bayern %>% mutate("DiffF"= cut(DiffExpFemale,
                                     breaks=QB, 
                                     labels=labels,
                                     include.lowest=F)) %>% 
  ggplot(data=.)+
  geom_sf(aes(fill=DiffF),alpha = 0.9, colour = "transparent", size = 0)+
  scale_fill_brewer(limits=labels, palette = "YlGnBu")+
  ggtitle("Female")+
  theme_void()+
  theme(legend.position = "none",
        plot.title = element_text(size = 15, face = "bold")) 

#plot for Females
gMDiff <- Bayern %>% mutate("DiffM"= cut(DiffExpMale,
                                     breaks=QB, 
                                     labels=labels,
                                     include.lowest=F)) %>% 
  ggplot(data=.)+
  geom_sf(aes(fill=DiffM),alpha = 0.9, colour = "transparent", size = 0)+
  scale_fill_brewer(limits=labels, palette = "YlGnBu",
                    name="")+
  ggtitle("Male")+
  theme_void()+
  theme(legend.position = "right",
        plot.title = element_text(size = 15, face = "bold"))

#extract legend
gLeg <- get_legend(gMDiff)
gMDiffNoLeg <- gMDiff + theme(legend.position = "none")


x11()
ggdraw() +
  draw_plot(gFDiff, x = 0.01, y = 0.1, width = 0.45, height = 0.95)+
  draw_plot(gMDiffNoLeg, x = 0.39, y = .1, width=0.45, height=0.95)+
  draw_plot(gLeg, x = 0.68, y = 0.1, width = 0.45, height = 0.95)



######## PLOT LE Estrimate of Region#######
#Which Region to use
RegUsed <- "Bamberg (Krfr.St)"
NumUsed <- InSampleData$Data$RegionNumber[which(InSampleData$Data$RegionName == RegUsed)[1]]

#Plot
x11() 
#combine by rbind in and Out of Sample
  rbind(LifeExpFrameFemaleStacking, #combine males and females
        LifeExpFrameMaleStacking) %>% 
  filter(RegNumber==NumUsed) %>%
    mutate("WUnique"=paste0(Width, substring(Sex,1,1))) %>% #Differentiate Mean Color of men and woman
    ggplot(., aes(group=Sex))+
    geom_lineribbon(aes(x = Year, y = Mean,lty=Type,fill=WUnique, ymin=PiLo, ymax=PiUp,
                        col=Sex),
                    alpha=0.5,size=0.8)+
  scale_color_manual(values=c("female"="#400040",
                              "male"="#011f4b"))+
  scale_fill_manual(values = c("0.8w"="#e5cce5","0.5w"="#bf7fbf",
                               "0.8m"="#d1e1ec", "0.5m"="#b3cde0"))+ #values of ribbon (second value is 50%PI)
  theme_bw()+
  ylab("Life Expectancy")+xlab("Year")+
  scale_x_discrete(breaks=c(2001,2010,2020,2030))+
  theme(legend.position = "none",
        axis.text = element_text(size = 13), #change size of axis text
        axis.title = element_text(size=13), #change size of axis title
        panel.grid = element_line(colour = "grey92"))





