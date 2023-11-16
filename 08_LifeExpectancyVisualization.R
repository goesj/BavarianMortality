## FILE FOR VISUALISATION OF RESULTS############################################
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2,RColorBrewer, cowplot, ggdist,
               ggridges, ggpubr)


### Load Functions
source("01_Functions.R")
#### Load Data ####
load("Data/TotalData.RData") #load Data

## Load Life Expectancy Stacking Frame Results (uploaded in Git)
load("Data/LifeExpFrameStacking.RData")

################################################################################
##### 1) Plot Mean LE Estimates for all Regions in 2017 in Map #################
################################################################################
InSampleData <- TrainingDataFun(TotalData, sex="male",LastYearObs = 2017)

#Save Results in Bavarian Data Frame
Bayern$LifeExpStack_Female<- LifeExpFrameFemaleStacking %>% 
                                          filter(Year == 2017 & 
                                                 Type == "InSample" & 
                                                 Width == 0.8) %>%
                                          select(Mean) %>% pull()
#calculation of Rank with highest value being 1                                                            
Bayern$RankLifeExpF <- rank(-Bayern$LifeExpStack_Female, 
                            ties.method = "random") 


#Save Results in Bavarian Data (for males)
Bayern$LifeExpStack_Male<- LifeExpFrameMaleStacking %>% 
                                        filter(Year == 2017 & 
                                                 Type == "InSample" & 
                                                 Width == 0.8) %>%
                                        select(Mean) %>% pull()
#Rank with highest value being 1
Bayern$RankLifeExpM <- rank(-Bayern$LifeExpStack_Male, 
                            ties.method = "random") 

#Create own Color Palette## For Females ###
Mypal <- RColorBrewer::brewer.pal(n = 11, name = "RdBu") 

#Create MAP first
g1F <- ggplot(data = Bayern) +
  geom_sf(aes(fill = LifeExpStack_Female), alpha = 0.9, 
          colour = "transparent", size = 0) +
  geom_sf_text(aes(label=RankLifeExpF, 
                   color=LifeExpStack_Female <83 | LifeExpStack_Female > 84.5 ),
               fontface=2, size=4.5)+
  scale_color_manual(guide = FALSE, values = c("black", "white"))+ 
  scale_fill_gradientn(colors=Mypal)+
  ggtitle("a. Female Life Expectancy")+
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_text(size = 25,hjust=0.5, face="bold")) 

# Crarte Kernel Density Estimate of LE
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


g3F <- ggpubr::get_legend(g2F) #extract the legend 
g2FDash <- g2F+ theme(legend.position = "none") #remove legend
g3Fplot <- ggpubr::as_ggplot(g3F) #legend as own plot


#for single picture, change font size text =3, theme text =15 in g1, 
#legend text = 10 in g2
#Add Map, Kernel Density and Legend together
gFemale <- 
  cowplot::ggdraw() +
  cowplot::draw_plot(g1F)+
  cowplot::draw_plot(g2FDash, x = 0.71, y = .8, width = .29, height = .15)+
  cowplot::draw_plot(g3Fplot, x = 0.80, y = .715, width = .1, height = .1)


### Same Plot for males
g1M <- ggplot(data = Bayern) +
  geom_sf(aes(fill = LifeExpStack_Male), alpha = 0.9, colour = "transparent", 
          size = 0) +
  geom_sf_text(aes(label=RankLifeExpM, 
                   color=LifeExpStack_Male <78.1 | LifeExpStack_Male > 80.5 ),
               fontface=2, size=4.5)+
  scale_color_manual(guide = FALSE, values = c("black", "white"))+ 
  scale_fill_gradientn(colors=Mypal)+
  ggtitle("b. Male Life Expectancy")+
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_text(size = 25, hjust=0.5, face="bold")) 


#Create Kernel Density Estimate
BayernDataM <- data.frame("Mean"=Bayern$LifeExpStack_Male,
                          "Group"=1)

g2M <- ggplot(BayernDataM, aes(x = Mean, y = Group, fill = stat(x))) +
  ggridges::geom_density_ridges_gradient(scale = 3, rel_min_height = 0.001) +
  scale_fill_gradientn(colours = Mypal)+
  theme_void()+
  theme(legend.position = "bottom",
        legend.title=element_blank(),
        legend.key.width = unit(0.8, 'cm'),
        legend.text=element_text(face = "bold", size=15))

g3M <- ggpubr::get_legend(g2M) #extract the legend 
g2MDash <- g2M+ theme(legend.position = "none") #remove legend
g3Mplot <- ggpubr::as_ggplot(g3M) #legend as own plot



#for single picture, change font size text =3, theme text =15 in g1, 
#legend text = 10 in g2
gMale <- ggdraw() +
  cowplot::draw_plot(g1M)+
  cowplot::draw_plot(g2MDash, x = 0.71, y = .8, width = .29, height = .15)+
  cowplot::draw_plot(g3Mplot, x = 0.80, y = .715, width = .1, height = .1) 


x11()
cowplot::plot_grid(gFemale, gMale) #save/look at in 16 to 9 for proper format




################################################################################
##### 2) Plot Difference in LE for Females and Males between 2001  and 2032 ####
################################################################################
DiffExpFemale <- (LifeExpFrameFemaleStacking %>% 
                    filter(Year == 2030 & #select Year
                             Width == 0.8) %>% #get only one width  
                    select(Mean) %>% pull()) - 
                  (LifeExpFrameFemaleStacking %>% 
                        filter(Year == 2001 & 
                             Width == 0.8) %>% #get only one width  
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
  scale_fill_brewer(limits=labels, palette = "Reds")+ 
  ggtitle("a. Females")+
  theme_void()+
  theme(legend.position = "none",
        plot.title = element_text(size = 13, face = "bold", hjust=0.5)) 

#plot for Females
gMDiff <- Bayern %>% mutate("DiffM"= cut(DiffExpMale,
                                         breaks=QB, 
                                         labels=labels,
                                         include.lowest=F)) %>% 
  ggplot(data=.)+
  geom_sf(aes(fill=DiffM),alpha = 0.9, colour = "transparent", size = 0)+
  scale_fill_brewer(limits=labels, palette = "Reds", 
                    name="")+
  ggtitle("b. Males")+
  theme_void()+
  theme(legend.position = "right",
        plot.title = element_text(size = 13, face = "bold", hjust=0.5))

#extract legend
gLeg <- ggpubr::get_legend(gMDiff)
gMDiffNoLeg <- gMDiff + theme(legend.position = "none")


### plot both together
cowplot::ggdraw() +
  cowplot::draw_plot(gFDiff, x = 0.01, y = 0.1, width = 0.45, height = 0.95)+
  cowplot::draw_plot(gMDiffNoLeg, x = 0.39, y = .1, width=0.45, height=0.95)+
  cowplot::draw_plot(gLeg, x = 0.68, y = 0.1, width = 0.45, height = 0.95)



################################################################################
########## 3.) PLOT LE Predicitons of Single Region ############################
################################################################################
#Which Region to use 
RegUsed <- "Bamberg (Krfr.St)"
NumUsed <- InSampleData$Data %>% 
  filter(RegionName==RegUsed) %>% 
  select(RegionNumber) %>% 
  head(., 1) %>% pull()

#Plot
#combine by rbind in and Out of Sample
rbind(LifeExpFrameFemaleStacking, #combine males and females
      LifeExpFrameMaleStacking) %>% 
  filter(RegNumber==NumUsed) %>%
  mutate("WUnique"=paste0(Width, substring(Sex,1,1))) %>% #Differentiate Mean Color of men and woman
  ggplot(., aes(group=Sex))+
  ggdist::geom_lineribbon(aes(x = Year, y = Mean,lty=Type,fill=WUnique, ymin=PiLo, ymax=PiUp,
                      col=Sex),
                  alpha=0.5,size=0.8)+
  scale_color_manual(values=c("weiblich"="#400040",
                              "männlich"="#011f4b"))+
  scale_fill_manual(values = c("0.8w"="#e5cce5","0.5w"="#bf7fbf",
                               "0.8m"="#d1e1ec", "0.5m"="#b3cde0"))+ #values of ribbon (second value is 50%PI)
  theme_bw()+
  ylab("Life Expectancy at Birth")+xlab("Year")+
  scale_x_discrete(breaks=c(2001,2010,2020,2030))+
  ggtitle("Life Expectancy Predictions for Bamberg")+
    theme(legend.position = "none",
        axis.text = element_text(size = 25), #change size of axis text
        axis.title = element_text(size=25 , face="bold"), #change size of axis title
        panel.grid = element_line(colour = "grey92"),
        plot.title = element_text(hjust=0.5, size=25, face="bold"),
        panel.border = element_blank())




