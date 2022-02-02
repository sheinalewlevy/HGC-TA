##This is the code for data preparation for the analysis of foraging children's participation in play and work. 
##NOTE: This code will not run without requesting the Tsimane dataset from Jonathan Stieglitz (jonathan.stieglitz@iast.fr)

##Make sure to set working diretory to 'data'

##IF you don't have the necessary packages, load them following these commands
##install.packages("dplyr")
##install.packages("tidyr")
##install.packages("reshape2")
##install.packages("splitstackshape")
##install.packages("ggplot2")
##install.packages("data.table")
##install.packages("ggpubr")

##Load necessary packages

library(dplyr)
library(tidyr)
library(reshape2)
library(splitstackshape)
library(ggplot2)
library(data.table)
library(ggpubr)
library(lme4)
library(usdm)

###################
##IMPORT DATA#####
##################

Agta_data <- read.csv("Agta_Hagen.csv")
Aka_data<-read.csv("Aka_Boyette.csv")
Baka_data<-read.csv("Baka_Sonoda.csv")
BaYaka_data<-read.csv("BaYaka_Lew-Levy.csv")
Dukha_data<-read.csv("Dukha_Randall_Obrien_Surovell.csv")
Hadza_data<-read.csv("Hadza_Crittenden.csv")
Camana_data<-read.csv("Matsigenka_Camana_Baksh.csv")
Shimaa_data<-read.csv("Matsigenka_Shimaa_Johnson.csv")
Maya_data<-read.csv("Maya_Kramer.csv")
Mayangna_data<-read.csv("Mayangna_Koster.csv")
Mikea_data<-read.csv("Mikea_Tucker.csv")
Pume_data<-read.csv("Pume_Kramer.csv")
Tsimane_data<-read.csv("Tsimane_Stieglitz.csv") ##This dataset is available by request from Jonathan Stieglitz

##Now we import the dataset with the group-level variables. For details on how div labor and prop non-foraged food was calculated, see end of R script (section entitled: additional information). For the environmental variables, see Enviro_variables.R in the supplementary materials.
group<-read.csv("group_variables.csv")

##Bind the rows together of the site datasets
data_full<-bind_rows(Agta_data,Aka_data,Baka_data,BaYaka_data,Camana_data,Dukha_data,Hadza_data,Maya_data,Mayangna_data,Mikea_data,Pume_data,Shimaa_data,Tsimane_data,.id=NULL)

##Now we make a new variable for age in round numbers. We choose 'floor' instead of 'round' because this is how age is usually reported
data_full$newAge<-floor(data_full$Age)

##Make sure the range is correct (3-20)
range(data_full$newAge)

##Double check that all the rows were included properly
sum(data_full$offset)- sum(Agta_data$offset,Aka_data$offset,Baka_data$offset,BaYaka_data$offset,Camana_data$offset,Dukha_data$offset,Hadza_data$offset,Maya_data$offset,Mayangna_data$offset,Mikea_data$offset,Pume_data$offset,Shimaa_data$offset,Tsimane_data$offset)##should equal 0

######################
###DATA DESCRIPTION###
######################

##Now we produce the information for the first part of Table S1. This includes number of kids (number of repeated kids), % female, Mean Age (SD), and mean number of observations/per child (SD). In datasets where children were observed >1 time, we pulled unique sex values, meaned age, and summed number of observations per individual for the means table. In cases where >1 observation could be recorded at a time, we also keep track of the percentage of observations reported.


##Agta
length(unique(Agta_data$newID))##15 kids, no repeats
round(table(Agta_data$Sex)/sum(table(Agta_data$Sex)),2) ##33% girls
round(mean(floor(Agta_data$Age)),digits=2) ##mean age=6.13 
round(sd(floor(Agta_data$Age)),digits=2) ##sd age=2.64
range(floor(Agta_data$Age)) ##3-12
round(mean(Agta_data$offset),digits=2) ##21.27 obs
round(sd(Agta_data$offset),digits=2) ##11.32 obs
(sum(Agta_data$doubled)/sum(Agta_data$offset))*100 ##1% of observations coded concurrently


##Aka
length(unique(Aka_data$newID))##50 kids, no repeats
round(table(Aka_data$Sex)/sum(table(Aka_data$Sex)),2) ##52% female
round(mean(floor(Aka_data$Age)),digits=2) ##mean age=9.44
round(sd(floor(Aka_data$Age)),digits=2) ##sd age=3.89
range(floor(Aka_data$Age)) ##4-16
round(mean(Aka_data$offset),digits=2) ##238.62 obs
round(sd(Aka_data$offset),digits=2) ##53.49 obs


##Baka
length(unique(Baka_data$newID))##14 kids, no repeats
round(table(Baka_data$Sex)/sum(table(Baka_data$Sex)),2) ##50% female
round(mean(floor(Baka_data$Age)),digits=2) ##mean age=9.21
round(sd(floor(Baka_data$Age)),digits=2) ##sd age=3.62
range(floor(Baka_data$Age)) ##5-15
round(mean(Baka_data$offset),digits=2) ##720 obs
round(sd(Baka_data$offset),digits=2) ##0 obs
(sum(Baka_data$doubled)/sum(Baka_data$offset))*100 ##7% of observations coded concurrently


##BaYaka
length(unique(BaYaka_data$newID))##53 kids, yes repeats
n_occur <- data.frame(table(BaYaka_data$newID))
BaYaka_rep<-n_occur[n_occur$Freq > 1,] 
length(unique(BaYaka_rep$Var1))##6 repeats

BaYaka_demo<-BaYaka_data%>%
  group_by(newID)%>%
  summarise(Sex=unique(Sex),Age=mean(floor(Age)),offset=sum(offset))
round(table(BaYaka_demo$Sex)/sum(table(BaYaka_demo$Sex)),2)##42% female

round(mean(BaYaka_demo$Age),digits=2) ##mean age=11.02
round(sd(BaYaka_demo$Age),digits=2) ##sd age=4.17
range(floor(BaYaka_data$Age))
round(mean(BaYaka_demo$offset),digits=2) ##253.87 obs
round(sd(BaYaka_demo$offset),digits=2) ##89.44 obs
(sum(BaYaka_data$doubled)/sum(BaYaka_data$offset))*100 ##1% of  observations coded concurrently

##Dukha
length(unique(Dukha_data$newID))##15 kids, yes repeats
n_occur <- data.frame(table(Dukha_data$newID))
Dukha_rep<-n_occur[n_occur$Freq > 1,] ##5 repeats
length(unique(Dukha_rep$Var1))
Dukha_demo<-Dukha_data%>%
  group_by(newID)%>%
  summarise(Sex=unique(Sex),Age=mean(floor(Age)),offset=sum(offset))
round(table(Dukha_demo$Sex)/sum(table(Dukha_demo$Sex)),2) ##53% female
round(mean(Dukha_demo$Age),digits=2) ##mean age=9.17
round(sd(Dukha_demo$Age),digits=2) ##sd age=5.57
range(floor(Dukha_data$Age)) ##3-18
round(mean(Dukha_demo$offset),digits=2) ##577 obs
round(sd(Dukha_demo$offset),digits=2) ##562.25 obs


##Hadza
length(unique(Hadza_data$newID)) ##18 kids, no repeats
round(table(Hadza_data$Sex)/sum(table(Hadza_data$Sex)),2) ##78% female
round(mean(floor(Hadza_data$Age)),digits=2) ##8.39
round(sd(floor(Hadza_data$Age)),digits=2) ##3.11
range(floor(Hadza_data$Age)) ##3-14
round(mean(Hadza_data$offset),digits=2) ##35.72 obs
round(sd(Hadza_data$offset),digits=2) ##26.81 obs
(sum(Hadza_data$doubled)/sum(Hadza_data$offset))*100 ##6% of observations coded concurrently

##Matsigenka
Matsigenka_data<-subset(data_full,Society=="Matsigenka")
length(unique(Matsigenka_data$newID)) ##119, no repeats
round(table(Matsigenka_data$Sex)/sum(table(Matsigenka_data$Sex)),2) ##61% female
round(mean(floor(Matsigenka_data$Age)),digits=2) ##9.28
round(sd(floor(Matsigenka_data$Age)),digits=2) ##4.62
range(floor(Matsigenka_data$Age)) ##3-18
round(mean(Matsigenka_data$offset),digits=2) ##24.22 obs
round(sd(Matsigenka_data$offset),digits=2) ##11.64 obs

##Maya
length(unique(Maya_data$newID)) ##49 kids, no repeats
round(table(Maya_data$Sex)/sum(table(Maya_data$Sex)),2)##59% female
round(mean(floor(Maya_data$Age)),digits=2) ##9.47
round(sd(floor(Maya_data$Age)),digits=2) ##4.98
range(floor(Maya_data$Age)) ##3-18
round(mean(Maya_data$offset),digits=2)##149.14 obs
round(sd(Maya_data$offset),digits=2) ##16.69 obs
sum(Maya_data$doubled)/sum(Maya_data$offset) ##2% of observations coded concurrently


##Mayangna
length(unique(Mayangna_data$newID)) ##114 kids, no repeats
round(table(Mayangna_data$Sex)/sum(table(Mayangna_data$Sex)),2)##46% female
round(mean(floor(Mayangna_data$Age)),digits=2) ##9.61
round(sd(floor(Mayangna_data$Age)),digits=2) ##4.89
range(floor(Mayangna_data$Age)) ##3-18
round(mean(Mayangna_data$offset),digits=2) ##67.56 obs
round(sd(Mayangna_data$offset),digits=2) ##17.86 obs


##Mikea
length(unique(Mikea_data$newID)) ##31 kids, yes repeats
n_occur <- data.frame(table(Mikea_data$newID))
Mikea_rep<-n_occur[n_occur$Freq > 1,] ##18 repeats
length(unique(Mikea_rep$Var1))
Mikea_demo<-Mikea_data%>%
  group_by(newID)%>%
  summarise(Sex=unique(Sex),Age=mean(Age),offset=sum(offset))
round(table(Mikea_demo$Sex)/sum(table(Mikea_demo$Sex)),2)##48% female
round(mean(Mikea_demo$Age),digits=2) ##mean age=11.51
round(sd(Mikea_demo$Age),digits=2) ##sd age=3.78
round(mean(Mikea_demo$offset),digits=2) ##150.13 obs
round(sd(Mikea_demo$offset),digits=2) ##118.67 obs


##Pume
length(unique(Pume_data$newID)) ##31 kids, yes repeats
n_occur <- data.frame(table(Pume_data$newID))
Pume_rep<-n_occur[n_occur$Freq > 1,] ##26 repeats
length(unique(Pume_rep$Var1))
Pume_demo<-Pume_data%>%
  group_by(newID)%>%
  summarise(Sex=unique(Sex),Age=mean(floor(Age)),offset=sum(offset))
round(table(Pume_demo$Sex)/sum(table(Pume_demo$Sex)),2)##52% female
round(mean(Pume_demo$Age),digits=2) ##mean age=9.32
round(sd(Pume_demo$Age),digits=2) ##sd age=4.41
range(floor(Pume_data$Age)) ##3-17
round(mean(Pume_demo$offset),digits=2) ##166.87 obs
round(sd(Pume_demo$offset),digits=2) ##81.44 obs


##Tsimane--this dataset is available by request to Jonathan Stieglitz
length(unique(Tsimane_data$newID)) ##181, no repeats
round(table(Tsimane_data$Sex)/sum(table(Tsimane_data$Sex)),2)##52% female
round(mean(floor(Tsimane_data$Age)),digits=2) ##8.53
round(sd(floor(Tsimane_data$Age)),digits=2) ##4.3
range(floor(Tsimane_data$Age)) ##3-18
round(mean(Tsimane_data$offset),digits=2) ##70.69 obs
round(sd(Tsimane_data$offset),digits=2) ##20.76 obs
(sum(Tsimane_data$doubled)/sum(Tsimane_data$offset))*100 ##1% of observations coded concurrently


##And for the full dataset 
length(unique(data_full$newID)) ##690 kids total
n_occur <- data.frame(table(data_full$newID))
full_rep<-n_occur[n_occur$Freq > 1,] ##55 repeats
length(unique(full_rep$Var1))
full_demo<-data_full%>%
  group_by(newID)%>%
  summarise(Sex=unique(Sex),Age=mean(newAge),offset=sum(offset))
round(table(full_demo$Sex)/sum(table(full_demo$Sex)),2)##52% female
round(mean(full_demo$Age),digits=2) ##mean age=9.29
round(sd(full_demo$Age),digits=2) ##sd age=4.48
round(mean(full_demo$offset),digits=2)##124.05 obs
round(sd(full_demo$offset),digits=2) ##160.88 obs

sum(full_demo$offset)##85597

(sum(Hadza_data$doubled,BaYaka_data$doubled,Maya_data$doubled,Tsimane_data$doubled,Agta_data$doubled,Baka_data$doubled)/sum(Hadza_data$offset,BaYaka_data$offset,Maya_data$offset,Tsimane_data$offset,Agta_data$offset,Baka_data$offset))*100 ##2.6% of those datasets which could could concurrently involved doubled codes

##############
###FIGURE 2###
##############

##To make Figure 2 from the text, first make the age groups
data_full$Age_cat<-ifelse(data_full$Age>12,3,ifelse(data_full$Age<7,1,2))

 ###we'll sum all variables by ID and age group
figure<-data_full%>%
  group_by(newID,Age_cat)%>%
  summarise(Society=unique(Society),Sex=unique(Sex),play=sum(play),food_production=sum(food_production),childcare=sum(childcare),household=sum(household),offset=sum(offset))

##We'll give sex names for the figure
figure$Sex_title<-ifelse(figure$Sex==0,"Girls","Boys")

##And give names to the age categories
figure$Age_cat_title<-ifelse(figure$Age_cat==1,"Early Childhood (3-6 years)",ifelse(figure$Age_cat==2,"Middle Childhood (7-12 years)","Adolescence (13-18 years)"))

##Now we'll make the proportions
figure$propplay<-figure$play/figure$offset
figure$propfood<-figure$food_production/figure$offset
figure$propchildcare<-figure$childcare/figure$offset
figure$prophousehold<-figure$household/figure$offset

##And now we make means across the proportions by soxiety, age, and sex
figure_Society<-figure%>%
  group_by(Society,Age_cat_title,Sex_title)%>%
  summarise(play=mean(propplay),food_production=mean(propfood),childcare=mean(propchildcare),household=mean(prophousehold))

##And we'll make one for the total sample as well, to serve as a reference, and bind them together
figuretotal<-figure%>%
  group_by(Sex_title,Age_cat_title)%>%
  summarise(play=mean(propplay),food_production=mean(propfood),childcare=mean(propchildcare),household=mean(prophousehold))
figuretotal$Society<-"Total"
figure2<-rbind(figure_Society,figuretotal)

##We'll releven Age categories and Societies so that age is in the proper order from oldest to youngest, and so that the totals appear last in the figure
figure2$Age_cat_title<-ordered(figure2$Age_cat_title,levels=c("Early Childhood (3-6 years)","Middle Childhood (7-12 years)","Adolescence (13-18 years)"))
figure2$Society<-ordered(figure2$Society,levels=c("Agta", "Aka", "Baka", "BaYaka", "Dukha", "Hadza", "Matsigenka","Maya", "Mayangna", "Mikea", "Pume", "Tsimane", "Total"))
levels(figure2$Society)[levels(figure2$Society)=="Pume"] <- "Pumé"

##Now we make Figure 2
childcare<-ggplot(figure2, aes(x=Society, y=childcare, fill=Age_cat_title)) + geom_bar(stat="identity",position="dodge",width=0.7) + facet_grid(Sex_title~.)
childcare<- childcare + theme_classic() + ylab("Childcare") + xlab("") + theme(legend.title = element_blank())
childcare

food<-ggplot(figure2, aes(x=Society, y=food_production, fill=Age_cat_title)) + geom_bar(stat="identity",position="dodge",width=0.7) + facet_grid(Sex_title~.)
food<- food+ theme_classic() + ylab(" Food production") + xlab("") 
food

household<-ggplot(figure2, aes(x=Society, y=household, fill=Age_cat_title)) + geom_bar(stat="identity",position="dodge",width=0.7) +  facet_grid(Sex_title~.) 
household<- household + theme_classic() + ylab("Domestic work") + xlab("") 
household

play<-ggplot(figure2, aes(x=Society, y=play, fill=Age_cat_title)) + geom_bar(stat="identity",position="dodge",width=0.7) + facet_grid(Sex_title~.)
play<- play + theme_classic() + ylab("Play") + xlab("") 
play

fig2<-ggarrange(
  childcare, food,household,play, 
  common.legend = TRUE, legend = "right",
  nrow=4,ncol=1
)

##Save the figure straight to the workspace
annotate_figure(fig2,left=text_grob("Proportion of Time",rot=90))%>%
  ggexport(filename = "fig2.pdf",width=12,height=10)

###################################
###SUPPLEMENTARY ACTIVITY TABLES###
###################################

##Now we're going to save the values for Table S2
##Relevel so that early childhood comes first
figure$Age_cat_title<-ordered(figure$Age_cat_title,levels=c("Early Childhood (3-6 years)","Middle Childhood (7-12 years)","Adolescence (13-18 years)"))

TableS2a<-figure %>% 
  group_by(Society,Age_cat_title,Sex_title) %>% 
  summarise(
    play_sum=sum(play),
    food_production_sum=sum(food_production),
    childcare_sum=sum(childcare),
    household_sum=sum(household),
    offset_sum=sum(offset),
    play=round(play_sum/offset_sum,digits=2),
    play_lwr=round((binom.test(play_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    play_upr=round((binom.test(play_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2),
    food_production=round(food_production_sum/offset_sum,digits=2),
    food_production_lwr=round((binom.test(food_production_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    food_production_upr=round((binom.test(food_production_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2),
    childcare=round(childcare_sum/offset_sum,digits=2),
    childcare_lwr=round((binom.test(childcare_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    childcare_upr=round((binom.test(childcare_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2),
    household=round(household_sum/offset_sum,digits=2),
    household_lwr=round((binom.test(household_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    household_upr=round((binom.test(household_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2)
  ) %>% 
  dplyr::select(-play_sum, -food_production_sum, -childcare_sum, -household_sum, -offset_sum)

write.csv(TableS2a,"TableS2a.csv")

TableS2b<-figure%>%
  group_by(Age_cat_title,Sex_title)%>%
  summarise(
    play_sum=sum(play),
    food_production_sum=sum(food_production),
    childcare_sum=sum(childcare),
    household_sum=sum(household),
    offset_sum=sum(offset),
    play=round(play_sum/offset_sum,digits=2),
    play_lwr=round((binom.test(play_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    play_upr=round((binom.test(play_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2),
    food_production=round(food_production_sum/offset_sum,digits=2),
    food_production_lwr=round((binom.test(food_production_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    food_production_upr=round((binom.test(food_production_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2),
    childcare=round(childcare_sum/offset_sum,digits=2),
    childcare_lwr=round((binom.test(childcare_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    childcare_upr=round((binom.test(childcare_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2),
    household=round(household_sum/offset_sum,digits=2),
    household_lwr=round((binom.test(household_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    household_upr=round((binom.test(household_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2)
  ) %>% 
  dplyr::select(-play_sum, -food_production_sum, -childcare_sum, -household_sum, -offset_sum)

write.csv(TableS2b,"TableS2b.csv")

TableS3a<-figure%>%
  group_by(Society,Age_cat_title)%>%  summarise(
    play_sum=sum(play),
    food_production_sum=sum(food_production),
    childcare_sum=sum(childcare),
    household_sum=sum(household),
    offset_sum=sum(offset),
    play=round(play_sum/offset_sum,digits=2),
    play_lwr=round((binom.test(play_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    play_upr=round((binom.test(play_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2),
    food_production=round(food_production_sum/offset_sum,digits=2),
    food_production_lwr=round((binom.test(food_production_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    food_production_upr=round((binom.test(food_production_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2),
    childcare=round(childcare_sum/offset_sum,digits=2),
    childcare_lwr=round((binom.test(childcare_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    childcare_upr=round((binom.test(childcare_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2),
    household=round(household_sum/offset_sum,digits=2),
    household_lwr=round((binom.test(household_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    household_upr=round((binom.test(household_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2)
  ) %>% 
  dplyr::select(-play_sum, -food_production_sum, -childcare_sum, -household_sum, -offset_sum)
write.csv(TableS3a,"TableS3a.csv")

TableS3b<-figure%>%
  group_by(Age_cat_title)%>%  summarise(
    play_sum=sum(play),
    food_production_sum=sum(food_production),
    childcare_sum=sum(childcare),
    household_sum=sum(household),
    offset_sum=sum(offset),
    play=round(play_sum/offset_sum,digits=2),
    play_lwr=round((binom.test(play_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    play_upr=round((binom.test(play_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2),
    food_production=round(food_production_sum/offset_sum,digits=2),
    food_production_lwr=round((binom.test(food_production_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    food_production_upr=round((binom.test(food_production_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2),
    childcare=round(childcare_sum/offset_sum,digits=2),
    childcare_lwr=round((binom.test(childcare_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    childcare_upr=round((binom.test(childcare_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2),
    household=round(household_sum/offset_sum,digits=2),
    household_lwr=round((binom.test(household_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    household_upr=round((binom.test(household_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2)
  ) %>% 
  dplyr::select(-play_sum, -food_production_sum, -childcare_sum, -household_sum, -offset_sum)

write.csv(TableS3b,"TableS3b.csv")

TableS4a<-figure%>%
  group_by(Society,Sex_title)%>%  summarise(
    play_sum=sum(play),
    food_production_sum=sum(food_production),
    childcare_sum=sum(childcare),
    household_sum=sum(household),
    offset_sum=sum(offset),
    play=round(play_sum/offset_sum,digits=2),
    play_lwr=round((binom.test(play_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    play_upr=round((binom.test(play_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2),
    food_production=round(food_production_sum/offset_sum,digits=2),
    food_production_lwr=round((binom.test(food_production_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    food_production_upr=round((binom.test(food_production_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2),
    childcare=round(childcare_sum/offset_sum,digits=2),
    childcare_lwr=round((binom.test(childcare_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    childcare_upr=round((binom.test(childcare_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2),
    household=round(household_sum/offset_sum,digits=2),
    household_lwr=round((binom.test(household_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    household_upr=round((binom.test(household_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2)
  ) %>% 
  dplyr::select(-play_sum, -food_production_sum, -childcare_sum, -household_sum, -offset_sum)
write.csv(TableS4a,"TableS4a.csv")

TableS4b<-figure%>%
  group_by(Sex_title)%>%  summarise(
    play_sum=sum(play),
    food_production_sum=sum(food_production),
    childcare_sum=sum(childcare),
    household_sum=sum(household),
    offset_sum=sum(offset),
    play=round(play_sum/offset_sum,digits=2),
    play_lwr=round((binom.test(play_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    play_upr=round((binom.test(play_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2),
    food_production=round(food_production_sum/offset_sum,digits=2),
    food_production_lwr=round((binom.test(food_production_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    food_production_upr=round((binom.test(food_production_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2),
    childcare=round(childcare_sum/offset_sum,digits=2),
    childcare_lwr=round((binom.test(childcare_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    childcare_upr=round((binom.test(childcare_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2),
    household=round(household_sum/offset_sum,digits=2),
    household_lwr=round((binom.test(household_sum,offset_sum,conf.level=0.9)$conf.int)[1],digits=2),
    household_upr=round((binom.test(household_sum,offset_sum,conf.level=0.9)$conf.int)[2],digits=2)
  ) %>% 
  dplyr::select(-play_sum, -food_production_sum, -childcare_sum, -household_sum, -offset_sum)

write.csv(TableS4b,"TableS4b.csv")

######################
###ANALYSIS DATASET###
######################

##First, we sum all the activities to subtract from the offset variable. We add 'doubled' to nonworkobs to avoid having children with negative values for their number of observations. Since only a very small proportion of the datasets include doubled observations (0.5% of all observations) this should have minimal impact on the final results.

data_full$sumact<-data_full$food_production+data_full$household+data_full$play+data_full$childcare
data_full$doubled[is.na(data_full$doubled)]<-0
data_full$znonworkobs<-data_full$offset-data_full$sumact+data_full$doubled

myvars<-c("Society","Sex","food_production","household","play","childcare","znonworkobs","newAge","newID")

##Subset data
newdata<-data_full[myvars]

##Now we use the melt function which takes data in wide format and stacks a set of columns into a single column of data (which for us is variable=activity category and value=number of observations for that activity category). In other words, each child (per observation year) is represented 5 times, once for food production, once for household, once for play, once for childcare, and once for non work observations.
melt<-reshape2::melt(newdata,id=c("Society","newID","newAge","Sex"))

##And now we expand these rows so that there is as many rows for each activity category as there are values
melt2<-expandRows(melt, "value")

##To double check if the dataset melted correctly, there should be as many rows in melt2 as the following sum
sum(data_full$offset,data_full$doubled) ##86741 including the Tsimane, whose data is available upon request from Jonathan Stieglitz

##And we merge melt2 and group together

d<-merge(melt2,group,by="Society")

##Now make age categories for the analysis
d$Age_cat<-ifelse(d$newAge>12,3,ifelse(d$newAge<7,1,2))
d$Middle<-ifelse(d$Age_cat==2,1,0)
d$Ado<-ifelse(d$Age_cat==3,1,0)

##And now we write d into a saved dataset. 
write.csv(d,"dataset.csv")

##Check multicollinearity
d$dens<-ifelse(d$sum_density>1,1,0)
d$NPP_z<-(d$NPP-mean(d$NPP))/sd(d$NPP)
d$nonforaged_z<-(d$nonforaged-mean(d$nonforaged))/sd(d$nonforaged)
d$temp_z<-(d$meanAnnualTemp-mean(d$meanAnnualTemp))/sd(d$meanAnnualTemp)
d$prec_z<-(d$totalAnnualPrec-mean(d$totalAnnualPrec))/sd(d$totalAnnualPrec)

M3<-d%>%dplyr::select(Middle,Ado,Sex,div,nonforaged_z)
c1<-as.data.frame(cor(M3))
v1<-as.data.frame(vif(M3))
c1<-cbind(c1,v1)
write.csv(c1,"c1.csv")

M4<-d%>%dplyr::select(Middle,Ado,Sex,dens,water_rating,nonforaged_z)
c2<-as.data.frame(cor(M4))
v2<-as.data.frame(vif(M4))
c2<-cbind(c2,v2)
write.csv(c2,"c2.csv")

M5<-d%>%dplyr::select(Middle,Ado,Sex,NPP_z,temp_z,prec_z,nonforaged_z)
c3<-as.data.frame(cor(M5))
v3<-as.data.frame(vif(M5))
c3<-cbind(c3,v3)
write.csv(c3,"c3.csv")

#################################################################
###PROPORTION OF NON-FORAGED FOOD AND SEXUAL DIVISION OF LABOR###
#################################################################

###calculation of the non-foraged foods and the division of labour variables

activities <- c(
  "plant foods",
  "large game",
  "small game",
  "sea food",
  "insects",
  "honey",
  "farming",
  "trading"
);

Society <- c(
  "Aka",
  "Tsimane",
  "Baka",
  "Agta",
  "Matsigenka",
  "Mayangna",
  "Hadza",
  "BaYaka",
  "Dukha",
  "Pume",
  "Maya",
  "Mikea"
);

# Recorded time spent in each activity, lower bound, for each Society.
# If no range (i.e. just a single value), use the single recorded value.
# If no value is recorded, replace with 0.
percentsbySocietyLow <- data.frame(
  "Aka"=c(25, 0, 0, 5, 0, 0, 10, 50),
  "Tsimane"=c(0.7, 6, 1.85, 15, 0, 0.02, 61.9, 15),
  "Baka"=c(15, 5, 15, 20, 10, 5, 25, 5),
  "Agta"=c(10, 7, 3, 25, 0, 5, 10, 40),
  "Matsigenka"=c(0.83, 0, 0.23, 2.23, 0.21, 0, 96.5, 0),
  "Mayangna"=c(0, 5, 5, 7, 0, 1, 60, 10),
  "Hadza"=c(60, 18, 13, 0, 0, 12, 0, 5),
  "BaYaka"=c(20, 15, 15, 5, 10, 5, 30, 0),
  "Dukha"=c(2,10,0,1,0,0,40,40),
  "Pume"=c(100,1,5,10,1,0.5,9,0),
  "Maya"=c(2,1,1,0,0,2,92,2),
  "Mikea"=c(40,0,10,0,0,5,40,5),
  row.names=activities
)

# Recorded time spent in each activity, lower bound, for each Society.
# If no range (i.e. just a single value), use the single recorded value.
# If no value is recorded, replace with 0.
percentsbySocietyHigh <- data.frame(
  "Aka"=c(50, 5, 25, 5, 60, 25, 80, 50),
  "Tsimane"=c(0.7, 6, 1.85, 15, 0, 0.02, 61.9, 15),
  "Baka"=c(15, 5, 15, 20, 10, 5, 25, 5),
  "Agta"=c(10, 7, 3, 25, 0, 5, 10, 40),
  "Matsigenka"=c(0.83, 0, 0.23, 2.23, 0.21, 0, 96.5, 0),
  "Mayangna"=c(5, 5, 5, 7, 1, 1, 60, 10),
  "Hadza"=c(60, 18, 13, 0, 0, 12, 0, 5),
  "BaYaka"=c(20, 15, 15, 5, 10, 5, 30, 0),
  "Dukha"=c(2,10,0,1,0,0,40,40),
  "Pume"=c(100,1,5,10,1,0.5,9,0),
  "Maya"=c(2,1,1,0,0,2,92,2),
  "Mikea"=c(40,0,15,0,1,5,40,10),
  row.names=activities
)

# Recorded division of labour for each activity, for each Society.
divbySociety <- data.frame(
  "Aka"=c(2, 4, 4, 2, 3, 4, 3, 3),
  "Tsimane"=c(2, 5, 4,4 ,0,4,3,4),
  "Baka"=c(2, 5, 4, 2, 2, 4, 3, 3),
  "Agta"=c(2, 4, 3, 3, 0, 4, 3, 3),
  "Matsigenka"=c(3, 0, 5, 4, 3, 0, 4, 0),
  "Mayangna"=c(3, 5, 5, 3, 0, 5, 4, 4),
  "Hadza"=c(2, 5, 5, 0, 0, 4, 0, 4),
  "BaYaka"=c(2, 5, 4, 3, 2, 4, 3, 3),
  "Dukha"=c(3,5,5,5,0,0,3,2),
  "Pume"=c(2,5,4,4,4,4,3,0),
  "Maya"=c(3,3,5,0,0,5,4,5),
  "Mikea"=c(3 ,0,3,2,3,4,3,2),
  
  row.names=activities
)

# Calculate the mean percentage of time spent in each activity, for each Society, as the midpoint of the lower and upper recorded bounds:
avgpercentsbySociety <- (percentsbySocietyHigh+percentsbySocietyLow)/2;

# Prepare a denominator to normalize the proportion of time spent in each activity, for each Society, to 1:
denom <- t(matrix(rep(colSums(avgpercentsbySociety),length(activities)),nrow=length(Society), ncol=length(activities)));

# Calculate the proportion of time spent in each activity, for each Society, ensuring that the proportions add up to 1:
avgpropsbySociety <- avgpercentsbySociety/denom;
write.csv(avgpropsbySociety,"avgpropbySociety.csv")
# Calculate the division of labour by Society: 
divbySociety = round(colSums(avgpropsbySociety*(divbySociety-3)),2);
write.csv(divbySociety,"divbySociety.csv")

