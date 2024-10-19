library(readxl)
library(dilutBMS2)
library(BMS)
library(fixest)
library(corrplot)

library(multiwayvcov)
library(foreign)
library(multcomp)
library(ggplot2)
library(dplyr)
library(forcats)
library(AER)
library(puniform)
library(Hmisc)
library(DescTools)
library(plm)
library(fwildclusterboot)
library(fixest)
library("RoBMA")
library("data.table")
library(MetaStudies)
library(phack)
library(tidyverse)


Data <- read_excel("C:/Users/USER PC/Documents/Data.xlsx", sheet = "All")
View(Data)



#Instrumented variance and standard errors
#Hike
Hike <- Data[Data$hike == 1, ]
Hike$Inv_Obs = 1/Hike$Obs
Hike$Var = Hike$se*Hike$se
Hike_model <-feols(Var ~ Inv_Obs, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", data = Hike)
summary(Hike_model)
Inv_Obs <- data.frame(Inv_Obs= Hike$Inv_Obs)
Hike$Var_instru <- abs(predict(Hike_model, newdata = Inv_Obs))
Hike$se_instru <- sqrt(Hike$Var_instru)

Hike <- Hike[, c("S|No.", "se_instru", "Var_instru")]

#Cut
Cut <- Data[Data$cut == 1, ]
Cut$Inv_Obs = 1/Cut$Obs
Cut$Var = Cut$se*Cut$se
Cut_model <-feols(Var ~ Inv_Obs, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", data = Cut)
summary(Cut_model)
Inv_Obs <- data.frame(Inv_Obs= Cut$Inv_Obs)
Cut$Var_instru <- abs(predict(Cut_model, newdata = Inv_Obs))
Cut$se_instru <- sqrt(Cut$Var_instru)

Cut <- Cut[, c("S|No.", "se_instru", "Var_instru")]


Var_instru_list <- list(Hike, Cut)
Var_instru_list %>% reduce(full_join)
Var_instru_list <- do.call(rbind.data.frame, Var_instru_list)
Var_instru_list <- Var_instru_list[order(Var_instru_list$`S|No.`), ]
Data$Var_instru <- Var_instru_list$Var_instru
Data$se_instru <- Var_instru_list$se_instru




#######################################
############### FAT-PET  ##############
#######################################


#Upward pass-through

Hike <- Data[Data$hike == 1, ]
Hike <- Hike[, c("study_id", "beta", "se", "Var_instru")]
Hike <- na.omit(Hike)

#The weights
Hike$Inv_SE <- 1/abs(Hike$Var_instru)

FAT_PET<-feols(beta ~ se, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", weights = ~Inv_SE, data = Hike)
summary(FAT_PET)

#Wild Bootstrap intervals
FAT_PET.B1<-boottest(FAT_PET,clustid = "study_id", param = c("se"), B=9999)
FAT_PET.B0<-boottest(FAT_PET,clustid = "study_id", param = c("(Intercept)"), B=9999)
summary(FAT_PET.B1)
summary(FAT_PET.B0)



#Downward pass-through

Cut <- Data[Data$cut == 1, ]
Cut <- Cut[, c("study_id", "beta", "se", "Var_instru")]
Cut <- na.omit(Cut)


#The weights
Cut$Inv_SE <- 1/abs(Cut$Var_instru)

FAT_PET<-feols(beta ~ se, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", weights = ~Inv_SE, data = Cut)
summary(FAT_PET)

#Wild Bootstrap intervals
FAT_PET.B1<-boottest(FAT_PET,clustid = "study_id", param = c("se"), B=9999)
FAT_PET.B0<-boottest(FAT_PET,clustid = "study_id", param = c("(Intercept)"), B=9999)
summary(FAT_PET.B1)
summary(FAT_PET.B0)



#Upward pass-through: Interbank

Hike <- Data[Data$hike == 1 & Data$interbank == 1, ]
Hike <- Hike[, c("study_id", "beta", "se", "Var_instru")]
Hike <- na.omit(Hike)

#The weights
Hike$Inv_SE <- 1/abs(Hike$Var_instru)

FAT_PET<-feols(beta ~ se, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", weights = ~Inv_SE, data = Hike)
summary(FAT_PET)

#Wild Bootstrap intervals
FAT_PET.B1<-boottest(FAT_PET,clustid = "study_id", param = c("se"), B=9999)
FAT_PET.B0<-boottest(FAT_PET,clustid = "study_id", param = c("(Intercept)"), B=9999)
summary(FAT_PET.B1)
summary(FAT_PET.B0)



#Downward pass-through: Interbank

Cut <- Data[Data$cut == 1 & Data$interbank == 1, ]
Cut <- Cut[, c("study_id", "beta", "se", "Var_instru")]
Cut <- na.omit(Cut)

#The weights
Cut$Inv_SE <- 1/abs(Cut$Var_instru)

FAT_PET<-feols(beta ~ se, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", weights = ~Inv_SE, data = Cut)
summary(FAT_PET)

#Wild Bootstrap intervals
FAT_PET.B1<-boottest(FAT_PET,clustid = "study_id", param = c("se"), B=9999)
FAT_PET.B0<-boottest(FAT_PET,clustid = "study_id", param = c("(Intercept)"), B=9999)
summary(FAT_PET.B1)
summary(FAT_PET.B0)



#Upward pass-through: Discount

Hike <- Data[Data$hike == 1 & Data$discount == 1, ]
Hike <- Hike[, c("study_id", "beta", "se", "Var_instru")]
Hike <- na.omit(Hike)

#The weights
Hike$Inv_SE <- 1/abs(Hike$Var_instru)

FAT_PET<-feols(beta ~ se, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", weights = ~Inv_SE, data = Hike)
summary(FAT_PET)

#Wild Bootstrap intervals
FAT_PET.B1<-boottest(FAT_PET,clustid = "study_id", param = c("se"), B=9999)
FAT_PET.B0<-boottest(FAT_PET,clustid = "study_id", param = c("(Intercept)"), B=9999)
summary(FAT_PET.B1)
summary(FAT_PET.B0)



#Downward pass-through: Discount

Cut <- Data[Data$cut == 1 & Data$discount == 1, ]
Cut <- Cut[, c("study_id", "beta", "se", "Var_instru")]
Cut <- na.omit(Cut)

#The weights
Cut$Inv_SE <- 1/abs(Cut$Var_instru)

FAT_PET<-feols(beta ~ se, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", weights = ~Inv_SE, data = Cut)
summary(FAT_PET)

#Wild Bootstrap intervals
FAT_PET.B1<-boottest(FAT_PET,clustid = "study_id", param = c("se"), B=9999)
FAT_PET.B0<-boottest(FAT_PET,clustid = "study_id", param = c("(Intercept)"), B=9999)
summary(FAT_PET.B1)
summary(FAT_PET.B0)




#######################################
########## Selection model  ###########
#######################################
# Codes for Selection Model by Andrews and Kasy (2019)

#Upward pass-through
Hike <- Data[Data$hike == 1, ]
Hike <- Hike[,c("beta", "se")]

ms = metastudies_estimation(X=Hike$beta, sigma=Hike$se, model = "t", cutoffs = c(1.96), symmetric = TRUE)
estimates_plot(ms)
ms$est_tab


#Downward pass-through
Cut <- Data[Data$cut == 1, ]
Cut <- Cut[,c("beta", "se")]

ms = metastudies_estimation(X=Cut$beta, sigma=Cut$se, model = "t", cutoffs = c(1.96), symmetric = TRUE)
estimates_plot(ms)
ms$est_tab



#Upward pass-through: interbank
Hike <- Data[Data$hike == 1 & Data$interbank == 1, ]
Hike <- Hike[,c("beta", "se")]

ms = metastudies_estimation(X=Hike$beta, sigma=Hike$se, model = "t", cutoffs = c(1.96), symmetric = TRUE)
estimates_plot(ms)
ms$est_tab


#Downward pass-through: interbank
Cut <- Data[Data$cut == 1 & Data$interbank == 1, ]
Cut <- Cut[,c("beta", "se")]

ms = metastudies_estimation(X=Cut$beta, sigma=Cut$se, model = "t", cutoffs = c(1.96), symmetric = TRUE)
estimates_plot(ms)
ms$est_tab



#Upward pass-through: discount
Hike <- Data[Data$hike == 1 & Data$discount == 1, ]
Hike <- Hike[,c("beta", "se")]

ms = metastudies_estimation(X=Hike$beta, sigma=Hike$se, model = "t", cutoffs = c(1.96), symmetric = TRUE)
estimates_plot(ms)
ms$est_tab


#Downward pass-through: discount
Cut <- Data[Data$cut == 1 & Data$discount == 1, ]
Cut <- Cut[,c("beta", "se")]

ms = metastudies_estimation(X=Cut$beta, sigma=Cut$se, model = "t", cutoffs = c(1.96), symmetric = TRUE)
estimates_plot(ms)
ms$est_tab




#FOR APPENDIX:
#Fixed effects model:  

Data$beta_wgt = Data$beta/Data$se_instru
Data$con_wgt = 1/Data$se_instru
Data$se_wgt = Data$se /Data$se_instru


#Upward pass-through
Hike <- Data[Data$hike == 1, ]
Hike <- Hike[, c("study_id", "beta_wgt", "con_wgt", "se_wgt")]
Hike <- na.omit(Hike)

FE <-feols(beta_wgt ~ con_wgt + se_wgt, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", fixef = "study_id", data = Hike)
summary(FE)




#Downward pass-through
Cut <- Data[Data$cut == 1, ]
Cut <- Cut[, c("study_id", "beta_wgt", "con_wgt", "se_wgt")]
Cut <- na.omit(Cut)

FE <-feols(beta_wgt ~ con_wgt + se_wgt, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", fixef = "study_id", data = Cut)
summary(FE)




#Upward pass-through: interbank
Hike <- Data[Data$hike == 1 & Data$interbank == 1, ]
Hike <- Hike[, c("study_id", "beta_wgt", "con_wgt", "se_wgt")]
Hike <- na.omit(Hike)

FE <-feols(beta_wgt ~ con_wgt + se_wgt, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", fixef = "study_id", data = Hike)
summary(FE)




#Downward pass-through: interbank  
Cut <- Data[Data$cut == 1 & Data$interbank == 1, ]
Cut <- Cut[, c("study_id", "beta_wgt", "con_wgt", "se_wgt")]
Cut <- na.omit(Cut)

FE <-feols(beta_wgt ~ con_wgt + se_wgt, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", fixef = "study_id", data = Cut)
summary(FE)




#Upward pass-through: discount
Hike <- Data[Data$hike == 1 & Data$discount == 1, ]
Hike <- Hike[, c("study_id", "beta_wgt", "con_wgt", "se_wgt")]
Hike <- na.omit(Hike)

FE <-feols(beta_wgt ~ con_wgt + se_wgt, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", fixef = "study_id", data = Hike)
summary(FE)




#Downward pass-through: discount  
Cut <- Data[Data$cut == 1 & Data$discount == 1, ]
Cut <- Cut[, c("study_id", "beta_wgt", "con_wgt", "se_wgt")]
Cut <- na.omit(Cut)

FE <-feols(beta_wgt ~ con_wgt + se_wgt, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", fixef = "study_id", data = Cut)
summary(FE)




#Random effects model (GLS):  

Data$beta_wgt = Data$beta/Data$se_instru
Data$con_wgt = 1/Data$se_instru
Data$se_wgt = Data$se /Data$se_instru


#Upward pass-through
Hike <- Data[Data$hike == 1, ]
Hike <- Hike[, c("study_id", "beta_wgt", "con_wgt", "se_wgt")]
Hike <- na.omit(Hike)

RE <- plm(beta_wgt ~ 0 + con_wgt + se_wgt, index = "study_id", data = Hike, model ="random", effect="time")
summary(RE)
coeftest(RE, vcov=vcovHC(RE,type="HC0",cluster="group")) 



#Downward pass-through
Cut <- Data[Data$cut == 1, ]
Cut <- Cut[, c("study_id", "beta_wgt", "con_wgt", "se_wgt")]
Cut <- na.omit(Cut)

RE <- plm(beta_wgt ~ 0 + con_wgt + se_wgt, index = "study_id", data = Cut, model ="random", effect="time")
summary(RE)
coeftest(RE, vcov=vcovHC(RE,type="HC0",cluster="group")) 



#Upward pass-through: interbank
Hike <- Data[Data$hike == 1 & Data$interbank == 1, ]
Hike <- Hike[, c("study_id", "beta_wgt", "con_wgt", "se_wgt")]
Hike <- na.omit(Hike)

RE <- plm(beta_wgt ~ 0 + con_wgt + se_wgt, index = "study_id", data = Hike, model ="random", effect="time")
summary(RE)
coeftest(RE, vcov=vcovHC(RE,type="HC0",cluster="group")) 



#Downward pass-through: interbank
Cut <- Data[Data$cut == 1 & Data$interbank == 1, ]
Cut <- Cut[, c("study_id", "beta_wgt", "con_wgt", "se_wgt")]
Cut <- na.omit(Cut)

RE <- plm(beta_wgt ~ 0 + con_wgt + se_wgt, index = "study_id", data = Cut, model ="random", effect="time")
summary(RE)
coeftest(RE, vcov=vcovHC(RE,type="HC0",cluster="group")) 




#Upward pass-through: discount
Hike <- Data[Data$hike == 1 & Data$discount == 1, ]
Hike <- Hike[, c("study_id", "beta_wgt", "con_wgt", "se_wgt")]
Hike <- na.omit(Hike)

RE <- plm(beta_wgt ~ 0 + con_wgt + se_wgt, index = "study_id", data = Hike, model ="random", effect="time")
summary(RE)
coeftest(RE, vcov=vcovHC(RE,type="HC0",cluster="group")) 



#Downward pass-through: discount
Cut <- Data[Data$cut == 1 & Data$discount == 1, ]
Cut <- Cut[, c("study_id", "beta_wgt", "con_wgt", "se_wgt")]
Cut <- na.omit(Cut)

RE <- plm(beta_wgt ~ 0 + con_wgt + se_wgt, index = "study_id", data = Cut, model ="random", effect="time")
summary(RE)
coeftest(RE, vcov=vcovHC(RE,type="HC0",cluster="group")) 





#WAAP method here: 
#Upward pass-through

Hike <- Data[Data$hike == 1, ]
Hike <- Hike[, c("study_id", "beta", "se", "Var_instru")]
Hike$se_instru <- sqrt(abs(Hike$Var_instru))
Hike <- na.omit(Hike)

#The weights
Hike$Inv_SE <- 1/abs(Hike$Var_instru)

WAAP_true_ef<-sum(Hike$beta)/nrow(Hike)
WAAP_data<-Hike[WAAP_true_ef/2.8 > Hike$se_instru,]
n_distinct(WAAP_data$study_id)

WAAP<-feols(beta ~ se, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", weights = ~Inv_SE, data = WAAP_data)
summary(WAAP)



#Downward pass-through

Cut <- Data[Data$cut == 1, ]
Cut <- Cut[, c("study_id", "beta", "se", "Var_instru")]
Cut$se_instru <- sqrt(abs(Cut$Var_instru))
Cut <- na.omit(Cut)

#The weights
Cut$Inv_SE <- 1/abs(Cut$Var_instru)

WAAP_true_ef<-sum(Cut$beta)/nrow(Cut)
WAAP_data<-Cut[WAAP_true_ef/2.8 > Cut$se_instru,]
n_distinct(WAAP_data$study_id)

WAAP<-feols(beta ~ se, ssc=ssc(adj = FALSE, cluster.adj=FALSE),  cluster = "study_id", weights = ~Inv_SE, data = WAAP_data)
summary(WAAP)




#Upward pass-through: interbank

Hike <- Data[Data$hike == 1 & Data$interbank == 1, ]
Hike <- Hike[, c("study_id", "beta", "se", "Var_instru")]
Hike$se_instru <- sqrt(abs(Hike$Var_instru))
Hike <- na.omit(Hike)

#The weights
Hike$Inv_SE <- 1/abs(Hike$Var_instru)

WAAP_true_ef<-sum(Hike$beta)/nrow(Hike)
WAAP_data<-Hike[WAAP_true_ef/2.8 > Hike$se_instru,]
n_distinct(WAAP_data$study_id)

WAAP<-feols(beta ~ se, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", weights = ~Inv_SE, data = WAAP_data)
summary(WAAP)



#Downward pass-through:interbank

Cut <- Data[Data$cut == 1 & Data$interbank == 1, ]
Cut <- Cut[, c("study_id", "beta", "se", "Var_instru")]
Cut$se_instru <- sqrt(abs(Cut$Var_instru))
Cut <- na.omit(Cut)

#The weights
Cut$Inv_SE <- 1/abs(Cut$Var_instru)

WAAP_true_ef<-sum(Cut$beta)/nrow(Cut)
WAAP_data<-Cut[WAAP_true_ef/2.8 > Cut$se_instru,]
n_distinct(WAAP_data$study_id)

WAAP<-feols(beta ~ se, ssc=ssc(adj = FALSE, cluster.adj=FALSE),  cluster = "study_id", weights = ~Inv_SE, data = WAAP_data)
summary(WAAP)




#Upward pass-through:discount

Hike <- Data[Data$hike == 1 & Data$discount == 1, ]
Hike <- Hike[, c("study_id", "beta", "se", "Var_instru")]
Hike$se_instru <- sqrt(abs(Hike$Var_instru))
Hike <- na.omit(Hike)

#The weights
Hike$Inv_SE <- 1/abs(Hike$Var_instru)

WAAP_true_ef<-sum(Hike$beta)/nrow(Hike)
WAAP_data<-Hike[WAAP_true_ef/2.8 > Hike$se_instru,]
n_distinct(WAAP_data$study_id)

WAAP<-feols(beta ~ se, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", weights = ~Inv_SE, data = WAAP_data)
summary(WAAP)



#Downward pass-through:discount

Cut <- Data[Data$cut == 1 & Data$discount == 1, ]
Cut <- Cut[, c("study_id", "beta", "se", "Var_instru")]
Cut$se_instru <- sqrt(abs(Cut$Var_instru))
Cut <- na.omit(Cut)

#The weights
Cut$Inv_SE <- 1/abs(Cut$Var_instru)

WAAP_true_ef<-sum(Cut$beta)/nrow(Cut)
WAAP_data<-Cut[WAAP_true_ef/2.8 > Cut$se_instru,]
n_distinct(WAAP_data$study_id)

WAAP<-feols(beta ~ se, ssc=ssc(adj = FALSE, cluster.adj=FALSE),  cluster = "study_id", weights = ~Inv_SE, data = WAAP_data)
summary(WAAP)





#Top 10%:

#Upward pass-through:
Data$prec <- 1/Data$se_instru
Hike <- Data[Data$hike == 1, ]
Hike$Inv_SE <- 1/abs(Hike$Var_instru)
Hike_sort <- Hike[order(-Hike$prec),]
First_10 <- 0.1 * nrow(Hike_sort)
Hike_top_10 <- head(Hike_sort,First_10)
n_distinct(Hike_top_10$study_id)

Top_10 <-feols(beta ~ se, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", weights = ~Inv_SE, data = Hike_top_10)
summary(Top_10)



#Downward pass-through:
Data$prec <- 1/Data$se_instru
Cut <- Data[Data$cut == 1, ]
Cut$Inv_SE <- 1/abs(Cut$Var_instru)
Cut_sort <- Cut[order(-Cut$prec),]
First_10 <- 0.1 * nrow(Cut_sort)
Cut_top_10 <- head(Cut_sort,First_10)
n_distinct(Cut_top_10$study_id)

Top_10 <-feols(beta ~ se, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", weights = ~Inv_SE, data = Cut_top_10)
summary(Top_10)




#Upward pass-through: interbank
Data$prec <- 1/Data$se_instru
Hike <- Data[Data$hike == 1 & Data$interbank == 1, ]
Hike$Inv_SE <- 1/abs(Hike$Var_instru)
Hike_sort <- Hike[order(-Hike$prec),]
First_10 <- 0.1 * nrow(Hike_sort)
Hike_top_10 <- head(Hike_sort,First_10)
n_distinct(Hike_top_10$study_id)

Top_10 <-feols(beta ~ se, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", weights = ~Inv_SE, data = Hike_top_10)
summary(Top_10)



#Downward pass-through: interbank
Data$prec <- 1/Data$se_instru
Cut <- Data[Data$cut == 1 & Data$interbank == 1, ]
Cut$Inv_SE <- 1/abs(Cut$Var_instru)
Cut_sort <- Cut[order(-Cut$prec),]
First_10 <- 0.1 * nrow(Cut_sort)
Cut_top_10 <- head(Cut_sort,First_10)
n_distinct(Cut_top_10$study_id)

Top_10 <-feols(beta ~ se, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", weights = ~Inv_SE, data = Cut_top_10)
summary(Top_10)




#Upward pass-through: discount
Data$prec <- 1/Data$se_instru
Hike <- Data[Data$hike == 1 & Data$discount == 1, ]
Hike$Inv_SE <- 1/abs(Hike$Var_instru)
Hike_sort <- Hike[order(-Hike$prec),]
First_10 <- 0.1 * nrow(Hike_sort)
Hike_top_10 <- head(Hike_sort,First_10)
n_distinct(Hike_top_10$study_id)

Top_10 <-feols(beta ~ se, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", weights = ~Inv_SE, data = Hike_top_10)
summary(Top_10)



#Downward pass-through: discount
Data$prec <- 1/Data$se_instru
Cut <- Data[Data$cut == 1 & Data$discount == 1, ]
Cut$Inv_SE <- 1/abs(Cut$Var_instru)
Cut_sort <- Cut[order(-Cut$prec),]
First_10 <- 0.1 * nrow(Cut_sort)
Cut_top_10 <- head(Cut_sort,First_10)
n_distinct(Cut_top_10$study_id)

Top_10 <-feols(beta ~ se, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", weights = ~Inv_SE, data = Cut_top_10)
summary(Top_10)




  
##########
#Set seed#
##########

set.seed(12345)




#######################################
########### P-hacking tests ###########
#######################################



#Create backup
Backup <- Data


#Upward pass-through
Data <- Backup[Backup$hike == 1, ]

# Derounding:
Data$added <- 0.005
Data$Ones <- 1

es_upper <- Data$beta + Data$added
es_lower <- Data$beta - Data$added

se_upper <- Data$se + Data$added
se_lower <- Data$se - Data$added

beta.deround = runif(Data$Ones, es_lower, es_upper)
se.deround = runif(Data$Ones, se_lower, se_upper)

t.deround = beta.deround / se.deround
p.deround = 2*(1-pnorm(t.deround))
study_id <- Data$study_id

Data$p.deround <- p.deround

phack_test_cox_shi(Data$p.deround, Data$study_id, p_min=0.0000, p_max=0.15, J=10, K=2, use_bound=TRUE)
phack_test_cox_shi(Data$p.deround, Data$study_id, p_min=0.0000, p_max=0.15, J=10)

Obs <- Data[Data$p.deround<=0.15,]
count(Obs)



#Downward pass-through
Data <- Backup[Backup$cut == 1, ]

Data$added <- 0.005
Data$Ones <- 1

es_upper <- Data$beta + Data$added
es_lower <- Data$beta - Data$added

se_upper <- Data$se + Data$added
se_lower <- Data$se - Data$added

beta.deround = runif(Data$Ones, es_lower, es_upper)
se.deround = runif(Data$Ones, se_lower, se_upper)

t.deround = beta.deround / se.deround
p.deround = 2*(1-pnorm(t.deround))
study_id <- Data$study_id

Data$p.deround <- p.deround

phack_test_cox_shi(Data$p.deround, Data$study_id, p_min=0.0000, p_max=0.15, J=10, K=2, use_bound=TRUE)
phack_test_cox_shi(Data$p.deround, Data$study_id, p_min=0.0000, p_max=0.15, J=10)

Obs <- Data[Data$p.deround<=0.15,]
count(Obs)





#Upward pass-through: interbank
Data <- Backup[Backup$hike == 1 & Backup$interbank == 1, ]

# Derounding:
Data$added <- 0.005
Data$Ones <- 1

es_upper <- Data$beta + Data$added
es_lower <- Data$beta - Data$added

se_upper <- Data$se + Data$added
se_lower <- Data$se - Data$added

beta.deround = runif(Data$Ones, es_lower, es_upper)
se.deround = runif(Data$Ones, se_lower, se_upper)

t.deround = beta.deround / se.deround
p.deround = 2*(1-pnorm(t.deround))
study_id <- Data$study_id

Data$p.deround <- p.deround

phack_test_cox_shi(Data$p.deround, Data$study_id, p_min=0.0000, p_max=0.15, J=5, K=2, use_bound=TRUE)
phack_test_cox_shi(Data$p.deround, Data$study_id, p_min=0.0000, p_max=0.15, J=5)

Obs <- Data[Data$p.deround<=0.15,]
count(Obs)



#Downward pass-through: interbank
Data <- Backup[Backup$cut == 1 & Backup$interbank == 1, ]

Data$added <- 0.005
Data$Ones <- 1

es_upper <- Data$beta + Data$added
es_lower <- Data$beta - Data$added

se_upper <- Data$se + Data$added
se_lower <- Data$se - Data$added

beta.deround = runif(Data$Ones, es_lower, es_upper)
se.deround = runif(Data$Ones, se_lower, se_upper)

t.deround = beta.deround / se.deround
p.deround = 2*(1-pnorm(t.deround))
study_id <- Data$study_id

Data$p.deround <- p.deround

phack_test_cox_shi(Data$p.deround, Data$study_id, p_min=0.0000, p_max=0.15, J=5, K=2, use_bound=TRUE)
phack_test_cox_shi(Data$p.deround, Data$study_id, p_min=0.0000, p_max=0.15, J=5)

Obs <- Data[Data$p.deround<=0.15,]
count(Obs)





#Upward pass-through: discount
Data <- Backup[Backup$hike == 1 & Backup$discount == 1, ]

# Derounding:
Data$added <- 0.005
Data$Ones <- 1

es_upper <- Data$beta + Data$added
es_lower <- Data$beta - Data$added

se_upper <- Data$se + Data$added
se_lower <- Data$se - Data$added

beta.deround = runif(Data$Ones, es_lower, es_upper)
se.deround = runif(Data$Ones, se_lower, se_upper)

t.deround = beta.deround / se.deround
p.deround = 2*(1-pnorm(t.deround))
study_id <- Data$study_id

Data$p.deround <- p.deround

phack_test_cox_shi(Data$p.deround, Data$study_id, p_min=0.0000, p_max=0.15, J=5, K=2, use_bound=TRUE)
phack_test_cox_shi(Data$p.deround, Data$study_id, p_min=0.0000, p_max=0.15, J=5)

Obs <- Data[Data$p.deround<=0.15,]
count(Obs)



#Downward pass-through: discount
Data <- Backup[Backup$cut == 1 & Backup$discount == 1, ]

Data$added <- 0.005
Data$Ones <- 1

es_upper <- Data$beta + Data$added
es_lower <- Data$beta - Data$added

se_upper <- Data$se + Data$added
se_lower <- Data$se - Data$added

beta.deround = runif(Data$Ones, es_lower, es_upper)
se.deround = runif(Data$Ones, se_lower, se_upper)

t.deround = beta.deround / se.deround
p.deround = 2*(1-pnorm(t.deround))
study_id <- Data$study_id

Data$p.deround <- p.deround

phack_test_cox_shi(Data$p.deround, Data$study_id, p_min=0.0000, p_max=0.15, J=5, K=2, use_bound=TRUE)
phack_test_cox_shi(Data$p.deround, Data$study_id, p_min=0.0000, p_max=0.15, J=5)

Obs <- Data[Data$p.deround<=0.15,]
count(Obs)


#Return backed-up data
Data <- Backup




## Bayesian model averaging

#Correlation matrix
bma_COR <- read_excel("C:/Users/USER PC/Documents/Data.xlsx", sheet = "CORR")
bma_COR <- na.omit(bma_COR)
BMA_COR<-cor(bma_COR)
corrplot(BMA_COR, type="lower", number.cex = 0.15, sig.level = 0.01, tl.col="black", addCoef.col = "black", method="color", tl.srt=45)


#set seed for replication
set.seed(12345)


IRPT <- Data[, c("beta", "se", "interbank", "hike", "T_span", "PreCrises", "PostCrises", "Gov", "IMF", "ECB", "Quarterly", "Monthly", "Weekly", "Time_Series", "Lag_length", "Long_run", "Short_run", "Infl_con", "OLS", "Johansen", "P_year", "A_Cite", "Peer", "IMPF", "Central_banker", "Firms", "Households", "Mortgage", "S_loan", "L_loan", "Imp_Open", "FDI_Open", "Dev", "Cap", "Infl", "GDP_gr", "Stock_turn", "CBI", "M_policy", "Ex_regime", "Eurozone")]


#Construct Weights:
Data$Weight_est <- 1/Data$est_no

#Weight the variables:
IRPT <- IRPT*Data$Weight_est

#running the model
bma_irpt = bms(IRPT, burn=1000000, iter=2000000, mprior="dilut", g="UIP")




#Before running this piece of code, save the work environment and restart R session:

library("BMS")

#Baseline BMA model inclusion graph
image(bma_irpt,order.by.pip = TRUE)


#Baseline BMA model size and convergence
plot(bma_irpt, include.legend=TRUE)


#Baseline BMA model fitted intercept
coef(bma_irpt, include.constant = TRUE)


#Baseline BMA model summary
summary(bma_irpt)


#Coefficient distributions: plotted sequentially
density(bma_irpt, reg = 30) #Openness to import trade
density(bma_irpt, reg = 31) #Openness to FDI inflows
density(bma_irpt, reg = 32) #Economic development status
density(bma_irpt, reg = 34) #Inflationary environment
density(bma_irpt, reg = 3)  #Upward pass-through


#Implied estimates:


#Overall:Short-run
#Hike
new.dat <- data.frame(se=0, interbank=0.7763, hike=1, T_span=12.0446, PreCrises=0, PostCrises=1, Gov=0.8249, IMF=0.0796, ECB=0.0835, Quarterly=1, Monthly=0, Weekly=0, Time_Series=0.9446, Lag_length=0.6275, Long_run=0, Short_run=1, Infl_con=0.0435, OLS=0, Johansen=1, P_year=19.1277, A_Cite=22.65, Peer=0.8332, IMPF=0.0756, Central_banker=0.3156, Firms=0.2697, Households=0.0877, Mortgage=0.0572, S_loan=0.3047, L_loan=0.0685, Imp_Open=0.2473, FDI_Open=0.0272, Dev=0.5203, Cap=0.7019, Infl=0.0538, GDP_gr=0.0396, Stock_turn=0.9999, CBI=0.5569, M_policy=0.4508, Ex_regime=0.7669, Eurozone=0.2124)
predict(bma_irpt, newdata = new.dat, interval = 'confidence')

#OLS confidence intervals
OLSCONF <- feols(beta ~ se + interbank + hike + T_span + PreCrises + PostCrises + Gov + IMF + ECB + Quarterly + Monthly + Weekly + Time_Series + Lag_length + Long_run + Short_run + Infl_con + OLS + Johansen + P_year + A_Cite + Peer + IMPF + Central_banker + Firms + Households + Mortgage + S_loan + L_loan + Imp_Open + FDI_Open + Dev + Cap + Infl + GDP_gr +  Stock_turn + CBI + M_policy + Ex_regime + Eurozone,  ssc=ssc(adj = FALSE, cluster.adj=FALSE), weights = ~Weight_est, cluster = "study_id", panel.id = "study_id",  data = Data)
predict(OLSCONF, newdata = new.dat, interval = 'confidence')

#Cut
new.dat <- data.frame(se=0, interbank=0.7763, hike=0, T_span=12.0446, PreCrises=0, PostCrises=1, Gov=0.8249, IMF=0.0796, ECB=0.0835, Quarterly=1, Monthly=0, Weekly=0, Time_Series=0.9446, Lag_length=0.6275, Long_run=0, Short_run=1, Infl_con=0.0435, OLS=0, Johansen=1, P_year=19.1277, A_Cite=22.65, Peer=0.8332, IMPF=0.0756, Central_banker=0.3156, Firms=0.2697, Households=0.0877, Mortgage=0.0572, S_loan=0.3047, L_loan=0.0685, Imp_Open=0.2473, FDI_Open=0.0272, Dev=0.5203, Cap=0.7019, Infl=0.0538, GDP_gr=0.0396, Stock_turn=0.9999, CBI=0.5569, M_policy=0.4508, Ex_regime=0.7669, Eurozone=0.2124)
predict(bma_irpt, newdata = new.dat, interval = 'confidence')

#OLS confidence intervals
OLSCONF <- feols(beta ~ se + interbank + hike + T_span + PreCrises + PostCrises + Gov + IMF + ECB + Quarterly + Monthly + Weekly + Time_Series + Lag_length + Long_run + Short_run + Infl_con + OLS + Johansen + P_year + A_Cite + Peer + IMPF + Central_banker + Firms + Households + Mortgage + S_loan + L_loan + Imp_Open + FDI_Open + Dev + Cap + Infl + GDP_gr +  Stock_turn + CBI + M_policy + Ex_regime + Eurozone,  ssc=ssc(adj = FALSE, cluster.adj=FALSE), weights = ~Weight_est, cluster = "study_id", panel.id = "study_id",  data = Data)
predict(OLSCONF, newdata = new.dat, interval = 'confidence')



#Overall:Long-run
#Hike
new.dat <- data.frame(se=0, interbank=0.7763, hike=1, T_span=12.0446, PreCrises=0, PostCrises=1, Gov=0.8249, IMF=0.0796, ECB=0.0835, Quarterly=1, Monthly=0, Weekly=0, Time_Series=0.9446, Lag_length=0.6275, Long_run=1, Short_run=0, Infl_con=0.0435, OLS=0, Johansen=1, P_year=19.1277, A_Cite=22.65, Peer=0.8332, IMPF=0.0756, Central_banker=0.3156, Firms=0.2697, Households=0.0877, Mortgage=0.0572, S_loan=0.3047, L_loan=0.0685, Imp_Open=0.2473, FDI_Open=0.0272, Dev=0.5203, Cap=0.7019, Infl=0.0538, GDP_gr=0.0396, Stock_turn=0.9999, CBI=0.5569, M_policy=0.4508, Ex_regime=0.7669, Eurozone=0.2124)
predict(bma_irpt, newdata = new.dat, interval = 'confidence')

#OLS confidence intervals
OLSCONF <- feols(beta ~ se + interbank + hike + T_span + PreCrises + PostCrises + Gov + IMF + ECB + Quarterly + Monthly + Weekly + Time_Series + Lag_length + Long_run + Short_run + Infl_con + OLS + Johansen + P_year + A_Cite + Peer + IMPF + Central_banker + Firms + Households + Mortgage + S_loan + L_loan + Imp_Open + FDI_Open + Dev + Cap + Infl + GDP_gr +  Stock_turn + CBI + M_policy + Ex_regime + Eurozone,  ssc=ssc(adj = FALSE, cluster.adj=FALSE), weights = ~Weight_est, cluster = "study_id", panel.id = "study_id",  data = Data)
predict(OLSCONF, newdata = new.dat, interval = 'confidence')

#Cut
new.dat <- data.frame(se=0, interbank=0.7763, hike=0, T_span=12.0446, PreCrises=0, PostCrises=1, Gov=0.8249, IMF=0.0796, ECB=0.0835, Quarterly=1, Monthly=0, Weekly=0, Time_Series=0.9446, Lag_length=0.6275, Long_run=1, Short_run=0, Infl_con=0.0435, OLS=0, Johansen=1, P_year=19.1277, A_Cite=22.65, Peer=0.8332, IMPF=0.0756, Central_banker=0.3156, Firms=0.2697, Households=0.0877, Mortgage=0.0572, S_loan=0.3047, L_loan=0.0685, Imp_Open=0.2473, FDI_Open=0.0272, Dev=0.5203, Cap=0.7019, Infl=0.0538, GDP_gr=0.0396, Stock_turn=0.9999, CBI=0.5569, M_policy=0.4508, Ex_regime=0.7669, Eurozone=0.2124)
predict(bma_irpt, newdata = new.dat, interval = 'confidence')

#OLS confidence intervals
OLSCONF <- feols(beta ~ se + interbank + hike + T_span + PreCrises + PostCrises + Gov + IMF + ECB + Quarterly + Monthly + Weekly + Time_Series + Lag_length + Long_run + Short_run + Infl_con + OLS + Johansen + P_year + A_Cite + Peer + IMPF + Central_banker + Firms + Households + Mortgage + S_loan + L_loan + Imp_Open + FDI_Open + Dev + Cap + Infl + GDP_gr +  Stock_turn + CBI + M_policy + Ex_regime + Eurozone,  ssc=ssc(adj = FALSE, cluster.adj=FALSE), weights = ~Weight_est, cluster = "study_id", panel.id = "study_id",  data = Data)
predict(OLSCONF, newdata = new.dat, interval = 'confidence')





#Developed:Short-run
#Hike
new.dat <- data.frame(se=0, interbank=0.7763, hike=1, T_span=12.0446, PreCrises=0, PostCrises=1, Gov=0.8249, IMF=0.0796, ECB=0.0835, Quarterly=1, Monthly=0, Weekly=0, Time_Series=0.9446, Lag_length=0.6275, Long_run=0, Short_run=1, Infl_con=0.0435, OLS=0, Johansen=1, P_year=19.1277, A_Cite=22.65, Peer=0.8332, IMPF=0.0756, Central_banker=0.3156, Firms=0.2697, Households=0.0877, Mortgage=0.0572, S_loan=0.3047, L_loan=0.0685, Imp_Open=0.2473, FDI_Open=0.0272, Dev=1, Cap=0.7019, Infl=0.0538, GDP_gr=0.0396, Stock_turn=0.9999, CBI=0.5569, M_policy=0.4508, Ex_regime=0.7669, Eurozone=0.2124)
predict(bma_irpt, newdata = new.dat, interval = 'confidence')

#OLS confidence intervals
OLSCONF <- feols(beta ~ se + interbank + hike + T_span + PreCrises + PostCrises + Gov + IMF + ECB + Quarterly + Monthly + Weekly + Time_Series + Lag_length + Long_run + Short_run + Infl_con + OLS + Johansen + P_year + A_Cite + Peer + IMPF + Central_banker + Firms + Households + Mortgage + S_loan + L_loan + Imp_Open + FDI_Open + Dev + Cap + Infl + GDP_gr +  Stock_turn + CBI + M_policy + Ex_regime + Eurozone,  ssc=ssc(adj = FALSE, cluster.adj=FALSE), weights = ~Weight_est, cluster = "study_id", panel.id = "study_id",  data = Data)
predict(OLSCONF, newdata = new.dat, interval = 'confidence')

#Cut
new.dat <- data.frame(se=0, interbank=0.7763, hike=0, T_span=12.0446, PreCrises=0, PostCrises=1, Gov=0.8249, IMF=0.0796, ECB=0.0835, Quarterly=1, Monthly=0, Weekly=0, Time_Series=0.9446, Lag_length=0.6275, Long_run=0, Short_run=1, Infl_con=0.0435, OLS=0, Johansen=1, P_year=19.1277, A_Cite=22.65, Peer=0.8332, IMPF=0.0756, Central_banker=0.3156, Firms=0.2697, Households=0.0877, Mortgage=0.0572, S_loan=0.3047, L_loan=0.0685, Imp_Open=0.2473, FDI_Open=0.0272, Dev=1, Cap=0.7019, Infl=0.0538, GDP_gr=0.0396, Stock_turn=0.9999, CBI=0.5569, M_policy=0.4508, Ex_regime=0.7669, Eurozone=0.2124)
predict(bma_irpt, newdata = new.dat, interval = 'confidence')

#OLS confidence intervals
OLSCONF <- feols(beta ~ se + interbank + hike + T_span + PreCrises + PostCrises + Gov + IMF + ECB + Quarterly + Monthly + Weekly + Time_Series + Lag_length + Long_run + Short_run + Infl_con + OLS + Johansen + P_year + A_Cite + Peer + IMPF + Central_banker + Firms + Households + Mortgage + S_loan + L_loan + Imp_Open + FDI_Open + Dev + Cap + Infl + GDP_gr +  Stock_turn + CBI + M_policy + Ex_regime + Eurozone,  ssc=ssc(adj = FALSE, cluster.adj=FALSE), weights = ~Weight_est, cluster = "study_id", panel.id = "study_id",  data = Data)
predict(OLSCONF, newdata = new.dat, interval = 'confidence')



#Developed:Long-run
#Hike
new.dat <- data.frame(se=0, interbank=0.7763, hike=1, T_span=12.0446, PreCrises=0, PostCrises=1, Gov=0.8249, IMF=0.0796, ECB=0.0835, Quarterly=1, Monthly=0, Weekly=0, Time_Series=0.9446, Lag_length=0.6275, Long_run=1, Short_run=0, Infl_con=0.0435, OLS=0, Johansen=1, P_year=19.1277, A_Cite=22.65, Peer=0.8332, IMPF=0.0756, Central_banker=0.3156, Firms=0.2697, Households=0.0877, Mortgage=0.0572, S_loan=0.3047, L_loan=0.0685, Imp_Open=0.2473, FDI_Open=0.0272, Dev=1, Cap=0.7019, Infl=0.0538, GDP_gr=0.0396, Stock_turn=0.9999, CBI=0.5569, M_policy=0.4508, Ex_regime=0.7669, Eurozone=0.2124)
predict(bma_irpt, newdata = new.dat, interval = 'confidence')

#OLS confidence intervals
OLSCONF <- feols(beta ~ se + interbank + hike + T_span + PreCrises + PostCrises + Gov + IMF + ECB + Quarterly + Monthly + Weekly + Time_Series + Lag_length + Long_run + Short_run + Infl_con + OLS + Johansen + P_year + A_Cite + Peer + IMPF + Central_banker + Firms + Households + Mortgage + S_loan + L_loan + Imp_Open + FDI_Open + Dev + Cap + Infl + GDP_gr +  Stock_turn + CBI + M_policy + Ex_regime + Eurozone,  ssc=ssc(adj = FALSE, cluster.adj=FALSE), weights = ~Weight_est, cluster = "study_id", panel.id = "study_id",  data = Data)
predict(OLSCONF, newdata = new.dat, interval = 'confidence')

#Cut
new.dat <- data.frame(se=0, interbank=0.7763, hike=0, T_span=12.0446, PreCrises=0, PostCrises=1, Gov=0.8249, IMF=0.0796, ECB=0.0835, Quarterly=1, Monthly=0, Weekly=0, Time_Series=0.9446, Lag_length=0.6275, Long_run=1, Short_run=0, Infl_con=0.0435, OLS=0, Johansen=1, P_year=19.1277, A_Cite=22.65, Peer=0.8332, IMPF=0.0756, Central_banker=0.3156, Firms=0.2697, Households=0.0877, Mortgage=0.0572, S_loan=0.3047, L_loan=0.0685, Imp_Open=0.2473, FDI_Open=0.0272, Dev=1, Cap=0.7019, Infl=0.0538, GDP_gr=0.0396, Stock_turn=0.9999, CBI=0.5569, M_policy=0.4508, Ex_regime=0.7669, Eurozone=0.2124)
predict(bma_irpt, newdata = new.dat, interval = 'confidence')

#OLS confidence intervals
OLSCONF <- feols(beta ~ se + interbank + hike + T_span + PreCrises + PostCrises + Gov + IMF + ECB + Quarterly + Monthly + Weekly + Time_Series + Lag_length + Long_run + Short_run + Infl_con + OLS + Johansen + P_year + A_Cite + Peer + IMPF + Central_banker + Firms + Households + Mortgage + S_loan + L_loan + Imp_Open + FDI_Open + Dev + Cap + Infl + GDP_gr +  Stock_turn + CBI + M_policy + Ex_regime + Eurozone,  ssc=ssc(adj = FALSE, cluster.adj=FALSE), weights = ~Weight_est, cluster = "study_id", panel.id = "study_id",  data = Data)
predict(OLSCONF, newdata = new.dat, interval = 'confidence')





#Developing:Short-run
#Hike
new.dat <- data.frame(se=0, interbank=0.7763, hike=1, T_span=12.0446, PreCrises=0, PostCrises=1, Gov=0.8249, IMF=0.0796, ECB=0.0835, Quarterly=1, Monthly=0, Weekly=0, Time_Series=0.9446, Lag_length=0.6275, Long_run=0, Short_run=1, Infl_con=0.0435, OLS=0, Johansen=1, P_year=19.1277, A_Cite=22.65, Peer=0.8332, IMPF=0.0756, Central_banker=0.3156, Firms=0.2697, Households=0.0877, Mortgage=0.0572, S_loan=0.3047, L_loan=0.0685, Imp_Open=0.2473, FDI_Open=0.0272, Dev=0, Cap=0.7019, Infl=0.0538, GDP_gr=0.0396, Stock_turn=0.9999, CBI=0.5569, M_policy=0.4508, Ex_regime=0.7669, Eurozone=0.2124)
predict(bma_irpt, newdata = new.dat, interval = 'confidence')

#OLS confidence intervals
OLSCONF <- feols(beta ~ se + interbank + hike + T_span + PreCrises + PostCrises + Gov + IMF + ECB + Quarterly + Monthly + Weekly + Time_Series + Lag_length + Long_run + Short_run + Infl_con + OLS + Johansen + P_year + A_Cite + Peer + IMPF + Central_banker + Firms + Households + Mortgage + S_loan + L_loan + Imp_Open + FDI_Open + Dev + Cap + Infl + GDP_gr +  Stock_turn + CBI + M_policy + Ex_regime + Eurozone,  ssc=ssc(adj = FALSE, cluster.adj=FALSE), weights = ~Weight_est, cluster = "study_id", panel.id = "study_id",  data = Data)
predict(OLSCONF, newdata = new.dat, interval = 'confidence')

#Cut
new.dat <- data.frame(se=0, interbank=0.7763, hike=0, T_span=12.0446, PreCrises=0, PostCrises=1, Gov=0.8249, IMF=0.0796, ECB=0.0835, Quarterly=1, Monthly=0, Weekly=0, Time_Series=0.9446, Lag_length=0.6275, Long_run=0, Short_run=1, Infl_con=0.0435, OLS=0, Johansen=1, P_year=19.1277, A_Cite=22.65, Peer=0.8332, IMPF=0.0756, Central_banker=0.3156, Firms=0.2697, Households=0.0877, Mortgage=0.0572, S_loan=0.3047, L_loan=0.0685, Imp_Open=0.2473, FDI_Open=0.0272, Dev=0, Cap=0.7019, Infl=0.0538, GDP_gr=0.0396, Stock_turn=0.9999, CBI=0.5569, M_policy=0.4508, Ex_regime=0.7669, Eurozone=0.2124)
predict(bma_irpt, newdata = new.dat, interval = 'confidence')

#OLS confidence intervals
OLSCONF <- feols(beta ~ se + interbank + hike + T_span + PreCrises + PostCrises + Gov + IMF + ECB + Quarterly + Monthly + Weekly + Time_Series + Lag_length + Long_run + Short_run + Infl_con + OLS + Johansen + P_year + A_Cite + Peer + IMPF + Central_banker + Firms + Households + Mortgage + S_loan + L_loan + Imp_Open + FDI_Open + Dev + Cap + Infl + GDP_gr +  Stock_turn + CBI + M_policy + Ex_regime + Eurozone,  ssc=ssc(adj = FALSE, cluster.adj=FALSE), weights = ~Weight_est, cluster = "study_id", panel.id = "study_id",  data = Data)
predict(OLSCONF, newdata = new.dat, interval = 'confidence')



#Developing:Long-run
#Hike
new.dat <- data.frame(se=0, interbank=0.7763, hike=1, T_span=12.0446, PreCrises=0, PostCrises=1, Gov=0.8249, IMF=0.0796, ECB=0.0835, Quarterly=1, Monthly=0, Weekly=0, Time_Series=0.9446, Lag_length=0.6275, Long_run=1, Short_run=0, Infl_con=0.0435, OLS=0, Johansen=1, P_year=19.1277, A_Cite=22.65, Peer=0.8332, IMPF=0.0756, Central_banker=0.3156, Firms=0.2697, Households=0.0877, Mortgage=0.0572, S_loan=0.3047, L_loan=0.0685, Imp_Open=0.2473, FDI_Open=0.0272, Dev=0, Cap=0.7019, Infl=0.0538, GDP_gr=0.0396, Stock_turn=0.9999, CBI=0.5569, M_policy=0.4508, Ex_regime=0.7669, Eurozone=0.2124)
predict(bma_irpt, newdata = new.dat, interval = 'confidence')

#OLS confidence intervals
OLSCONF <- feols(beta ~ se + interbank + hike + T_span + PreCrises + PostCrises + Gov + IMF + ECB + Quarterly + Monthly + Weekly + Time_Series + Lag_length + Long_run + Short_run + Infl_con + OLS + Johansen + P_year + A_Cite + Peer + IMPF + Central_banker + Firms + Households + Mortgage + S_loan + L_loan + Imp_Open + FDI_Open + Dev + Cap + Infl + GDP_gr +  Stock_turn + CBI + M_policy + Ex_regime + Eurozone,  ssc=ssc(adj = FALSE, cluster.adj=FALSE), weights = ~Weight_est, cluster = "study_id", panel.id = "study_id",  data = Data)
predict(OLSCONF, newdata = new.dat, interval = 'confidence')

#Cut
new.dat <- data.frame(se=0, interbank=0.7763, hike=0, T_span=12.0446, PreCrises=0, PostCrises=1, Gov=0.8249, IMF=0.0796, ECB=0.0835, Quarterly=1, Monthly=0, Weekly=0, Time_Series=0.9446, Lag_length=0.6275, Long_run=1, Short_run=0, Infl_con=0.0435, OLS=0, Johansen=1, P_year=19.1277, A_Cite=22.65, Peer=0.8332, IMPF=0.0756, Central_banker=0.3156, Firms=0.2697, Households=0.0877, Mortgage=0.0572, S_loan=0.3047, L_loan=0.0685, Imp_Open=0.2473, FDI_Open=0.0272, Dev=0, Cap=0.7019, Infl=0.0538, GDP_gr=0.0396, Stock_turn=0.9999, CBI=0.5569, M_policy=0.4508, Ex_regime=0.7669, Eurozone=0.2124)
predict(bma_irpt, newdata = new.dat, interval = 'confidence')

#OLS confidence intervals
OLSCONF <- feols(beta ~ se + interbank + hike + T_span + PreCrises + PostCrises + Gov + IMF + ECB + Quarterly + Monthly + Weekly + Time_Series + Lag_length + Long_run + Short_run + Infl_con + OLS + Johansen + P_year + A_Cite + Peer + IMPF + Central_banker + Firms + Households + Mortgage + S_loan + L_loan + Imp_Open + FDI_Open + Dev + Cap + Infl + GDP_gr +  Stock_turn + CBI + M_policy + Ex_regime + Eurozone,  ssc=ssc(adj = FALSE, cluster.adj=FALSE), weights = ~Weight_est, cluster = "study_id", panel.id = "study_id",  data = Data)
predict(OLSCONF, newdata = new.dat, interval = 'confidence')






#  FMA  ###################################################################

#================================
# Mallows Model Averaging Program 
#================================

# Loading libraries
library(foreign)
library(xtable)
library(LowRankQP)
#==================================================
set.seed(42)
IRPT <- Data[, c("beta", "se", "interbank", "hike", "T_span", "PreCrises", "PostCrises", "Gov", "IMF", "ECB", "Quarterly", "Monthly", "Weekly", "Time_Series", "Lag_length", "Long_run", "Short_run", "Infl_con", "OLS", "Johansen", "P_year", "A_Cite", "Peer", "IMPF", "Central_banker", "Firms", "Households", "Mortgage", "S_loan", "L_loan", "Imp_Open", "FDI_Open", "Dev", "Cap", "Infl", "GDP_gr", "Stock_turn", "CBI", "M_policy", "Ex_regime", "Eurozone")]


#Construct Weights:
Data$Weight_est <- 1/Data$est_no

#Weight the variables:
IRPT <- IRPT*Data$Weight_est
IRPT <- na.omit(IRPT)

mydata <- IRPT

#adding constant
x.data <- mydata[,-1]
const_<-c(1)
x.data <-cbind(const_,x.data)


x <- sapply(1:ncol(x.data),function(i){x.data[,i]/max(x.data[,i])})   
scale.vector <- as.matrix(sapply(1:ncol(x.data),function(i){max(x.data[,i])}))        
Y <- as.matrix(mydata[,1])
output.colnames <- colnames(x.data)
full.fit <- lm(Y~x-1)
beta.full <- as.matrix(coef(full.fit))
M <- k <- ncol(x)
n <- nrow(x)
beta <- matrix(0,k,M)
e <- matrix(0,n,M)
K_vector <- matrix(c(1:M))
var.matrix <- matrix(0,k,M)
bias.sq <- matrix(0,k,M)             


# MMA Estimator using orthogonalization 
for(i in 1:M)
{
  X <- as.matrix(x[,1:i])
  ortho <- eigen(t(X)%*%X)
  Q <- ortho$vectors ; lambda <- ortho$values 
  x.tilda <- X%*%Q%*%(diag(lambda^-0.5,i,i))
  beta.star <- t(x.tilda)%*%Y
  beta.hat <- Q%*%diag(lambda^-0.5,i,i)%*%beta.star
  beta[1:i,i] <- beta.hat
  e[,i] <- Y-x.tilda%*%as.matrix(beta.star)
  bias.sq[,i] <- (beta[,i]-beta.full)^2
  var.matrix.star <- diag(as.numeric(((t(e[,i])%*%e[,i])/(n-i))),i,i)
  var.matrix.hat <- var.matrix.star%*%(Q%*%diag(lambda^-1,i,i)%*%t(Q))
  var.matrix[1:i,i] <- diag(var.matrix.hat)
  var.matrix[,i] <- var.matrix[,i]+ bias.sq[,i]
} # End loop over i

e_k <- e[,M]
sigma_hat <- as.numeric((t(e_k)%*%e_k)/(n-M))
G <- t(e)%*%e
a <- ((sigma_hat)^2)*K_vector
A <- matrix(1,1,M)
b <- matrix(1,1,1)
u <- matrix(1,M,1)
optim <- LowRankQP(Vmat=G,dvec=a,Amat=A,bvec=b,uvec=u,method="LU",verbose=FALSE)
weights <- as.matrix(optim$alpha)
beta.scaled <- beta%*%weights
final.beta <- beta.scaled/scale.vector
std.scaled <- sqrt(var.matrix)%*%weights
final.std <- std.scaled/scale.vector
results.reduced <- as.matrix(cbind(final.beta,final.std))
rownames(results.reduced) <- output.colnames; colnames(results.reduced) <- c("Coefficient", "Sd. Err")
#MMA.fls <- round(results.reduced,4)[-1,]
MMA.fls <- round(results.reduced,4)
list(MMA.fls)


MMA.fls <- data.frame(MMA.fls)
t <- as.data.frame(MMA.fls$Coefficient/MMA.fls$Sd..Err)
MMA.fls$pv <-round( (1-apply(as.data.frame(apply(t,1,abs)), 1, pnorm))*2,3)
MMA.fls$names <- rownames(MMA.fls)
names <- c(colnames(mydata))
names <- c(names,"const_")
MMA.fls <- MMA.fls[match(names, MMA.fls$names),]
MMA.fls$names <- NULL
MMA.fls



xtable(MMA.fls,digits=c(0,3,3,3))




#Frequentist check
IRPT <- Data[, c("beta", "se", "interbank", "hike", "T_span", "PreCrises", "PostCrises", "Gov", "IMF", "ECB", "Quarterly", "Monthly", "Weekly", "Time_Series", "Lag_length", "Long_run", "Short_run", "Infl_con", "OLS", "Johansen", "P_year", "A_Cite", "Peer", "IMPF", "Central_banker", "Firms", "Households", "Mortgage", "S_loan", "L_loan", "Imp_Open", "FDI_Open", "Dev", "Cap", "Infl", "GDP_gr", "Stock_turn", "CBI", "M_policy", "Ex_regime", "Eurozone")]


#Construct Weights:
IRPT$Weight_est <- 1/(Data$est_no*Data$est_no)
IRPT$study_id <- Data$study_id

freq_ols <-feols(beta ~ se + hike + Quarterly + IMF + PreCrises + Long_run + Short_run + OLS + Johansen + A_Cite + Central_banker + Imp_Open + FDI_Open + Dev + Infl, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster = "study_id", weights = ~Weight_est, data = IRPT)

summary(freq_ols)