#***********************************************************************************
#BIOR90: Phylogenetic Comparative Analyses ----
#Trait evolution by Charlie Cornwallis
#PCM3 Script to accompany lecture
#***********************************************************************************

#***********************************************************************************
#Packages ----
#***********************************************************************************
pacman::p_load(tidyverse,RColorBrewer,ape,MCMCglmm,picante,geiger,phytools,ggtree)

#***********************************************************************************
#Read in test data again: Primates ----
#***********************************************************************************
#Read in data
load("Data/primate.RData")

#***********************************************************************************
#1. Analysing multiple variables simultaneously - more on linear models ----
#Does sociality and sexual dimorphism predict variation in longevity? 
#***********************************************************************************
#Visualise relationships
ggplot(data=primate, aes(x=sex_dimorphism, y=longevity))+
  geom_boxplot()

#***********************************************************************************     
#* Fit Linear model with continuous and categorical explanatory variables ----
#***********************************************************************************
gls2<-gls(longevity ~ social_group_size+sex_dimorphism, data=primate,method="ML")

#Check for heteroscedasticity
plot(gls2, resid(., type="n")~fitted(.), main="Normalized Residuals vs. Fitted Values", abline=c(0,0))

#Check for departures from normal distribution of residuals
res2 <- resid(gls2, type="n")
qqnorm(res2)
qqline(res2)

#Check for heteroscedasticity in relation to factors
hist(gls2$residuals[primate$sex_dimorphism == "0"],col="blue",xlim=c(-450,450),ylim=c(0,20),xlab="",main="")
par(new=T)
hist(gls2$residuals[primate$sex_dimorphism == "1"],col="red",xlim=c(-450,450), ylim=c(0,20),xlab="residuals",main="")

#***********************************************************************************
#* Log-transform response and re-run model ----
#***********************************************************************************
gls3<-gls(log(longevity) ~ social_group_size+sex_dimorphism, data=primate,method="ML")

#Check for heteroscedasticity
plot(gls3, resid(., type="n")~fitted(.), main="Normalized Residuals vs. Fitted Values", abline=c(0,0))

#Check for departures from normal distribution of residuals
res3 <- resid(gls3, type="n")
qqnorm(res3)
qqline(res3)

#Check for heteroscedasticity in relation to factors
hist(gls3$residuals[primate$sex_dimorphism == "0"],col="blue",xlim=c(-1,1),ylim=c(0,20),xlab="",main="")
par(new=T)
hist(gls3$residuals[primate$sex_dimorphism == "1"],col="red",xlim=c(-1,1), ylim=c(0,20),xlab="residuals",main="")

summary(gls3)

#Hypothesis testing
anova(gls3, type ="sequential")
anova(gls3, type ="marginal")

#***********************************************************************************
#* Remove the intercept to get absolute estimates and SEs ----
#***********************************************************************************
gls4<-gls(log(longevity) ~ sex_dimorphism-1+social_group_size, data=primate,method="ML")
summary(gls4)

#***********************************************************************************
#* Linear model: interactions ----
# Is the relationship between group size and longevity different for sexually monomorphic and dimorphic species?
#***********************************************************************************
#Visualise the relationship
ggplot(data=primate, aes(x=social_group_size, y=log(longevity), group=sex_dimorphism,fill=sex_dimorphism,color=sex_dimorphism)) + 
  geom_point() + geom_smooth(method='lm')+theme(legend.position="top")

#Run the model
gls5<-gls(log(longevity) ~ social_group_size+sex_dimorphism+social_group_size:sex_dimorphism, data=primate,method="ML")
summary(gls5)

anova(gls5, type ="marginal")
anova(gls5, type ="sequential")

#***********************************************************************************
#* Estimate separate slopes for sexually dimorphic and monomorphic species ----
#***********************************************************************************
gls6<-gls(log(longevity) ~ sex_dimorphism+social_group_size:sex_dimorphism-1, data=primate,method="ML")
summary(gls6)

#***********************************************************************************
#* An aside about scaling and intercepts ----
#***********************************************************************************
#Convert social group size to standard normal
primate$Zsocial_group_size<-scale(primate$social_group_size)
ggplot(data=primate, aes(x=Zsocial_group_size, y=longevity, group=sex_dimorphism,fill=sex_dimorphism,color=sex_dimorphism)) + 
  geom_point() + geom_smooth(method='lm')+theme(legend.position="right")

gls7 <- gls(log(longevity) ~ Zsocial_group_size+sex_dimorphism+Zsocial_group_size:sex_dimorphism, data=primate,method="ML")
summary(gls7)

anova(gls7, type ="marginal")

#***********************************************************************************
#* PGLS: continuous and categorical explanatory variables & interactions ----
#***********************************************************************************
pgls5 <- gls(log(longevity) ~ Zsocial_group_size+sex_dimorphism+Zsocial_group_size:sex_dimorphism, data=primate, correlation=corMartins(0.5, primatetree,fixed=FALSE, form = ~species),method="ML")
summary(pgls5)

#Check diagnostics
#Check for heteroscedasticity
plot(pgls5, resid(., type="n")~fitted(.), main="Normalized Residuals vs. Fitted Values", abline=c(0,0))
respgls5 <- resid(pgls5, type="n")
qqnorm(respgls5)
qqline(respgls5)

#Hypothesis testing
anova(pgls5)

#***********************************************************************************
#END ####
#***********************************************************************************
