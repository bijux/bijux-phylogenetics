#***********************************************************************************
#BIOR90: Phylogenetic Comparative Analyses ----
#Trait evolution by Charlie Cornwallis
#PCM4 Script to accompany lecture
#***********************************************************************************

#***********************************************************************************
#Packages ----
#***********************************************************************************
pacman::p_load(tidyverse,RColorBrewer,ape,MCMCglmm,picante,geiger,phytools,ggtree)

#***********************************************************************************
#* Bayesian Linear Models ----
#***********************************************************************************
#A new data example - Cichlids
fish <- read.csv("Data/cichlid.csv")
fishtree<-read.nexus("Data/cichlidtree.nex")

#Check how the data looks
head(fish)

#Are all species in the dataset in the tree?
table(fish$treetip %in% fishtree$tip.label)

#Are all tree tips in the dataset?
table(fishtree$tip.label %in% fish$treetip)

#tree pruning
missingfish<-fishtree$tip.label[fishtree$tip.label %in% fish$treetip == "FALSE"]
missingfish
fishtree<-drop.tip(fishtree,missingfish)

#Make a summary dataset per species
fishSum<- fish %>% group_by(treetip,species,mating_system) %>% summarise(across(c(mass,gut_length,habitat_complexity:prop_spp),
                                                                                .fns = ~mean(.x,na.rm=T), #fns = function = mean. na.rm=T removes missing values
                                                                                .names = "{col}"))#This could be used to rename columns if you wanted
#Lets save these fish files in one RData file that we can use later on
save(fishtree,fish, fishSum,file = "Data/fishdata.RData")

#***********************************************************************************
#* Non-phylogenetic linear regression fit with gls and MCMCglmm ----
#Does gut length predict how many parasites abundance?
#***********************************************************************************
hist(fishSum$parasite_abundance)
hist(log(fishSum$parasite_abundance+1))
ggplot(data=fishSum, aes(x=gut_length, y=log(parasite_abundance+1)))+
       geom_point()+
       geom_smooth(method="lm")


#gls fit
glsFish <- gls(log(parasite_abundance+1) ~ gut_length, data=fishSum,method="ML")
summary(glsFish)
plot(glsFish)

#MCMCglmm fit
#Setting the priors: these are the defaults
p1 = list(B=list(mu=rep(0,2), V=diag(2)*1e+8), R=list(V=1,nu=0.002))

#Specify the model
m1<-MCMCglmm(log(parasite_abundance+1) ~ gut_length, family ="gaussian", data = fishSum,prior=p1, nitt=110000, burnin=10000, thin=100,verbose=F)
#returns an error about 'family' - what does this mean?

#MCMCglmm does not like tibble tables so put as.data.frame() around the dataset
fishSum = as.data.frame(fishSum)
m1<-MCMCglmm(log(parasite_abundance+1) ~ gut_length, family ="gaussian", data = as.data.frame(fishSum),prior=p1, nitt=110000, burnin=10000, thin=100,verbose=F)


#***********************************************************************************
#* Parameter estimation & Hypothesis testing ----
#***********************************************************************************
#What is the intercept and slope of the regression line of gut length vs parasite abundance?
summary(m1a)
posterior.mode(m1a$Sol)

#And the credible intervals around both estimates?
HPDinterval(m1a$Sol)

#Can do the same for random effects
posterior.mode(m1a$VCV)
HPDinterval(m1a$VCV)

#Hypothesis testing - do the credible intervals span 0, if not then significantly different from 0
##Can be backed up with summary() which gives a Bayesian version of a p value
summary(m1a)

##Can also do similar custom tests using % of iterations different from a specified value / total iterations e.g.
table(m1a$Sol[,'gut_length']>0)/1000


#***********************************************************************************
#* Diagnostics ----
#***********************************************************************************

##Checking the run of a model
###Check the convergence of the fixed effects ('Sol' object of MCMCglmm output)
plot(m1$Sol)
autocorr(m1$Sol)

###Check the convergence of the random effects ('VCV' object of MCMCglmm output)
plot(m1$VCV)
autocorr(m1$VCV)

##Checking consistency across runs
m1a<-MCMCglmm(log(parasite_abundance+1) ~ gut_length, family ="gaussian", data = fishSum,prior=p1, nitt=110000, burnin=10000, thin=100,verbose=F)
m1b<-MCMCglmm(log(parasite_abundance+1) ~ gut_length, family ="gaussian", data = fishSum,prior=p1, nitt=110000, burnin=10000, thin=100,verbose=F)
m1c<-MCMCglmm(log(parasite_abundance+1) ~ gut_length, family ="gaussian", data = fishSum,prior=p1, nitt=110000, burnin=10000, thin=100,verbose=F)

###Plot the fixed effects
plot(mcmc.list(m1a$Sol,m1b$Sol,m1c$Sol))
###And random effects
plot(mcmc.list(m1a$Sol,m1b$Sol,m1c$Sol))

###Test whether between chain variance is less than within chain variance = convergence (PSR <1.1 acceptable))
gelman.diag(mcmc.list(m1a$Sol,m1b$Sol,m1c$Sol))

#***********************************************************************************
# Bayesian Phylogenetic Mixed Models ----
#***********************************************************************************
#Setting the priors: these are the defaults
p2 = list(B=list(mu=rep(0,2), V=diag(2)*1e+8), G=list(G1=list(V=1,nu=0.002)),R=list(V=1,nu=0.002))

#Specify the model
##Need to convert the phylogeny into a covariance matrix to feed to MCMCglmm 
inv.fishphylo<-inverseA(fishtree,nodes="TIPS",scale=TRUE)

##Run the model
m2<-MCMCglmm(log(parasite_abundance+1) ~ gut_length, random=~treetip, ginverse=list(treetip=inv.fishphylo$Ainv), family ="gaussian", data = fishSum,prior=p2, nitt=110000, burnin=10000, thin=100,verbose=F)

#Check the model
plot(m2$Sol)
plot(m2$VCV)

#Parameter estimates for the fixed effects
summary(m2)
posterior.mode(m2$Sol)
HPDinterval(m2$Sol)

#Parameter estimates for the random effects
posterior.mode(m2$VCV)
HPDinterval(m2$VCV)

#*********************
#* Compare to PGLS ----
#*********************
pgls6 <- gls(log(parasite_abundance+1) ~ gut_length, data=fishSum, correlation=corBrownian(1, fishtree, form = ~treetip), method="ML")
summary(pgls6)


#***********************************************************************************
#END ####
#***********************************************************************************
