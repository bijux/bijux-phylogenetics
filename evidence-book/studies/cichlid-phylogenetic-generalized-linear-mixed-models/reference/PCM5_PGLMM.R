#***********************************************************************************
#BIOR90: Phylogenetic Comparative Analyses ----
#Trait evolution by Charlie Cornwallis
#Day5 Script to accompany lecture
#***********************************************************************************

#***********************************************************************************
#Packages ----
#***********************************************************************************
pacman::p_load(tidyverse,RColorBrewer,ape,MCMCglmm,picante,geiger,phytools,ggtree,boot)

#***********************************************************************************
#Read in fish data that you saved yesterday ----
#***********************************************************************************
#Read in data
load("Data/fishdata.RData")

#remove formatting as MCMCglmm doesn't like it
fish<-as.data.frame(fish)
fishSum<-as.data.frame(fishSum)

#***********************************************************************************
#1. Phylogenetic heritability (means the same thing as phylogenetic signal)----
#***********************************************************************************
#MCMCglmm fit
#Setting the priors: these are the defaults
p1 = list(R=list(V=1,nu=0.002))

#Non-phylogenetic model from yesterday
m1a<-MCMCglmm(log(parasite_abundance+1) ~ gut_length, family ="gaussian", data = fishSum,prior=p1, nitt=110000, burnin=10000, thin=100,verbose=F)

#Phylogenetic model from yesterday
#Priors
p2 = list(G=list(G1=list(V=1,nu=0.002)),R=list(V=1,nu=0.002))

inv.fishphylo<-inverseA(fishtree,nodes="TIPS",scale=TRUE)

m2<-MCMCglmm(log(parasite_abundance+1) ~ gut_length, random=~treetip, ginverse=list(treetip=inv.fishphylo$Ainv), family ="gaussian", data = fishSum,prior=p2, nitt=110000, burnin=10000, thin=100,verbose=F)

#Does model fit better than without phylogenetic effect
m1a$DIC
m2$DIC

m2PhyloSig<-m2$VCV[,'treetip']/(m2$VCV[,'treetip']+m2$VCV[,'units'])

#Check it's estimate converged
plot(m2PhyloSig)

#Calculate it's point estimate and credible interval
posterior.mode(m2PhyloSig)
HPDinterval(m2PhyloSig)

#***********************************************************************************#
#* Ancestral state reconstructions with MCMCglmm ----
#***********************************************************************************#
#Specify the model
##Need to convert the phylogeny into a covariance matrix to feed to MCMCglmm 
inv.fishphylo<-inverseA(fishtree)

#Priors
p3 = list(B=list(mu=rep(0,1), V=diag(1)*1e+8), G=list(G1=list(V=1,nu=0.002)),R=list(V=1,nu=0.002))

m3<-MCMCglmm(log(parasite_abundance+1) ~ 1, random=~treetip,ginverse=list(treetip=inv.fishphylo$Ainv), pr=TRUE,family ="gaussian", data = fishSum,prior=p3, nitt=110000, burnin=10000, thin=100,verbose=F)

plot(m3$Sol)
plot(m3$VCV)

#Take a look at the estimates - what do they mean?
posterior.mode(m3$Sol)

#***********************************************************************************
#* Calculating and exporting ancestral values to dataset ----
#***********************************************************************************
#There are 33 nodes - 1 for the root which is the intercept 
#so the first 32 are the deviations (residuals) of the nodes from the mean (also known as BLUPs)
#To get estimates of ancestral states add intercept to node BLUPs
blupsm3<-data.frame(effect=colnames(m3$Sol), estimate=posterior.mode(m3$Sol[,'(Intercept)'])+posterior.mode(m3$Sol),CI=HPDinterval(m3$Sol[,'(Intercept)']+m3$Sol))

#We added the intercept value to all values and so "Node1" value is incorrect - lets replace it
blupsm3['(Intercept)','estimate']<-posterior.mode(m3$Sol[,'(Intercept)'])
blupsm3['(Intercept)',c('CI.lower','CI.upper')]<-HPDinterval(m3$Sol[,'(Intercept)'])

#The intercept is the value of the root ("Node1") so lets rename it
blupsm3$effect<-as.character(blupsm3$effect)
blupsm3['(Intercept)','effect']<-"Node1"

#and remove the text "treetip." from the effect
blupsm3$effect<-gsub("treetip.","",blupsm3$effect)

#***********************************************************************************
#* Plot estimates for nodes and species values ----
#***********************************************************************************
# Convert 'effect' to a factor with levels in the same order as the dataframe
blupsm3$effect <- factor(blupsm3$effect, levels = unique(blupsm3$effect))

ggplot(blupsm3, aes(x=effect, y=estimate)) +
  geom_point() +  # Add points for each estimate
  geom_hline(yintercept=blupsm3['(Intercept)','estimate']) +  # Add vertical line at node1 / intercept
  geom_segment(aes(yend=blupsm3['(Intercept)','estimate']), 
               x=blupsm3$effect, y=blupsm3$estimate)+ 
  theme(axis.text.x = element_text(angle = 65, hjust = 1, vjust = 1))
  

#***********************************************************************************
#* Plotting Ancestral state reconstructions on tree ----
#***********************************************************************************
#Plot values on ancestral reconstructions on the phylogeny

#MCMCglmm uses the same naming procedure of nodes as ape so plotting will be in the same order 
#Here's a messy way to check this
#Put node labels on the tree so that we can make up values
fishParaAbund<-log(fish$parasite_abundance+1)[match(fishtree$tip.label,fish$treetip)]

plot(fishtree, cex=0.8, no.margin =T, label.offset = 0.7)
nodelabels(pch=21, cex=blupsm3$estimate[1:fishtree$Nnode])
tiplabels(pch=21,cex=fishParaAbund/2,bg="black")

#***********************************************************************************
#2. Bayesian Generalised Linear Mixed Models: parasite abundance should be modelled as a Poisson process ----
#***********************************************************************************
#* Poisson ----
#Set the priors
p5 = list(B=list(mu=rep(0,1), V=diag(1)*1e+8), G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.002))

#We are going to model this as a Poisson process (suitable for count data)
#These need to be whole numbers and as we've taken the mean of of parasite_abundance in the fishSum data we need to round it
fishSum$parasite_abundance<-round(fishSum$parasite_abundance,0)

m5<-MCMCglmm(parasite_abundance ~ 1, random=~treetip, ginverse=list(treetip=inv.fishphylo$Ainv), pr=TRUE,family ="poisson", data = fishSum,prior=p5, nitt=110000, burnin=10000, thin=100,verbose=F)


#***********************************************************************************
#* Binomial traits ----
#***********************************************************************************
#MCMCglmm models proportions as number of events that A happen versus number of events where B happened (must be whole numbers)
##We will model the proportion of parasite species (prop_spp) infecting species
##To do this we need to convert the data into number of species infecting and number of species not infecting
##Number of species infecting =
fishSum$infected<-round(fishSum$prop_spp*100,0)
fishSum$not_infected<-100-fishSum$infected

#Setting the priors: lets keep with the defaults
p6 = list(B=list(mu=rep(0,1), V=diag(1)*1e+8), G=list(G1=list(V=1,nu=0.002)),R=list(V=1,nu=0.002))

#Specify the model
#The term 'family' allows the error distribution to be defined. If there is a variable called family in your dataset you will get an error message
#family= multinomial2 allows percentage data e.g. (10 out of 50) not to be confused with family=categorical which allows both binary (0,1) and mutliple discrete categories (red, blue, yellow) to be modelled
m6<-MCMCglmm(cbind(infected,not_infected) ~ 1, random=~treetip, ginverse=list(treetip=inv.fishphylo$Ainv),family ="multinomial2", data = fishSum,prior=p6, nitt=110000, burnin=10000, thin=100,verbose=F)

#Check the model
par(mar=c(2,2,2,2))
plot(m6$Sol)
plot(m6$VCV)

#Parameter estimates for the fixed effects: these are on the logit scale
posterior.mode(m6$Sol)
HPDinterval(m6$Sol)

#We can backtransform these to get estimates on the probability scale
inv.logit(posterior.mode(m6$Sol))
inv.logit(HPDinterval(m6$Sol))

#How does the intercept compare to the mean % of parasite species
mean(fish$prop_spp)

#Parameter estimates for the random effects
posterior.mode(m6$VCV)
HPDinterval(m6$VCV)

#For non-Gaussian traits we need to calculate the intra-class correlation coefficient (ICC) = associated between observations.
#Same as before but with distribution variance. For logit link function (binomial data) this = pi^2/3
m6PhyloSig<-m6$VCV[,'treetip']/((m6$VCV[,'treetip']+m6$VCV[,'units'])+pi^2/3)
plot(m6PhyloSig)

posterior.mode(m6PhyloSig)
HPDinterval(m6PhyloSig)

#***********************************************************************************
#Back-transforming continuous fixed effects 
#***********************************************************************************
#Lets fit a model with a continous fixed effects: Is having a broader diet associated with more parasite species?
#Setting the priors
p7 = list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.002))

m7<-MCMCglmm(cbind(infected,not_infected) ~ diet_breadth, random=~treetip, ginverse=list(treetip=inv.fishphylo$Ainv),family ="multinomial2", data = fishSum,prior=p7, nitt=110000, burnin=10000, thin=100,verbose=F)

#Parameter estimates on the logit scale
posterior.mode(m7$Sol)
HPDinterval(m7$Sol)

#We can backtransform these to get estimates on the probability scale

# Posterior samples
intercept <- m7$Sol[, "(Intercept)"]
diet_breadth <- m7$Sol[, "diet_breadth"]

# For the intercept

# Predicted intercept when diet_breadth = 0
x0 <- 0

# Predicted probability at diet_breadth = 0
p_intercept <- boot::inv.logit(intercept + diet_breadth * x0)

# Since x0 = 0, this is equivalent to:
# p_intercept <- boot::inv.logit(intercept)

# Posterior mode and HPD interval on the probability scale
posterior.mode(as.mcmc(p_intercept))
HPDinterval(as.mcmc(p_intercept))

# For the slope of diet breadth

# Choose the value of diet_breadth where you want the slope
# If diet_breadth is centred/scaled, x0 = 0 is usually the mean
x0 <- 0

# Predicted probability at x0
p_x <- boot::inv.logit(intercept + diet_breadth * x0)

# Predicted probability after a 1-unit increase in diet_breadth
p_x_plus_1 <- inv.logit(intercept + diet_breadth * (x0 + 1))

# Probability-scale slope: change in probability for 1-unit increase
slope_prob <- p_x_plus_1 - p_x

posterior.mode(as.mcmc(slope_prob))
HPDinterval(as.mcmc(slope_prob))

#***********************************************************************************#***********************************************************************************##
#3. Fitting additional random effects ----
#***********************************************************************************#***********************************************************************************##
#Before we had to summarise the dataset (1 observation per species), but now we can analyse repeated measures within species as well
fish$infected<-round(fish$prop_spp*100,0)
fish$not_infected<-100-fish$infected

#But in mixed models we don't have to do as they are more flexible than pgls
#We can fit an extra random effect of species to account for non-independence of observations
p10 = list(G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002)),R=list(V=1,nu=0.002))

m10<-MCMCglmm(cbind(infected,not_infected) ~ diet_breadth, random=~species+treetip, ginverse=list(treetip=inv.fishphylo$Ainv), family ="multinomial2", data = fish,prior=p10, nitt=110000, burnin=10000, thin=100,verbose=F)

#Fixed effect estimates
posterior.mode(m10$Sol)
HPDinterval(m10$Sol)

#Random effect estimates
posterior.mode(m10$VCV)
HPDinterval(m10$VCV)

#***********************************************************************************
#* Calculating phylogenetic heritability in models with multiple random effects
#***********************************************************************************
m10PhyloSig<-m10$VCV[,'treetip']/(m10$VCV[,'species']+m10$VCV[,'treetip']+m10$VCV[,'units'])

#Check it's estimate converged
plot(m10PhyloSig)

#Calculate it's point estimate and credible interval
posterior.mode(m10PhyloSig)
HPDinterval(m10PhyloSig)

#***********************************************************************************
#* Can also calculate amount of variation explained by species
#***********************************************************************************
m10spSig<-m10$VCV[,'species']/(m10$VCV[,'species']+m10$VCV[,'treetip']+m10$VCV[,'units'])

#Check it's estimate converged
plot(m10spSig)

#Calculate it's point estimate and credible interval
posterior.mode(m10spSig)
HPDinterval(m10spSig)

#***********************************************************************************
#END ####
#***********************************************************************************
