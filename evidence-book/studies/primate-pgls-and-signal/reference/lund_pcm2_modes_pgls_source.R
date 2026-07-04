#***********************************************************************************
#BIOR90: Phylogenetic Comparative Analyses ----
#Trait evolution by Charlie Cornwallis
#PCM2 Script to accompany lecture
#***********************************************************************************

#***********************************************************************************
#Packages ----
#***********************************************************************************
pacman::p_load(tidyverse,RColorBrewer,ape,MCMCglmm,picante,geiger,phytools,ggtree)

#***********************************************************************************
#Read in test data again: Primates ----
#***********************************************************************************
#Read in data 
load("Data/primate.RData") #The advantage of RData objects - everything we need with one line of code

#***********************************************************************************
#1. Modes of evolution: BM, OU, ACDC ----
#***********************************************************************************
#See how the tree transformation change the tree
##Transform the trees
primatetreeOU1<-rescale(primatetree, "OU", alpha=1)
primatetreeOU10<-rescale(primatetree, "OU", alpha=10)
primatetreeEB<-rescale(primatetree, "EB", a=2)
primatetreeLB<-rescale(primatetree, "EB", a=-2)

##Plot the trees
par(mfrow=c(1,5))
plot(primatetree, cex=0.6, show.tip.label = F);title("Real phylogeny", cex.main=1)
plot(primatetreeOU1, cex=0.6, show.tip.label = F);title("OU 1", cex.main=1)
plot(primatetreeOU10,cex=0.6, show.tip.label = F);title("OU 10", cex.main=1)
plot(primatetreeEB,cex=0.6, show.tip.label = F);title("EB (a = 2)", cex.main=1)
plot(primatetreeLB,cex=0.6, show.tip.label = F);title("LB (a = -2)", cex.main=1)

#Fitting different models
#Formate data for Geiger?
primateLL<-as.numeric(primate$longevity[match(primatetree$tip.label,primate$species)])
names(primateLL)<-primatetree$tip.label

BM<-fitContinuous(phy= primatetree, dat = primateLL,model = "BM")
OU<-fitContinuous(phy= primatetree, dat = primateLL,model = "OU")
OU<-fitContinuous(phy= primatetree, dat = primateLL,model = "OU")
OU$opt$alpha

#Restrict the bounds of alpha
OU<-fitContinuous(phy= primatetree, dat = primateLL,model = "OU", bounds=list(alpha=c(0,10)))
OU$opt$alpha

EB<-fitContinuous(phy= primatetree, dat = primateLL,model = "EB")

#Restrict the bounds of a
EB<-fitContinuous(phy= primatetree, dat = primateLL,model = "EB",bounds=list(a=c(0,50)))
EB$opt$a

#Compare model fit statistics
BM$opt$aic
OU$opt$aic
EB$opt$aic

#Formally test differences between models
##BM versus OU
BMvsOU <- -2*(BM$opt$lnL - OU$opt$lnL)
BMvsOU
pchisq(BMvsOU, df=1,lower.tail = FALSE)

##BM versus EB
BMvsEB <- -2*(BM$opt$lnL - EB$opt$lnL)
BMvsEB
pchisq(BMvsEB, df=1,lower.tail = FALSE)

##OU versus EB
OUvsEB <- -2*(OU$opt$lnL - EB$opt$lnL)
OUvsEB
pchisq(OUvsEB, df=1,lower.tail = FALSE)

#Plot best fitting trees and compare
primatetreeOU5.3<-rescale(primatetree, "OU", alpha=5.335655)
primatetreeEB10.6<-rescale(primatetree, "EB", a=10.67126)

par(mfrow=c(1,3),mar=c(2,2,2,2))
plot(primatetree, show.tip.label=F)
tiplabels(pch=21, bg="grey38",cex=primateLL/300);title("Real Phylogeny", cex.main=2)
plot(primatetreeOU5.3,cex=0.6, show.tip.label = F)
tiplabels(pch=21, bg="grey38",cex=primateLL/300);title("OU (alpha = 5.3)", cex.main=2)
plot(primatetreeEB10.6,cex=0.6, show.tip.label = F)
tiplabels(pch=21, bg="grey38",cex=primateLL/300);title("EB (a = 10.6)", cex.main=2)

#***********************************************************************************
#2. Ancestral state reconstructions in light of modes of evolution ----
#***********************************************************************************
#Reconstruct states on transformed tree
long_ace_EB<-ace(primate$longevity,primatetreeEB10.6,type="continuous",method="pic")

#Reconstructions assuming BM
long_ace<-ace(primate$longevity,primatetree,type="continuous",method="pic")

par(mfrow=c(1,2))
plot(primatetree, cex=0.4);title("Real phylogeny", cex=28)
nodelabels(pch=21, bg="grey38",cex = long_ace$ace/300)

plot(primatetreeEB10.6,cex=0.4);title("EB", cex=28)
nodelabels(pch=21, bg="grey38",cex = long_ace_EB$ace/300)

#How this influence some of the basal nodes
long_ace_EB$ace[1:5]
long_ace$ace[1:5]

#And the more recent nodes
long_ace_EB$ace[70:74]
long_ace$ace[70:74]

#***********************************************************************************
#3. Estimating Relationships Between Variables – linear models ----
#***********************************************************************************
#*******************************************
#Does sociality predict variation in longevity?
#Linear regression
#*******************************************
ggplot(data=primate, aes(x=social_group_size, y=longevity)) + geom_point() + geom_smooth(method='lm')#Looks a bit ropey

#Fit a non-phylogenetic regression
gls1 <- gls(longevity ~ social_group_size, data=primate,method="ML")
summary(gls1)

#Check for heteroscedasticity
plot(gls1, resid(., type="n")~fitted(.), main="Normalized Residuals vs. Fitted Values", abline=c(0,0))

#Check for departures from normal distribution of residuals
par(mfrow=c(1,1))
res1 <- resid(gls1, type="n")
qqnorm(res1)
qqline(res1)

#Hypothesis testing
anova(gls1)

#***********************************************************************************
#4. Relationships Between Variables – phylogenetic linear models ----
#***********************************************************************************
#*******************************************
#* Lets first see if residuals are correlated with respect to phylogeny (one of the main motivations of doing pgls) -----
#*******************************************
primate$predicted <- predict(gls1)   # Save the predicted values
primate$residuals <- residuals(gls1) # Save the residual values

ggplot(data=primate, aes(x=social_group_size, y=longevity)) + 
  geom_smooth(method = "lm", se = T) +  # Plot regression slope
  geom_segment(aes(xend = social_group_size, yend = predicted)) + 
  geom_point()

#Visualise residuals on phylogeny
primateResid<-primate$residuals[match(primatetree$tip.label,primate$species)]#make sure residuals are in the same order as the tips
names(primateResid)<-primatetree$tip.label
plotTree.wBars(primatetree,primateResid,scale = 0.001)

#*******************************************
#* Fit a phylogenetic regression ----
#*******************************************
#Use correlated residuals statement to fit phylogenetically correlated errors. First we can show that fixing correlations at 0 is the same as gls
pgls1 <- gls(longevity ~ social_group_size, data=primate, correlation=corPagel(0, primatetree,fixed=TRUE,form = ~species),method="ML")
summary(pgls1)

#Next we can fit a phylogenetic regression where lambda is estimated in the response variable where fixed = FALSE
pgls2 <- gls(longevity ~ social_group_size, data=primate, correlation=corPagel(0, primatetree, fixed=FALSE,form = ~species),method="ML")
summary(pgls2)

#Check for heteroscedasticity
plot(pgls2, resid(., type="n")~fitted(.), main="Normalized Residuals vs. Fitted Values", abline=c(0,0))

#Check for departures from normal distribution of residuals
respgls2 <- resid(pgls2, type="n")
qqnorm(respgls2)
qqline(respgls2)

#Plot the fit of the gls and pgls
plot(primate$social_group_size, primate$longevity, pch=19, col="black", xlab="Social group size", ylab="Longevity", main="Compare the fits: PGLS (blue) vs. OLS (red)")
abline(gls1, col="red")
abline(pgls2, col="blue")

#***********************************************************************************
#5. Using PGLS to examine phylogenetic signal ----
#***********************************************************************************
pgls3 <- gls(longevity ~ 1, data=primate, correlation=corPagel(0, primatetree, fixed=FALSE, form = ~species),method="ML")
summary(pgls3)

#*******************************************
#* Testing if Pagel's Lambda is different from 0 ----
#*******************************************
pgls3 <- gls(longevity ~ 1, data=primate, correlation=corPagel(0, primatetree, fixed=FALSE, form = ~species),method="ML")
pgls4 <- gls(longevity ~ 1, data=primate, correlation=corPagel(0, primatetree, fixed=TRUE, form = ~species),method="ML")
anova(pgls3,pgls4)

#*******************************************
#* Fitting models with different modes of evolution ----
#*******************************************
#Lets fit equivalent model to BM<-fitContinuous(phy= primatetree, dat = primateLL,model = "BM")
pglsBM <- gls(longevity ~ 1, data=primate, correlation=corBrownian(1, primatetree, form = ~species),method="ML")
summary(pglsBM)

#Lets fit equivalent model to OU<-fitContinuous(phy= primatetree, dat = primateLL,model = "OU")
pglsOU <- gls(longevity ~ 1, data=primate, correlation=corMartins(1, primatetree,fixed=FALSE, form = ~species),method="ML")
summary(pglsOU)

#Test whether OU is better fit than BM model
anova(pglsBM,pglsOU)

#Lets fit equivalent model to EB<-fitContinuous(phy= primatetree, dat = primateLL,model = "EB")
pglsEB <- gls(longevity ~ 1, data=primate, correlation=corBlomberg(1, primatetree,fixed=FALSE, form = ~species),method="ML")


#*******************************************
#* Finding the maximum likelihood by fixing the value of delta ----
#*******************************************
#define the range for delta and at which intervals we are estimating the likelihood
a <- seq(0.1, 20, 0.1)

#Calculate likelihood over the entire range of delta
likEB <- sapply(a, function(a) logLik(gls(longevity ~ 1, data=primate, correlation = corBlomberg(value=a, phy=primatetree, fixed = TRUE, form = ~species))))
plot(likEB ~ a, type = "l", main = "Likelihood Plot for a", ylab = "Log Likelihood", xlab = "a")


#***********************************************************************************
#END ####
#***********************************************************************************


