#***********************************************************************************
#Packages
#***********************************************************************************
pacman::p_load(ggplot2,openxlsx,tidyverse,RColorBrewer,ape,MCMCglmm,picante,geiger,phytools,ggtree, treeio,ggimage,cowplot,lmerTest,boot)

#Data
birddata = read.csv("Data/birds.csv", header = TRUE, sep = ";")
birdtree<-read.nexus("Data/birdtree.nex")

#***********************************************************************************
#An aside on log transformations and back-transformations ----
#***********************************************************************************
# Create a dummy dataset
x <- 1:100  # A sequence from 1 to 100
# Generate y values with an exponential relationship to x
y <- exp(0.05 * x)
# Calculate the natural logarithm of y
log_y <- log(y)
# Combine into a data frame
data <- data.frame(x, y, log_y)

#The model of testes mass
pgls_mp_vs_tm = gls(log(testes_mass) ~multiple_paternity, data = birddata, correlation=corPagel
                    (0, birdtree,fixed=FALSE,form = ~species),method="ML")

summary(pgls_mp_vs_tm)

#Data on log scale
#Change in testes mass from 0 to 10%
logy1 = -0.113 + 0.01*0
logy2 = -0.113 + 0.01*10
logy2-logy1

logy1 = -0.113 + 0.01*10
logy2 = -0.113 + 0.01*20
logy2-logy1

# Plot log(y) vs. x
ggplot(data, aes(x, log_y)) +
  geom_line() +
  ggtitle("Linear Relationship: log(y) vs. x") +
  xlab("x") +
  ylab("log(y)")

#Conclusion because the relationship is linear on the log scale a 10% increase is the same across the whole range of x

#Backtransforming to raw data scale
#Backtransform coefficient
exp(0.01)

#Working out a 10% increase when backtransforming is not as straight forward
log(y) = intercept + 0.01*multiple_paternity
y = exp(intercept + 0.01*multiple_paternity)

#Change in testes mass from 0 to 10%
y1 = exp(-0.113 + 0.01*0)
y2 = exp(-0.113 + 0.01*10)
y2-y1

#Change in testes mass from 10 to 20%
y1 = exp(0.01*50)
y2 = exp(0.01*60)
y2-y1


# Plot y vs. x
ggplot(data, aes(x, y)) +
  geom_line() +
  ggtitle("Exponential Relationship: y vs. x") +
  xlab("x") +
  ylab("y")

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



