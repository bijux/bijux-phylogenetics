#***********************************************************************************
#BIOR90: Phylogenetic Comparative Analyses ----
#Trait evolution by Charlie Cornwallis
#PCM1 Script to accompany lecture
#***********************************************************************************

#***********************************************************************************
#Packages ----
#***********************************************************************************
#pacman::p_load makes loading packages easier
pacman::p_load(openxlsx,tidyverse,RColorBrewer,ape,MCMCglmm,picante,geiger,phytools,ggtree, treeio,ggimage)

#If packages don't install, google for solutions. 
#ggtree often has this problem that can usually be solved with
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree",force = TRUE)

#If you need to citing R and R packages in your work you can check how to do this
citation() #R
citation("ape") #Specific package

#***********************************************************************************
#1. Data importing, checking and manipulation ----
#***********************************************************************************
#If you started R studio using the Rproject file then the correct working directory should be set automatically and work on any computer

#Data on Primates - read in from excel
primate_raw <- readWorkbook("Data/primate_raw.xlsx")

#***************************
#Checking how the data looks and variable classifications
#***************************
#What are the column names?
names(primate_raw)

#How does the dataset look?
View(primate_raw)

#How are variables classified?
str(primate_raw) #Oh dear many are incorrectly classified

#Solution: - Reclassify manually e.g.
primate_raw$sex_dimorphism<-as.factor(primate_raw$sex_dimorphism)

#Or can change several columns at once
primate_raw <- primate_raw %>% mutate(across(body_mass:social_group_size, ~as.numeric(.)))

str(primate_raw) #Has it worked?

#***************************
#Check for missing data (NAs) and multiple measurements per species?
#***************************
#How many NAs are there per variable?
primate_raw %>% summarise(across(order:mating_system, ~sum(is.na(.))))

#Show rows with missing data (not (!) complete cases)
primate_raw %>% filter(!complete.cases(.))

#Are their multiple measurements for species?
primate_raw %>% group_by(species) %>% 
                summarise(n=n())  %>% 
                dplyr::filter(n>1)

#Lets check out this species with multiple measurements and missing values
primate_raw %>% filter(species == "Galago_alleni")

#Solution - summarise data and remove NAs at the same time
primate <- primate_raw %>% group_by(family, species) %>% 
                           summarise(across(body_mass:social_group_size, 
                                      .fns = ~mean(.x,na.rm=T), #fns stands for function which in this case is the "mean" and "na.rm=T" removes missing values
                                      .names = "{col}"),#This could be used to rename columns if you wanted
                                      sex_dimorphism=sex_dimorphism[!is.na(sex_dimorphism)], #Sexual dimorphism stays the same but removes NAs (!is.na())
                                      mating_system=mating_system[!is.na(mating_system)]) #Mating system stays the same but removes NAs

primate %>% group_by(species) %>% 
  summarise(n=n())  %>% 
  filter(n>1)

#***********************************************************************************
#2. Tree importing, checking and manipulation ----
#***********************************************************************************
#Read in Tree using Package ape, read.tree multipurpose function, for nexus, newick and most tree formats
primatetree<-read.tree("Data/primatetree.nex")

#See what the file contains
primatetree

#You can access these with the variable operator $
primatetree$tip.label

#Is the tree rooted?
is.rooted(primatetree)

#Does it have polytomies
is.binary(primatetree)

#Is the tree ultrametric? Ultrametric means that all tips have same distance to the root (genetic distances between all species are relative rather than absolute), which is important for some analyses
is.ultrametric(primatetree)

#Does it have node labels?
primatetree$node.label

#Lets label the nodes so we can keep track of ancestors
primatetree <- makeNodeLabel(primatetree)

#Check correspondence between data and tree

#Are all species in the dataset in the tree?
table(primate$species %in% primatetree$tip.label)

#Are all tree tips in the dataset?
table(primatetree$tip.label %in% primate$species)

#Which are missing?
missingtips<-primatetree$tip.label[primatetree$tip.label %in% primate$species == "FALSE"]
missingtips

#Remove these species from the tree using drop.tip
primatetree<-ape::drop.tip(primatetree,missingtips)

#***********************************************************************************
#* Saving processed files ----
#***********************************************************************************
#After spending time sorting out the files you may want to save the processed files for using / importing later on for analyses

#We can easily write & read data and trees to specific files, for example:
write.csv(primate,"Data/primate.csv")
write.tree(primatetree,"Data/trimmed_primatetree.nex")

#Or you can save packages of files in one RData object
save(primate,primatetree, file = "Data/primate.RData") # I prefer this option as it creates less files, less coding and is simpler to keep track of

save.image("Data/primate.RData") #Or all files present in your environment. Simple but trash can accummulate.

#***********************************************************************************
#3. Plotting ----
#***********************************************************************************
#* Package ape - quite straight forward and good for checking things out. Becomes a bit limited if you want to do more complicated figures ----
#************************************
par(mar=c(0,0,0,0))
plot(primatetree)

#Change size of tip labels
plot(primatetree,cex=0.5)

#Add node numbers
nodelabels(cex=0.5) 

#Add tip numbers
tiplabels(cex=0.5) 

#Offset labels
plot(primatetree,cex=0.5,label.offset = 0.1)
tiplabels(cex=0.5) 

#Can also plot in different ways
plot(primatetree,cex=0.5,type="cladogram",main="cladogram")
plot(unroot(primatetree),cex=0.5,type="unrooted",main="unrooted")
plot(primatetree,cex=0.5,type="fan",main="fan")

#************************************
#* Package Phytools - another package with some useful features for exploring trees ----
#************************************
plotTree(primatetree,fsize=0.5)
nodelabels()

#now extract the clade descended from node #77
tt77<-extract.clade(primatetree,77)
plotTree(tt77,fsize=0.5)

#rotate the tree around node #130
rt.130<-rotateNodes(primatetree,130)
plotTree(rt.130,fsize=0.5)

#Compare trees
par(mfrow=c(1,2))
plotTree(primatetree,fsize=0.5)
nodelabels()
plotTree(rt.130,fsize=0.5)
nodelabels()

# now rotate all nodes
rt.all<-rotateNodes(primatetree,"all")
plotTree(rt.all,fsize=0.5)
plotTree(primatetree,fsize=0.5)

# check if tree & rt.all are equal
all.equal(primatetree,rt.130)

# check if tree & tt77 are equal
all.equal(primatetree,tt77)

#************************************
#* Package ggtree - alot of useful plotting functions, handy that it is similar to ggplot, but can be frustrating and a bit buggy ----
#************************************
ggtree(primatetree)
ggtree(primatetree,layout="circular")
ggtree(primatetree, branch.length="none") #You can easily turn your tree into a cladogram with the branch.length = “none” parameter.

#add tiplabels
ggtree(primatetree) + geom_tiplab()
ggtree(primatetree,layout="circular") + geom_tiplab()

#add scale
ggtree(primatetree) + geom_treescale(offset=-2,linesize=1)

#add nodelabels
ggtree(primatetree) + geom_text(aes(label=node),size=3)

#add titles
ggtree(primatetree) + ggtitle("A great title")

#annotate selected clades in trees.
ggtree(primatetree) + xlim(NA, 2)+ geom_cladelabel(node=130, label="test label", align=T, color='red') +
geom_cladelabel(node=78, label="another clade", align=T, color='blue') 

#Can add layers similar to ggplot
p <- ggtree(primatetree) + geom_text(aes(label=node),size=3)
p + geom_cladelabel(node=130, label="test label", align=T, color='red') +geom_cladelabel(node=78, label="another clade", align=T, color='blue')

#************************************
#* Combine data and trees ----
#************************************
#** ape ----
#Lets do this with ape - need to make sure data is in the same order as tree
#Make a dataset with order of tips
tip_order <- data.frame(tip=primatetree$tip.label,
                        order=1:length(primatetree$tip.label))
head(tip_order)

#Add this to the primate dataframe
primate$tip_order <- tip_order$order[match(primate$species,tip_order$tip)]

#Now you can sort data frame according to tip order
primate <- primate %>% arrange(tip_order)

#Has it worked?
head(primate)
primatetree$tip.label[1:6]

#Now we can plot the data on the tree
par(mar=c(0,0,0,0))
plot(primatetree, show.tip.label=T, label.offset=0.05,cex=0.6)
tiplabels(pch=21, bg="grey38",cex=primate$longevity/400) # cex = size in plot function. need to /300 otherwise symbols are too big


#** ggtree ----
#Lets do this with ggtree and the treeio package
primate$node<-nodeid(primatetree,primate$species)

#Lets see what this has done - in one step it has done what we did above
View(primate)

#can now use this variable to combine tree and data into one object
p_tree_data <- treeio::full_join(primatetree, primate, by = 'node')

#Create an S4 object = different slots for different types of data
p_tree_data

#You can still access the data and tree e.g.
p_tree_data@phylo$tip.label
View(p_tree_data@extraInfo)

#Now you have this combined object you can do more complicated plots
ggtree(p_tree_data) +
  geom_tiplab(aes(colour="red"),offset=0.1) +
  geom_tippoint(aes(x=x+0.,size=longevity),fill="grey",colour="black",shape=21)

#Plot area needs adjusting and perhaps some other adjustments could make it nicer
ggtree(p_tree_data) +
  geom_tiplab(colour="red",offset=0.1,size=2) +
  geom_tippoint(aes(x=x+0.05,size=longevity),fill="grey",colour="black",shape=21)+
  xlim(0,1.5)+
  ylim(0,75)+
  guides(size=guide_legend(title="Longevity (months)"),color="none")

#Much more on ggtree can be found here: https://guangchuangyu.github.io/software/ggtree/documentation/
#and here: http://www.randigriffin.com/2017/05/11/primate-phylogeny-ggtree.html

....Back to the slides

#***********************************************************************************
#4. Phylogenetic signal in continuous traits ----
#***********************************************************************************
#Lets simulate some data to see how traits look when there values are are product of evolutionary history
random_tree<-rcoal(30)
random_data<-rTraitCont(random_tree, model="BM", sigma = 0.5,root.value=1)

plot(random_tree, show.tip.label=F)
tiplabels(pch=21, cex=random_data*2, bg="orange")

#Have a play around with the root value and sigma (standarde deviation in trait)
random_data2<-rTraitCont(random_tree, model="BM", sigma = 0.5,root.value=10)
random_data3<-rTraitCont(random_tree, model="BM", sigma = 5,root.value=1)

par(mfrow=c(1,3))
plot(random_tree, show.tip.label=F)
tiplabels(pch=21, cex=random_data, bg="orange")

plot(random_tree, show.tip.label=F)
tiplabels(pch=21, cex=random_data2, bg="orange")

plot(random_tree, show.tip.label=F)
tiplabels(pch=21, cex=random_data3, bg="orange")

#Different traits may evolve in different ways aka 'the mode of evolution may be different' 
par(mfrow=c(1,2))
random_data4<-rTraitCont(random_tree, model="OU", sigma = 5,root.value=1,alpha=0)#alpha is the strength of attraction to an optimum. 0 = no attractionso should be the same as the BM model.

plot(random_tree, show.tip.label=F)
tiplabels(pch=21, cex=random_data4, bg="orange")

#Lets increase the alpha and see what happens
random_data5<-rTraitCont(random_tree, model="OU", sigma = 5,root.value=1,alpha=5)

plot(random_tree, show.tip.label=F)
tiplabels(pch=21, cex=random_data5, bg="orange")

#*******************************************
#* Quantifying phylogenetic signal: Pagels Lamda using fitContinuous in Geiger ----
#*******************************************

Lambda_random_data<-fitContinuous(phy=random_tree, dat = random_data,model = "lambda")
Lambda_random_data

#Try with random_data2 and random_data3. What do you find?

#What happens if you do it for random_data4 and random_data5?
Lambda_random_data5<-fitContinuous(phy=random_tree, dat = random_data5,model = "lambda",control = list(niter = 1000))
Lambda_random_data5

#*******************************************
#* Example data: Primate longevity ----
#*******************************************
#How is the data distributed?
hist(primate$longevity)

#Format data for Geiger tests and see how it looks on the phylogeny?
primateLL<-as.numeric(primate$longevity[match(primatetree$tip.label,primate$species)])
names(primateLL)<-primatetree$tip.label

par(mfrow=c(1,1))
plot(primatetree, show.tip.label=T, label.offset=0.05,cex=0.6)
tiplabels(pch=21, bg="grey38",cex=primateLL/300)

#*******************************************
#* Pagels Lamda using fitContinuous in Geiger ----
#*******************************************

#Estimate Lambda
Lambda_LL<-fitContinuous(phy= primatetree, dat = primateLL,model = "lambda")
Lambda_LL
Lambda_LL$opt$lnL


#*******************************************
#* Testing if Pagels Lamda  is significantly different from 0 (start phylogeny) ----
#*******************************************
##First transform the tree to Lambda = 0 
primatetreeL0<-rescale(primatetree, "lambda", 0)

##See how they look
par(mfrow=c(1,2),mar=c(0,0,2,0))
plot(primatetree,cex=0.4,label.offset=0.1)
title("Real Phylogeny", cex=28)
tiplabels(pch=21, bg="grey38",cex=primateLL/300)

plot(primatetreeL0, cex=0.4,label.offset=0.1)
title("Lambda = 0", cex=28)
tiplabels(pch=21, bg="grey38",cex=primateLL/300)


#You can see how this influences the covariance-variance matrix between species
vcv.phylo(primatetreeL0)[1:3,1:3]
vcv.phylo(primatetree)[1:3,1:3]

#Calculate the loglikelihood when Lambda is set to 0
Lambda_LL0<-fitContinuous(phy= primatetreeL0, dat = primateLL,model = "lambda")
Lambda_LL0$opt$lnL

#Use a log-likelihood ratio test to calculate if Lambda is different from 0
##Calculate the  log likelihood ratio
LLDiff0 <- -2*(Lambda_LL0$opt$lnL - Lambda_LL$opt$lnL)
LLDiff0

##Perform Chi sq test
pchisq(LLDiff0, df=1,lower.tail = FALSE)


#***********************************************************************************
#5. Ancestral state reconstructions ----
#***********************************************************************************
#* Continuous traits ----
long_ace<-ace(primate$longevity,primatetree,type="continuous",method="pic")

#Lets look at the output
names(long_ace$ace) #nodes names
long_ace$ace #ancestral estimates
long_ace$CI95 #)95% CIs

#We can also plot it
#Put the information from the ace's together with tip data in a data frame 
primate_ace_long<-data.frame(species=c(primate$species, #species names
                                       names(long_ace$ace)), #nodes names
                             longevity=c(primate$longevity, #longevity values
                                         long_ace$ace)) #estimated ancestral longevity

#Combine data and tree
primate_ace_long$node<-nodeid(primatetree,primate_ace_long$species)#Convert tips labels to numbers
View(primate_ace_long)#how does it look?
primate_ace_long_tree <- full_join(primatetree,primate_ace_long, by = 'node')#Combine with tree in ggtree object

ggtree(primate_ace_long_tree)+
  geom_tree(aes(color=primate_ace_long_tree@extraInfo$longevity),size=1) + # changes colour of branches to indicate ancestral estimates
  geom_tiplab(aes(x=x+0.08),size=2) + #plot tiplabels
  xlim(0,1.6)+
  guides(color=guide_legend(title="Lifespan"),size=guide_legend(title="Lifespan"))+
  theme(legend.position = c(.1, .85))

#Lets examine the estimated lifespans of Bonobos and Gibbons and their the most recent common ancestor
primate_ace_long %>% filter(species == "Pan_paniscus")
primate_ace_long %>% filter(species == "Hylobates_lar")

#Node number of most recent common ancestor
getMRCA(primatetree, tip=c("Pan_paniscus", "Hylobates_lar"))

#Estimate of node
primate_ace_long %>% filter(node == "107")

#***********************************************************************************
#* How many times has lifespan increased across primates? ----
#***********************************************************************************

#Create a dataset of ancestors and descendants
tree_df<-tidytree::as_tibble(primatetree)

#Add in ancestral estimates matching by node id's
tree_df <- tree_df %>% mutate(ancestor_long=primate_ace_long$longevity[match(parent,primate_ace_long$node)],
                              descendant_long=primate_ace_long$longevity[match(node,primate_ace_long$node)])

tree_df

#How many increases are there?
tree_df <- tree_df %>% mutate(diff_long = descendant_long - ancestor_long) #Calculate the difference between ancestor and descendant
tree_df %>% filter(diff_long > 0) %>% summarise(n=n()) #Count the number of increases

tree_df %>% filter(diff_long > 12) %>% summarise(n=n()) #Count the number of increases

#You can save all your work in a RData object
save.image("./Results/primate_results.RData")

#***********************************************************************************
#END ####
#***********************************************************************************

