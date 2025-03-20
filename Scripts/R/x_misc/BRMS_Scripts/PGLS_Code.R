## basic phylogenetic comparative analysis 
## module 10
## 12/01/2022

## clean environment & plots
rm(list=ls()) 
graphics.off()

## libraries
library(visreg)
library(ggplot2)
library(ape)
library(plyr)
library(ggtree)
library(caper)
library(phytools)
library(MuMIn)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree", force = TRUE)

library(ggtree)

## load mammal dataset
setwd("~/OneDrive - University of Oklahoma/R statistics workshop Fall 2022/Module 10")
data=read.csv("mammals_latin.csv",header=T)

## do larger species sleep less?
data=data[!is.na(data$TotalSleep),]
data=data[!is.na(data$BodyWt),]

## visualize
ggplot(data,aes(BodyWt,TotalSleep))+
  geom_point()+
  scale_x_log10()+
  geom_smooth(method="lm")

## check linear model
lmod=lm(TotalSleep~log10(BodyWt),data=data)
summary(lmod)

## load in a mammal-wide phylogeny
tree=read.nexus("mammal phylogeny.tre")

## the phylogeny tips contain extra information
tree$tip.label[1:10]

## let's simplify
tree$tip.label=sapply(strsplit(tree$tip.label,"_"),function(x) paste(x[1],x[2]))

## do all the tips match our data?
setdiff(data$Latin,tree$tip.label)

## let's fix mismatch names
data$tip.label=revalue(data$Latin,
                       c("Spermophilus parryii"="Urocitellus parryii",
                         "Pan troglodytes troglodytes"="Pan troglodytes",
                         "Equus asinus"="Equus africanus",
                         "Myotis lucifugus lucifugus"="Myotis lucifugus",
                         "Sus scrofa domesticus"="Sus scrofa"))

## check if that worked
setdiff(data$tip.label,tree$tip.label) 

## use keep.tip to simplify the tree to our species
tree=keep.tip(tree,data$tip.label)

## base plot our tree
par(oma=c(0,0,0,0),mar=c(0,0,0,0))
plot(tree,cex=0.75)

## ggplot
ggtree(tree)+
  geom_tiplab(size=2.25,offset=2)+
  xlim(0,250)

## we can use the caper package to combine a tree with meta-data
## save our data in the same order as the phylogey

## save species names
data$label=data$tip.label

## merge data and tree
cdata=comparative.data(phy=tree,data=data,names.col=label,vcv=T,na.omit=F,warn.dropped=T)

## pagel's lambda in totalsleep
pl=pgls(TotalSleep~1,data=cdata,lambda="ML")
summary(pl)
summary(pl)$param["lambda"]

## phylogenetic generalized least squares: phylogenetically controlled regression
mod=pgls(TotalSleep~log10(BodyWt),data=cdata,lambda="ML")
summary(mod) ## still significant, but much weaker effect

## let's compare model coefficients with a LM
coef(lmod)
coef(mod)

## compare R2
summary(lmod)$adj.r.squared
summary(mod)$r.squared

library(tidytree)

## combine tree and data with treeio
cdata$data$label=cdata$data$tip.label
dtree=treeio::full_join(as.treedata(cdata$phy),cdata$data,by="label")
dtree

## ggtree
ggtree(dtree)+
  xlim(0,250)+
  geom_tiplab(aes(label=Species),
              size=2.5,offset=2)+
  geom_tippoint(aes(size=TotalSleep,fill=log10(BodyWt)),shape=21)+
  scale_size_continuous(range=c(0.2,3))+
  scale_fill_viridis_c()+
  theme(legend.position="bottom")

## do species spend over 25% of their sleep dreaming?
cdata$data$dtime=cdata$data$Dreaming/cdata$data$TotalSleep
cdata$data$dlot=ifelse(cdata$data$dtime>=0.25,1,0)
cdata$data$dlot=factor(cdata$data$dlot)

## remove NA
cdata=cdata[!is.na(cdata$data$dlot),]

## merge back with phylo
dtree=treeio::full_join(as.treedata(cdata$phy),cdata$data,by="label")

## plot it
ggtree(dtree)+
  xlim(0,250)+
  geom_tiplab(aes(label=Species),
              size=2.5,offset=2)+
  geom_tippoint(aes(fill=dlot),shape=21)+
  theme(legend.position="bottom")

## ancestral trait reconstruction: we need the trait states
state=as.character(cdata$data$dlot)
names(state)=rownames(cdata$data)

## estimate ancestral states under an equal-rate model and all-rates-different model
fitER=ace(state,cdata$phy,model="ER",type="discrete",method="ML",CI=T) 
fitARD=ace(state,cdata$phy,model="ARD",type="discrete",method="ML",CI=T) 

## use AIC to compare models
AIC(fitER)
AIC(fitARD)
## equal rates model fits slightly better

## stochastic character map
set.seed(1)
ssm=make.simmap(cdata$phy,state,model="ER",nsim=100)

## summarize
summary(ssm,plot=T)

## tasks!
## (1) estimate phylogenetic signal in predation risk
## (2) use a PGLS to ask if predation risk is correlated with total sleep
## (3) make a new binary trait for low/high predation risk, reconstruct its evoluion
