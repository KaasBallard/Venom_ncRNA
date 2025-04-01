library(ape)
library(brms)

######################################
#Island size PGLMM

tree <- read.nexus("PHL_Reptile_Tree_IncludingTurtles_CorrectOrder")

data <- read.csv("2024_02_26_PGLMM_Shannon_DendroLuzon.csv")

#Check that the names match

setdiff(data$tip,tree$tip.label)

#No difference in the tip labels in the metadata and the tip labels in the tree...good to go!

#Calculate correlation matrix
cmatrix=vcv.phylo(tree,cor=T)

#make phylo variable (tip lables)
data$tree=data$tip

#Shannon Island Size regression
chain=4
iter=20000
thin=50
warmup=7000

mod=brm(shannon~IslandSizeSqKm+(1|tip)+(1|gr(tree,cov=A)),data=data,family=gaussian,data2=list(A=cmatrix), 
        chain=chain,iter=iter,thin=thin,warmup=warmup, control = list(adapt_delta = 0.85))

mcmc_plot(mod,type="trace")

## check Rhat and summary values
summary(mod)

plot(mod, N = 2, ask = FALSE)

## pulling the posterior means
fixef(mod)

## easy plotting
plot(conditional_effects(mod), points = TRUE)

## Running on log-transformed island size
tree <- read.nexus("PHL_Reptile_Tree_IncludingTurtles_CorrectOrder")

data <- read.csv("2024_02_26_PGLMM_Shannon_DendroLuzon.csv")

data$logsize = log(data$IslandSizeSqKm)
head(data)

#Check that the names match

setdiff(data$tip,tree$tip.label)

#No difference in the tip labels in the metadata and the tip labels in the tree...good to go!

#Calculate correlation matrix
cmatrix=vcv.phylo(tree,cor=T)

#make phylo variable (tip lables)
data$tree=data$tip

chain=4
iter=20000
thin=50
warmup=7000

mod=brm(shannon~logsize+(1|tip)+(1|gr(tree,cov=A)),data=data,family=gaussian,data2=list(A=cmatrix), 
        chain=chain,iter=iter,thin=thin,warmup=warmup, control = list(adapt_delta = 0.95))

## check Rhat and summary values
summary(mod)

plot(mod, N = 2, ask = FALSE)

## pulling the posterior means
fixef(mod)

# and confidence interval
library(bayestestR)
hdi(mod)

## easy plotting
plot(conditional_effects(mod), points = TRUE)








