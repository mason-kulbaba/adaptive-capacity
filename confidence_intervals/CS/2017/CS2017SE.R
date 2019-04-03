#### This code was taken from the G x E analsyis, and repurposed to generate CIs
#### for the second Chamaecrista adaptive capacity paper

#########################################################################
### 1. PRELIMINARIES: set working directory, load packages & data
#########################################################################


 setwd("C:/Users/Mason Kulbaba/Dropbox/Rscripts/chamaecrista-adaptive-capacity/bootstrap/CS/2017")

# load packages
library(aster)

# load raster output
load("rout2017.RData")

# load redata file
load("redata2017.RData")

# load modmat.siredam matrix 
load("modmat.siredam2017.RData")

modmat.siredamCS<- as.matrix(modmat.siredam2017)

#set random seed
set.seed(1729)

#########################################################################
### 2. PREPARE FOR BOOTSTRAPPING
#########################################################################

### set up graphical model
pred<- c(0,1,2,3,4) #make sure this matches your model

# designate family for each node in graphical model
fam<- c(1,1,2,2,2) #make sure this matches your model

# generate "hat" estimates from data
alpha.hat <- rout2017$alpha
sigma.hat <- rout2017$sigma
nu.hat <- rout2017$nu
b.hat <- rout2017$b
c.hat <- rout2017$c
sout <- summary(rout2017)
se.alpha.hat <- sout$alpha[ , "Std. Error"]
se.sigma.hat <- sout$sigma[ , "Std. Error"]
se.nu.hat <- sout$nu[ , "Std. Error"]

fixed <- rout2017$fixed
random <- rout2017$random
modmat.tot <- cbind(fixed, Reduce(cbind, random))

nfix <- ncol(fixed)
nrand <- sapply(random, ncol)
a.hat <- rep(sigma.hat, times = nrand)

### prepare for boots function that performs a single bootstrap

#alpha.star <- matrix(NaN, nboot, length(alpha.hat))
#sigma.star <- matrix(NaN, nboot, length(sigma.hat))
#nu.star <- matrix(NaN, nboot, length(nu.hat))
#se.alpha.star <- alpha.star
#se.sigma.star <- sigma.star
#se.nu.star <- nu.star

# generate function to compute VaW estimate
VaW <- function (rout2017){
  bhat<- rout2017$b
  bhat.sire<- bhat[grep("sireID", names(bhat))]
  hoom.star <- predict(rout2017$obj,  newcoef = rout2017$alpha)
  hoom.star<- matrix(hoom.star, ncol =5)
  hoom.star<- hoom.star[ , 5]  
  # mapping function
  map <- function(b) {
    stopifnot(length(b) == 1)
    stopifnot(is.finite(b))
    alpha <- rout2017$alpha
    alpha[8] <- alpha[8] + b # adding fixed effect for fit:block3 (e.g. alpha[8]), make sure this matches your model
    hoom.star <- predict(rout2017$obj, newcoef = alpha)
    hoom.star<- matrix(hoom.star, ncol = 5)
    return(hoom.star[971, 5]) # return value of 5th node (seed #) for 971st indiv (1st indiv in block 3). Change this to whatever your terminal fitness node is (e.g. fruit #)
  }
  fred<- Vectorize(map)
  bhat.sire.mu<- fred(bhat.sire)  
  hoom.star2<- predict(rout2017$obj, newcoef = rout2017$alpha, se.fit=TRUE, info.tol = 1e-13)
  goom.star <- hoom.star2$gradient  
  moom.star<- goom.star[,5]
  moom.star<- matrix(moom.star, ncol=5)  
  # calcualtion for Va(w) 
  est <- 4*moom.star[971 , 5]^2 * rout2017$nu[1]/map(0) # final calcuation of additive genetic variance for fitness
  # SE
  se <- 4*moom.star[971,5]^2 * sout$nu["parental", "Std. Error"]/map(0)
  ci <- c(est - 1.96*se, est, est + 1.96*se)
  return(ci)
}

#########################################################################
### 3. Get SE
#########################################################################

CI <- VaW(rout2017) # generates Va(W) from data
print(CI)

