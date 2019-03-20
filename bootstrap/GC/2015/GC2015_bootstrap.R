#### This code was taken from the G x E analsyis, and repurposed to generate CIs
#### for the second Chamaecrista adaptive capacity paper

#########################################################################
### 1. PRELIMINARIES: set working directory, load packages & data
#########################################################################


setwd("C:/Users/Mason Kulbaba/Dropbox/Rscripts/chamaecrista-adaptive-capacity/bootstrap/GC/2015")

# load packages
library(aster)


# load raster output
load("rout2015.RData")

# load redata file
load("redata2015.RData")

# load modmat.siredam matrix 
load("modmat.siredam2015.RData")

modmat.siredamGC<- as.matrix(modmat.siredam2015)

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
alpha.hat <- rout2015$alpha
sigma.hat <- rout2015$sigma
nu.hat <- rout2015$nu
b.hat <- rout2015$b
c.hat <- rout2015$c
sout <- summary(rout2015)
se.alpha.hat <- sout$alpha[ , "Std. Error"]
se.sigma.hat <- sout$sigma[ , "Std. Error"]
se.nu.hat <- sout$nu[ , "Std. Error"]

fixed <- rout2015$fixed
random <- rout2015$random
modmat.tot <- cbind(fixed, Reduce(cbind, random))

nfix <- ncol(fixed)
nrand <- sapply(random, ncol)
a.hat <- rep(sigma.hat, times = nrand)

# set outer (nboot) and inner (ndoubleboot)
nboot <- 100
ndoubleboot <- 31

### prepare for boots function that performs a single bootstrap

alpha.star <- matrix(NaN, nboot, length(alpha.hat))
sigma.star <- matrix(NaN, nboot, length(sigma.hat))
nu.star <- matrix(NaN, nboot, length(nu.hat))
se.alpha.star <- alpha.star
se.sigma.star <- sigma.star
se.nu.star <- nu.star

# make boots function
boots <- function(alpha, sigma) {
  
  a.hat <- rep(sigma, times = nrand)
  
  c.star <- rnorm(sum(nrand))
  b.star <- a.hat * c.star
  eff.star <- c(alpha, b.star)
  phi.star <- as.numeric(as.vector(rout2015$obj$origin) +
                           modmat.tot %*% eff.star)
  theta.star <- astertransform(phi.star, rout2015$obj,
                               to.cond = "conditional", to.mean = "canonical")
  y.star <- raster(theta.star, pred, fam, rout2015$obj$root)
  y.star <- as.vector(y.star)
  
  rout.star <- reaster(y.star ~ varb + fit : (block), # make sure this matches your model
                       list(parental=~0 + fit:modmat.siredam2015),
                       pred, fam, varb, id, root, data = redata2015,
                       effects = c(alpha, c.star), sigma = sigma)
}

# generate function to compute VaW estimate
VaW<- function (rout2015){
  bhat<- rout2015$b
  bhat.sire<- bhat[grep("sireID", names(bhat))]
  hoom.star <- predict(rout2015$obj,  newcoef = rout2015$alpha)
  hoom.star<- matrix(hoom.star, ncol =5)
  hoom.star<- hoom.star[ , 5]  
  # mapping function
  map <- function(b) {
    stopifnot(length(b) == 1)
    stopifnot(is.finite(b))
    alpha <- rout2015$alpha
    alpha[11] <- alpha[11] + b # adding fixed effect for fit:block6 (e.g. alpha[11]), make sure this matches your model
    hoom.star <- predict(rout2015$obj, newcoef = alpha)
    hoom.star<- matrix(hoom.star, ncol = 5)
    return(hoom.star[1875, 5]) # return value of 5th node (seed #) for 1875th indiv (1st indiv in block 6). Change this to whatever your terminal fitness node is (e.g. fruit #)
  }
  fred<- Vectorize(map)
  bhat.sire.mu<- fred(bhat.sire)  
  hoom.star2<- predict(rout2015$obj, newcoef = rout2015$alpha, se.fit=TRUE)
  goom.star <- hoom.star2$gradient  
  moom.star<- goom.star[,5]
  moom.star<- matrix(moom.star, ncol=5)  
  # calcualtion for Va(w) 
  GC_Va<- 4*moom.star[1875 , 5]^2 * rout2015$nu[1]/map(0) # final calcuation of additive genetic variance for fitness
  
}

#########################################################################
### 3. DOUBLE BOOTSTRAPPING
#########################################################################

GC_Va.star <- rep(-1, nboot)
GC_Va.star.star<- rep(-1, ndoubleboot)
se.GC_Va.star <- rep(-1, nboot)
se.GC_Va.star_sd <- rep(-1, nboot)
for (iboot in 1:nboot) {
  
  tryCatch({
  rout.star <- boots(alpha.hat, sigma.hat)
  if("try-error" %in% class(rout.star)) stop("reaster error (Nelder-Mead with pickle fail?)")
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  alpha.star.star <- matrix(NaN, ndoubleboot, length(alpha.hat))
  sigma.star.star <- matrix(NaN, ndoubleboot, length(sigma.hat))
  nu.star.star <- matrix(NaN, ndoubleboot, length(nu.hat))
  
  tryCatch({
  GC_Va.star[iboot]<- VaW(rout.star) # function for VaW calcualtion at bootstrap level 1
  if("try-error" %in% class( GC_Va.star[iboot])) stop("aster.predict error (direction of recession error?)")
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  print(GC_Va.star[iboot]) # will print outer bootstrap estimate for each nboot interation
  print(iboot)
  
  for (idoubleboot in 1:ndoubleboot) {
    
    tryCatch({
    rout.star.star <- boots(rout.star$alpha, rout.star$sigma)
    if("try-error" %in% class(rout.star.star)) stop("reaster error (Nelder-Mead with pickle fail?)")
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    tryCatch({
    GC_Va.star.star[idoubleboot]<- VaW(rout.star.star) # same function for VaW calcualtion at bootstrap level 2
    if("try-error" %in% class( GC_Va.star.star[iboot])) stop("aster.predict error (direction of recession error?)")
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    alpha.star.star[idoubleboot, ] <- rout.star.star$alpha
    sigma.star.star[idoubleboot, ] <- rout.star.star$sigma
    nu.star.star[idoubleboot, ] <- rout.star.star$nu
    
    print(GC_Va.star.star[idoubleboot]) # will print each inner bootstrap estimate for each ndboubleboot iteration
    print(idoubleboot)
  }
  
  # This code uses the interquartile range as our measure of scale, when we set the inner bootstrap
  # to 31 (difference between the 25'th and 8'th elements).  We used IQR during the testing of this code
  # as it is a more robust measure of scale with small samples sizes (we developed this code with ndoubleboot = 7).
  # However, at larger sample sizes standard deviation is suitable. Therefore, the current code generates
  # both IQR and standard deviation scale metrics
  
  se.GC_Va.star[iboot] <- diff(sort(GC_Va.star.star)[c(8,25)]) # calculate IQR
  se.GC_Va.star_sd[iboot]<- sd(GC_Va.star.star) # calcualte standard deviation
}

#########################################################################
### 4. GENERATE BOOTSTRAP t CONFIDENCE INTERVALS
#########################################################################

# Save Rdata for entire bootstrap run
save.image("doubleboot_GC2015_run.RData")

# Generate Va(W) from data
GC_Va.hat <- VaW(rout2015) # generates Va(W) from data

# calculate t for IQR-based confidence interval
t <- (GC_Va.star - GC_Va.hat) / se.GC_Va.star

# set confidence level
conf.level <- 0.95

# generate critical values with t for IQR-based confidence interval
crit <- quantile(t, probs = c((1 - conf.level) / 2, (1 + conf.level) / 2))  # note this uses the t just calculated

# generate IQR and standard deviation measures of scale
foo <- GC_Va.star[1:nboot]    # the outer loop Va's
se.GC_Va.hat <- diff(sort(foo)[c(25,76)])  # middle 50 percent of 100 or so

se.GC_Va.hat_sd <- sd(foo)  # standard deviation calculation of scale

# Conf. Int. for GC_Va with IQR scale metric
GC_Va.hat - rev(crit) * se.GC_Va.hat

# calculate t for standard deviation-based confidence interval
t <- (GC_Va.star - GC_Va.hat) / se.GC_Va.star_sd

# generate critical values with t for standard deviation-based confidence interval
crit <- quantile(t, probs = c((1 - conf.level) / 2, (1 + conf.level) / 2))

# Conf. Int. for GC_Va with standard deviation scale metric
GC_Va.hat - rev(crit) * se.GC_Va.hat_sd