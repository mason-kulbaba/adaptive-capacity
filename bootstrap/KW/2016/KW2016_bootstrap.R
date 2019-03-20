#### This code was taken from the G x E analsyis, and repurposed to generate CIs
#### for the second Chamaecrista adaptive capacity paper

#########################################################################
### 1. PRELIMINARIES: set working directory, load packages & data
#########################################################################


setwd("C:/Users/Mason Kulbaba/Dropbox/Rscripts/chamaecrista-adaptive-capacity/bootstrap/KW/2016")

# load packages
library(aster)

# load raster output
load("rout2016.RData")

# load redata file
load("redata2016.RData")

# load modmat.siredam matrix 
load("modmat.siredam2016.RData")

modmat.siredamKW<- as.matrix(modmat.siredam2016)

#set seed
set.seed(1729)

#########################################################################
### 2. PREPARE FOR BOOTSTRAPPING
#########################################################################

### set up graphical model
pred<- c(0,1,2,3,4) #make sure this matches your model

# designate family for each node in graphical model
fam<- c(1,1,2,2,2) #make sure this matches your model

# generate "hat" estimates from data
alpha.hat <- rout2016$alpha
sigma.hat <- rout2016$sigma
nu.hat <- rout2016$nu
b.hat <- rout2016$b
c.hat <- rout2016$c
sout <- summary(rout2016)
se.alpha.hat <- sout$alpha[ , "Std. Error"]
se.sigma.hat <- sout$sigma[ , "Std. Error"]
se.nu.hat <- sout$nu[ , "Std. Error"]

fixed <- rout2016$fixed
random <- rout2016$random
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
  phi.star <- as.numeric(as.vector(rout2016$obj$origin) +
                           modmat.tot %*% eff.star)
  theta.star <- astertransform(phi.star, rout2016$obj,
                               to.cond = "conditional", to.mean = "canonical")
  y.star <- raster(theta.star, pred, fam, rout2016$obj$root)
  y.star <- as.vector(y.star)
  
  rout.star <- reaster(y.star ~ varb + fit : (block), # make sure this matches your model
                       list(parental=~0 + fit:modmat.siredam2016),
                       pred, fam, varb, id, root, data = redata2016,
                       effects = c(alpha, c.star), sigma = sigma)
}

# generate function to compute VaW estimate
VaW<- function (rout2016){
  bhat<- rout2016$b
  bhat.sire<- bhat[grep("sireID", names(bhat))]
  hoom.star <- predict(rout2016$obj,  newcoef = rout2016$alpha)
  hoom.star<- matrix(hoom.star, ncol =5)
  hoom.star<- hoom.star[ , 5]  
  # mapping function
  map <- function(b) {
    stopifnot(length(b) == 1)
    stopifnot(is.finite(b))
    alpha <- rout2016$alpha
    alpha[10] <- alpha[10] + b # adding fixed effect for fit:block5 (e.g. alpha[10]), make sure this matches your model
    hoom.star <- predict(rout2016$obj, newcoef = alpha)
    hoom.star<- matrix(hoom.star, ncol = 5)
    return(hoom.star[2251, 5]) # return value of 5th node (seed #) for 2251st indiv (1st indiv in block 5). Change this to whatever your terminal fitness node is (e.g. fruit #)
  }
  fred<- Vectorize(map)
  bhat.sire.mu<- fred(bhat.sire)  
  hoom.star2<- predict(rout2016$obj, newcoef = rout2016$alpha, se.fit=TRUE, info.tol=1e-13)
  goom.star <- hoom.star2$gradient  
  moom.star<- goom.star[,5]
  moom.star<- matrix(moom.star, ncol=5)  
  # calcualtion for Va(w) 
  KW_Va<- 4*moom.star[2251 , 5]^2 * rout2016$nu[1]/map(0) # final calcuation of additive genetic variance for fitness
  
}

#########################################################################
### 3. DOUBLE BOOTSTRAPPING
#########################################################################

KW_Va.star <- rep(-1, nboot)
KW_Va.star.star<- rep(-1, ndoubleboot)
se.KW_Va.star <- rep(-1, nboot)
se.KW_Va.star_sd <- rep(-1, nboot)
for (iboot in 1:nboot) {
  
  tryCatch({
  rout.star <- boots(alpha.hat, sigma.hat)
  if("try-error" %in% class(rout.star)) stop("reaster error (Nelder-Mead with pickle fail?)")
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  alpha.star.star <- matrix(NaN, ndoubleboot, length(alpha.hat))
  sigma.star.star <- matrix(NaN, ndoubleboot, length(sigma.hat))
  nu.star.star <- matrix(NaN, ndoubleboot, length(nu.hat))
  
  tryCatch({
  KW_Va.star[iboot]<- VaW(rout.star) # function for VaW calcualtion at bootstrap level 1
  if("try-error" %in% class( KW_Va.star[iboot])) stop("reaster error (Nelder-Mead with pickle fail?)")
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  print(KW_Va.star[iboot]) # will print outer bootstrap estimate for each nboot interation
  print(iboot)
  
  for (idoubleboot in 1:ndoubleboot) {
    
    tryCatch({
    rout.star.star <- boots(rout.star$alpha, rout.star$sigma)
    if("try-error" %in% class(rout.star.star)) stop("reaster error (Nelder-Mead with pickle fail?)")
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    tryCatch({
    KW_Va.star.star[idoubleboot]<- VaW(rout.star.star) # same function for VaW calcualtion at bootstrap level 2
    if("try-error" %in% class( KW_Va.star.star[iboot])) stop("reaster error (Nelder-Mead with pickle fail?)")
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    
    alpha.star.star[idoubleboot, ] <- rout.star.star$alpha
    sigma.star.star[idoubleboot, ] <- rout.star.star$sigma
    nu.star.star[idoubleboot, ] <- rout.star.star$nu
    
    print(KW_Va.star.star[idoubleboot]) # will print each inner bootstrap estimate for each ndboubleboot iteration
    print(idoubleboot)
  }
  
  # This code uses the interquartile range as our measure of scale, when we set the inner bootstrap
  # to 31 (difference between the 25'th and 8'th elements).  We used IQR during the testing of this code
  # as it is a more robust measure of scale with small samples sizes (we developed this code with ndoubleboot = 7).
  # However, at larger sample sizes standard deviation is suitable. Therefore, the current code generates
  # both IQR and standard deviation scale metrics
  
  se.KW_Va.star[iboot] <- diff(sort(KW_Va.star.star)[c(8,25)]) # calculate IQR
  se.KW_Va.star_sd[iboot]<- sd(KW_Va.star.star) # calcualte standard deviation
}

#########################################################################
### 4. GENERATE BOOTSTRAP t CONFIDENCE INTERVALS
#########################################################################

# Save Rdata for entire bootstrap run
save.image("doubleboot_KW2016_run.RData")

#load combined resutls from Seema

kw2016<- read.csv("doubleboot_KW2016_all.csv")

KW_Va.star<- kw2016$KW_Va.star
se.KW_Va.star<- kw2016$se.KW_Va.star
se.KW_Va.star_sd<- kw2016$se.KW_Va.star_sd

# Generate Va(W) from data
KW_Va.hat <- VaW(rout2016) # generates Va(W) from data

# calculate t for IQR-based confidence interval
t <- (KW_Va.star - KW_Va.hat) / se.KW_Va.star

# set confidence level
conf.level <- 0.95

# generate critical values with t for IQR-based confidence interval
crit <- quantile(t, probs = c((1 - conf.level) / 2, (1 + conf.level) / 2))  # note this uses the t just calculated

# generate IQR and standard deviation measures of scale
foo <- KW_Va.star[1:nboot]    # the outer loop Va's
se.KW_Va.hat <- diff(sort(foo)[c(25,76)])  # middle 50 percent of 100 or so

se.KW_Va.hat_sd <- sd(foo)  # standard deviation calculation of scale

# Conf. Int. for KW_Va with IQR scale metric
KW_Va.hat - rev(crit) * se.KW_Va.hat

# calculate t for standard deviation-based confidence interval
t <- (KW_Va.star - KW_Va.hat) / se.KW_Va.star_sd

# generate critical values with t for standard deviation-based confidence interval
crit <- quantile(t, probs = c((1 - conf.level) / 2, (1 + conf.level) / 2))

# Conf. Int. for KW_Va with standard deviation scale metric
KW_Va.hat - rev(crit) * se.KW_Va.hat_sd # 0.702887, 5.140464
