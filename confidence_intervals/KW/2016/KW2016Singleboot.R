setwd()
# load packages
library(aster)

# load raster output
load("rout2016.RData")

# load redata file
load("redata2016.RData")

# load modmat.siredam matrix 
load("modmat.siredam2016.RData")

modmat.siredamKW<- as.matrix(modmat.siredam2016)

#set random seed
set.seed(51729)


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

# set outer (nboot) only
nboot <- 100

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
    alpha[10] <- alpha[10] + b # adding fixed effect for fit:block6 (e.g. alpha[11]), make sure this matches your model
    hoom.star <- predict(rout2016$obj, newcoef = alpha)
    hoom.star<- matrix(hoom.star, ncol = 5)
    return(hoom.star[2251, 5]) # return value of 5th node (seed #) for 2274th indiv (1st indiv in block 6). Change this to whatever your terminal fitness node is (e.g. fruit #)
  }
  fred<- Vectorize(map)
  bhat.sire.mu<- fred(bhat.sire)  
  hoom.star2<- predict(rout2016$obj, newcoef = rout2016$alpha, se.fit=TRUE, info.tol = 1e-13)
  goom.star <- hoom.star2$gradient  
  moom.star<- goom.star[,5]
  moom.star<- matrix(moom.star, ncol=5)  
  # calcualtion for Va(w) 
  KW_Va<- 4*moom.star[2251 , 5]^2 * rout2016$nu[1]/map(0) # final calcuation of additive genetic variance for fitness
  soutstar <- summary(rout2016)
  KW_SE<- 4 * moom.star[2251, 5]^2 * soutstar$nu["parental", "Std. Error"] / map(0) 
  vaandse <- c(KW_Va,KW_SE)
  return(vaandse)
}

# SINGLE BOOTSTRAPPING

KW_Va.star <- rep(-1, nboot)
se.KW_Va.star <- rep(-1, nboot)
se.KW_Va.star_sd <- rep(-1, nboot)
for (iboot in 1:nboot) {
  
  tryCatch({
    rout.star <- boots(alpha.hat, sigma.hat)
    if("try-error" %in% class(rout.star)) stop("reaster error (Nelder-Mead with pickle fail?)")
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  
  tryCatch({
    vandse <- VaW(rout.star) # function for VaW calcualtion at bootstrap level 1
    KW_Va.star[iboot] <- vandse[1] # function for VaW calcualtion at bootstrap level 1
    se.KW_Va.star_sd[iboot]<- vandse[2] # calcualte standard deviation
    if("try-error" %in% class(vandse)) stop("aster.predict error (direction of recession error?)") 
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  
  print(KW_Va.star[iboot]) # will print outer bootstrap estimate for each nboot interation
  print(se.KW_Va.star_sd[iboot]) # will print outer bootstrap estimate for each nboot interation
  print(iboot)
  
  # The current code generates a standard deviation scale metric
  
}

# GENERATE BOOTSTRAP t CONFIDENCE INTERVALS


# Save Rdata for entire bootstrap run
#save.image("singleboot_KW2016_run.RData")

load("singleboot_KW2016_run.RData")

# Generate Va(W) from data
KW_Va.hat <- VaW(rout2016)[1] # generates va from data
KW_Va.hat_sd <- VaW(rout2016)[2] # using the empirical sd.

# set confidence level
conf.level <- 0.95

# calculate t for standard deviation-based confidence interval
t <- (KW_Va.star - KW_Va.hat) / se.KW_Va.star_sd

# generate critical values with t for standard deviation-based confidence interval
crit <- quantile(t, probs = c((1 - conf.level) / 2, (1 + conf.level) / 2))

# Conf. Int. for KW_Va with standard deviation scale metric
KW_Va.hat - rev(crit) * KW_Va.hat_sd