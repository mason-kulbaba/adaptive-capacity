
# The following code calculates the additive genetic effects on total lifetime fitness
# from the three reaster output files generated in: GC_W_reaster_analyses.R

setwd()

#Load reaster 

#2015 reaster data
 load(file="rout2015.RData")

#2016 reaster data
load(file="rout2016.RData")

#2017 reaster data
load(file="rout2017.RData")

library(aster)


#########################
#Estimate Va(W) for 2015#
#########################

#extract bhat - the estimates of random effects
bhat<- rout2015$b

#extract the sire effects
bhat.sire<- bhat[grep("paternalID", names(bhat))]

#stem plot of bhat estimates
stem(bhat.sire)#the canonical estimates "look" somewhat normal...as they should

hoom <- predict(rout2015$obj, newcoef = rout2015$alpha)
hoom<- matrix(hoom, ncol =5)
hoom<- hoom[ , 5]


#add a component of b_hat that we want to map to the mean value
#parameter scale (from canoonical scale) to fit:block1 fixed effect (alpha[11]= block6),
#predict, and take value of fourth node for first indiv.
map <- function(b) {
  stopifnot(length(b) == 1)
  stopifnot(is.finite(b))
  alpha <- rout2015$alpha
  alpha[11] <- alpha[11] + b #block 6 effects
  hoom <- predict(rout2015$obj, newcoef = alpha)
  hoom <- matrix(hoom, ncol = 5)
  return(hoom[1874, 5])#first individual in block 6, fifth node
}

#vectorize mapping function
fred<- Vectorize(map)

#plot the curve
curve(fred, from = -1 / 2, to = 1/2, xlab="b", ylab=expression(mu(b))) 

#use the mapping function "fred" to convert the sire effects from the 
#canonical to mean-value parameter scale
bhat.sireGC.mu<- fred(bhat.sire)

GC_b2015<- bhat.sireGC.mu

stem(bhat.sireGC.mu)#no longer normal as expected

#generate probability density distribution
gc2015den<-density(bhat.sireGC.mu)

GC2015_den<- cbind(gc2015den[[1]], gc2015den[[2]])

#write.csv(GC2015_den, "gc2015den.csv", quote = FALSE, row.names = FALSE)


hoom<- predict(rout2015$obj, newcoef = rout2015$alpha, se.fit=TRUE)
goom <- hoom$gradient
moom<- goom[,5]
moom<- matrix(moom, ncol=5)

#this is additive genetic variation for fitness!
GC_Va<- 4*moom[1874 ,5]^2 * rout2015$nu[1]/map(0) 
GC_Va #1.8576

#recall mean fitness in 2015 = 1.11

#predicted change in fitness from FFTNS:

GC_Va/1.11

#alternate way to get Va(W)
hoom <- predict(rout2015$obj, newcoef = rout2015$alpha, gradient = TRUE)
goom <- hoom$gradient



moom <- goom[ , 5]
moom <- matrix(moom, ncol = 5)
moom[1874, 5]

#very close
4 * moom[1874, 5]^2 * rout2015$nu[1]




#Standard Error for predicted change in fitness

sout<- summary(rout2015)

sout

4*moom[1874 ,5]^2 * sout$nu["parental", "Std. Error"]/1.11


#OK, very nice. Now see how this standard error compares with that generated from a
#more "rigorous" method

#inv. Fisher info.

fisher <- sout$fisher
eout <- eigen(fisher, symmetric = TRUE)

fisher.inv <- eout$vectors %*% diag(1/eout$values) %*% t(eout$vectors)

dim(fisher)

dim(fisher.inv)

length(rout2015$alpha)

length(rout2015$nu)

#check - NOTE: changed from [6, 6] that was in Charlie's origional code
sqrt(fisher.inv[13, 13])

#good are the same
sout$nu["parental", "Std. Error"]




ioom <- seq(along = hoom$fit)
ioom <- matrix(ioom, ncol = 5)
idx <- ioom[1874, 5]
idx

g <- hoom$gradient[idx, ]


#check
predict(rout2015$obj, newcoef = rout2015$alpha, se.fit = TRUE)$se.fit[idx]

sqrt(t(g) %*% solve(rout2015$obj$fisher) %*% g)

#apply delta method for predicted change in fitness
g <- c(g, 0, 0)


xi.hat <- rout2015$nu[1]
eta.hat <- foom[1874, 5]


v.nu <- fisher.inv[13, 13]


v.meanfit <- t(g) %*% fisher.inv %*% g
v.nu.meanfit <- (fisher.inv %*% g)[13]







#########################
#Estimate Va(W) for 2016#
#########################

#extract bhat - the estimates of random effects
bhat<- rout2016$b

#extract the sire effects
bhat.sire<- bhat[grep("paternalID", names(bhat))]

#stem plot of bhat estimates
stem(bhat.sire)#the canonical estimates "look" somewhat normal...as they should

hoom <- predict(rout2016$obj, newcoef = rout2016$alpha)
hoom<- matrix(hoom, ncol =5)
hoom<- hoom[ , 5]


#add a component of b_hat that we want to map to the mean value
#parameter scale (from canoonical scale) to fit:block1 fixed effect (alpha[11]= block6),
#predict, and take value of fourth node for first indiv.
map <- function(b) {
  stopifnot(length(b) == 1)
  stopifnot(is.finite(b))
  alpha <- rout2016$alpha
  alpha[11] <- alpha[11] + b #block 6 effects
  hoom <- predict(rout2016$obj, newcoef = alpha)
  hoom <- matrix(hoom, ncol = 5)
  return(hoom[2274, 5])#first individual in block 6, fifth node
}

#vectorize mapping function
fred<- Vectorize(map)

#plot the curve
curve(fred, from = -1 / 2, to = 1/2, xlab="b", ylab=expression(mu(b)))

#use the mapping function "fred" to convert the sire effects from the 
#canonical to mean-value parameter scale
bhat.sireGC.mu<- fred(bhat.sire)

GC_b2016<- bhat.sireGC.mu

stem(bhat.sireGC.mu)#no longer normal as expected

#generate probability density distribution
gc2016den<- density(bhat.sireGC.mu)

GC2016_den<- cbind(gc2016den[[1]], gc2016den[[2]])

#write.csv(GC2016_den, "gc2016den.csv", quote = FALSE, row.names = FALSE)



hoom<- predict(rout2016$obj, newcoef = rout2016$alpha, se.fit=TRUE)
goom <- hoom$gradient
moom<- goom[,5]
moom<- matrix(moom, ncol=5)

#this is additive genetic variation for fitness!
GC_Va<- 4*moom[2274 ,5]^2 * rout2016$nu[1]/map(0)
GC_Va #0.83021

#recall mean fitness in 2016 = 0.64

#Predicted change in fitness

GC_Va/0.64

sout<-summary(rout2016)

4*moom[2274 ,5]^2 * sout$nu["parental", "Std. Error"]/0.64


#########################
#Estimate Va(W) for 2017#
#########################

#extract bhat - the estimates of random effects
bhat<- rout2017$b

#extract the sire effects
bhat.sire<- bhat[grep("paternalID", names(bhat))]

#stem plot of bhat estimates
stem(bhat.sire)#the canonical estimates "look" somewhat normal...as they should

hoom <- predict(rout2017$obj, newcoef = rout2017$alpha)
hoom<- matrix(hoom, ncol =5)
hoom<- hoom[ , 5]


#add a component of b_hat that we want to map to the mean value
#parameter scale (from canoonical scale) to fit:block1 fixed effect (alpha[12]= block7),
#predict, and take value of fourth node for first indiv.
map <- function(b) {
  stopifnot(length(b) == 1)
  stopifnot(is.finite(b))
  alpha <- rout2017$alpha
  alpha[12] <- alpha[12] + b #block 7 effects
  hoom <- predict(rout2017$obj, newcoef = alpha)
  hoom <- matrix(hoom, ncol = 5)
  return(hoom[3126, 5])#first individual in block 7, fourth node
}

#vectorize mapping function
fred<- Vectorize(map)

#plot the curve
curve(fred, from = -1 / 2, to = 1/2, xlab="b", ylab=expression(mu(b)))  

#use the mapping function "fred" to convert the sire effects from the 
#canonical to mean-value parameter scale
bhat.sireGC.mu<- fred(bhat.sire)

GC_b2017<- bhat.sireGC.mu

GCb<- rbind(GC_b2015, GC_b2016, GC_b2017)

#write.table(GCb, "GCb_threeyears.csv", sep=",")

stem(bhat.sireGC.mu)#no longer normal as expected

#generate probability density distribution
gc2017den<- density(bhat.sireGC.mu)

GC2017_den<- cbind(gc2017den[[1]], gc2017den[[2]])

#write.csv(GC2017_den, "gc2017den.csv", quote = FALSE, row.names = FALSE)


hoom<- predict(rout2017$obj, newcoef = rout2017$alpha, se.fit=TRUE)
goom <- hoom$gradient
moom<- goom[,5]
moom<- matrix(moom, ncol=5)

#this is additive genetic variation for fitness!
GC_Va<- 4*moom[3126 ,5]^2 * rout2017$nu[1]/map(0) 
GC_Va #6.49131

#recall mean fitness in 2017 = 1.06

#predicted change in fitness

GC_Va/1.06


sout<-summary(rout2017)

4*moom[3126 ,5]^2 * sout$nu["parental", "Std. Error"]/1.06
