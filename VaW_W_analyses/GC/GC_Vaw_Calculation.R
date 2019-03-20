
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

