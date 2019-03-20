
# The following code calculates the additive genetic effects on total lifetime fitness
# from the three reaster output files generated in: CS_W_reaster_analyses.R

setwd("C:/Users/Mason Kulbaba/Dropbox/Rscripts/Final_VaW_Wbar/CS")

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
#parameter scale (from canoonical scale) to fit:block1 fixed effect (alpha[6]= block1),
#predict, and take value of fourth node for first indiv.
map <- function(b) {
  stopifnot(length(b) == 1)
  stopifnot(is.finite(b))
  alpha <- rout2015$alpha
  alpha[6] <- alpha[6] + b #block 1 effects
  hoom <- predict(rout2015$obj, newcoef = alpha)
  hoom <- matrix(hoom, ncol = 5)
  return(hoom[1, 5])# individual in block 1, fifth node
}

#vectorize mapping function
fred<- Vectorize(map)

#plot the curve
curve(fred, from = -1 / 2, to = 1/2, xlab="b", ylab=expression(mu(b))) 

#use the mapping function "fred" to convert the sire effects from the 
#canonical to mean-value parameter scale
bhat.sireGC.mu<- fred(bhat.sire)

cs_b2015<- bhat.sireGC.mu

stem(bhat.sireGC.mu)#no longer normal as expected

#generate probability density distribution
cs2015den<-density(bhat.sireGC.mu)

cs2015den<- cbind(cs2015den[[1]], cs2015den[[2]])

write.csv(cs2015den, "cs2015den.csv", quote = FALSE, row.names = FALSE)


hoom<- predict(rout2015$obj, newcoef = rout2015$alpha, se.fit=TRUE)
goom <- hoom$gradient
moom<- goom[,5]
moom<- matrix(moom, ncol=5)

#this is additive genetic variation for fitness!
CS_Va2015<- 4*moom[1 ,5]^2 * rout2015$nu[1]/map(0) 
CS_Va2015 #3.204069

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
#parameter scale (from canoonical scale) to fit:block1 fixed effect (alpha[7]= block2),
#predict, and take value of fourth node for first indiv.
map <- function(b) {
  stopifnot(length(b) == 1)
  stopifnot(is.finite(b))
  alpha <- rout2016$alpha
  alpha[7] <- alpha[7] + b #block 2 effects
  hoom <- predict(rout2016$obj, newcoef = alpha)
  hoom <- matrix(hoom, ncol = 5)
  return(hoom[251, 5])# individual in block 2, fifth node
}

#vectorize mapping function
fred<- Vectorize(map)

#plot the curve
curve(fred, from = -1 / 2, to = 1/2, xlab="b", ylab=expression(mu(b)))

#use the mapping function "fred" to convert the sire effects from the 
#canonical to mean-value parameter scale
bhat.sireGC.mu<- fred(bhat.sire)

cs_b2016<- bhat.sireGC.mu

stem(bhat.sireGC.mu)#no longer normal as expected

#generate probability density distribution
cs2016den<- density(bhat.sireGC.mu)

cs2016den<- cbind(cs2016den[[1]], cs2016den[[2]])

write.csv(cs2016den, "cs2016den.csv", quote = FALSE, row.names = FALSE)



hoom<- predict(rout2016$obj, newcoef = rout2016$alpha, se.fit=TRUE)
goom <- hoom$gradient
moom<- goom[,5]
moom<- matrix(moom, ncol=5)

#this is additive genetic variation for fitness!
CS_Va2016<- 4*moom[251 ,5]^2 * rout2016$nu[1]/map(0)
CS_Va2016 #0.8562831


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
#parameter scale (from canoonical scale) to fit:block1 fixed effect (alpha[8]= block3),
#predict, and take value of fourth node for first indiv.
map <- function(b) {
  stopifnot(length(b) == 1)
  stopifnot(is.finite(b))
  alpha <- rout2017$alpha
  alpha[8] <- alpha[8] + b #block 3 effects
  hoom <- predict(rout2017$obj, newcoef = alpha)
  hoom <- matrix(hoom, ncol = 5)
  return(hoom[971, 5])# individual in block 3, fifth node
}

#vectorize mapping function
fred<- Vectorize(map)

#plot the curve
curve(fred, from = -1 / 2, to = 1/2, xlab="b", ylab=expression(mu(b)))  

#use the mapping function "fred" to convert the sire effects from the 
#canonical to mean-value parameter scale
bhat.sireGC.mu<- fred(bhat.sire)

cs_b2017<- bhat.sireGC.mu

cs_bs<- rbind(cs_b2015, cs_b2016, cs_b2017)

write.table(cs_bs, file="CSb_threeyears.csv", sep=",")

stem(bhat.sireGC.mu)#no longer normal as expected

#generate probability density distribution
cs2017den<- density(bhat.sireGC.mu)

cs2017den<- cbind(cs2017den[[1]], cs2017den[[2]])

write.csv(cs2017den, "cs2017den.csv", quote = FALSE, row.names = FALSE)


hoom<- predict(rout2017$obj, newcoef = rout2017$alpha, se.fit=TRUE)
goom <- hoom$gradient
moom<- goom[,5]
moom<- matrix(moom, ncol=5)

#this is additive genetic variation for fitness!
CS_Va2017<- 4*moom[971 ,5]^2 * rout2017$nu[1]/map(0) 
CS_Va2017 #1.572956