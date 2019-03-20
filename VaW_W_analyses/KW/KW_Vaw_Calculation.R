
# The following code calculates the additive genetic effects on total lifetime fitness
# from the three reaster output files generated in: KW_W_reaster_analyses.R

setwd("C:/Users/Mason Kulbaba/Dropbox/Rscripts/chamaecrista-adaptive-capacity/VaW_W_analyses/KW")

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
#parameter scale (from canoonical scale) to fit:block1 fixed effect (alpha[9]= block4),
#predict, and take value of fourth node for first indiv.
map <- function(b) {
  stopifnot(length(b) == 1)
  stopifnot(is.finite(b))
  alpha <- rout2015$alpha
  alpha[9] <- alpha[9] + b #block 4 effects
  hoom <- predict(rout2015$obj, newcoef = alpha)
  hoom <- matrix(hoom, ncol = 5)
  return(hoom[1570, 5])# individual in block 4, fourth node
}

#vectorize mapping function
fred<- Vectorize(map)

#plot the curve
curve(fred, from = -1 / 2, to = 1/2, xlab="b", ylab=expression(mu(b))) 

#use the mapping function "fred" to convert the sire effects from the 
#canonical to mean-value parameter scale
bhat.sireGC.mu<- fred(bhat.sire)

kw_b2015<- bhat.sireGC.mu

stem(bhat.sireGC.mu)#no longer normal as expected

#generate probability density distribution
kw2015den<-density(bhat.sireGC.mu)

kw2015den<- cbind(kw2015den[[1]], kw2015den[[2]])

write.csv(kw2015den, "kw2015den.csv", quote = FALSE, row.names = FALSE)


hoom<- predict(rout2015$obj, newcoef = rout2015$alpha, se.fit=TRUE)
goom <- hoom$gradient
moom<- goom[,5]
moom<- matrix(moom, ncol=5)

#this is additive genetic variation for fitness!
KW_Va2015<- 4*moom[1570 ,5]^2 * rout2015$nu[1]/map(0) 
KW_Va2015 #3.527548


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
#parameter scale (from canoonical scale) to fit:block1 fixed effect (alpha[10]= block5),
#predict, and take value of fourth node for first indiv.
map <- function(b) {
  stopifnot(length(b) == 1)
  stopifnot(is.finite(b))
  alpha <- rout2016$alpha
  alpha[10] <- alpha[10] + b #block 5 effects
  hoom <- predict(rout2016$obj, newcoef = alpha)
  hoom <- matrix(hoom, ncol = 5)
  return(hoom[2251, 5])# individual in block 5, fourth node
}

#vectorize mapping function
fred<- Vectorize(map)

#plot the curve
curve(fred, from = -1 / 2, to = 1/2, xlab="b", ylab=expression(mu(b)))

#use the mapping function "fred" to convert the sire effects from the 
#canonical to mean-value parameter scale
bhat.sireGC.mu<- fred(bhat.sire)


kw_b2016<- bhat.sireGC.mu

stem(bhat.sireGC.mu)#no longer normal as expected

#generate probability density distribution
kw2016den<- density(bhat.sireGC.mu)

kw2016den<- cbind(kw2016den[[1]], kw2016den[[2]])

write.csv(kw2016den, "kw2016den.csv", quote = FALSE, row.names = FALSE)



hoom<- predict(rout2016$obj, newcoef = rout2016$alpha, se.fit=TRUE)
goom <- hoom$gradient
moom<- goom[,5]
moom<- matrix(moom, ncol=5)

#this is additive genetic variation for fitness!
KW_Va2016<- 4*moom[2251 ,5]^2 * rout2016$nu[1]/map(0)
KW_Va2016 #1.571878


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
#parameter scale (from canoonical scale) to fit:block1 fixed effect (alpha[11]= block6),
#predict, and take value of fourth node for first indiv.
map <- function(b) {
  stopifnot(length(b) == 1)
  stopifnot(is.finite(b))
  alpha <- rout2017$alpha
  alpha[11] <- alpha[11] + b #block 6 effects
  hoom <- predict(rout2017$obj, newcoef = alpha)
  hoom <- matrix(hoom, ncol = 5)
  return(hoom[2302, 5])# individual in block 6, fourth node
}

#vectorize mapping function
fred<- Vectorize(map)

#plot the curve
curve(fred, from = -1 / 2, to = 1/2, xlab="b", ylab=expression(mu(b)))  

#use the mapping function "fred" to convert the sire effects from the 
#canonical to mean-value parameter scale
bhat.sireGC.mu<- fred(bhat.sire)

kw_b2017<- bhat.sireGC.mu

kw_bs<- rbind(kw_b2015, kw_b2016, kw_b2017)

write.table(kw_bs, file="KWb_threeyears.csv", sep=",")

stem(bhat.sireGC.mu)#no longer normal as expected

#generate probability density distribution
kw2017den<- density(bhat.sireGC.mu)

kw2017den<- cbind(kw2017den[[1]], kw2017den[[2]])

write.csv(kw2017den, "kw2017den.csv", quote = FALSE, row.names = FALSE)


hoom<- predict(rout2017$obj, newcoef = rout2017$alpha, se.fit=TRUE)
goom <- hoom$gradient
moom<- goom[,5]
moom<- matrix(moom, ncol=5)

#this is additive genetic variation for fitness!
KW_Va2017<- 4*moom[2303 ,5]^2 * rout2017$nu[1]/map(0) 
KW_Va2017 #1.409605
