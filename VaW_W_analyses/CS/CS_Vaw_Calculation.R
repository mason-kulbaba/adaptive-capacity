
# The following code calculates the additive genetic effects on total lifetime fitness
# from the three reaster output files generated in: CS_W_reaster_analyses.R

setwd("C:/Users/Mason Kulbaba/Dropbox/git/adaptive-capacity/VaW_W_analyses/CS")

#Load reaster 

#2015 reaster data
load(file="rout2015b.RData")

#2016 reaster data
load(file="rout2016b.RData")

#2017 reaster data
load(file="rout2017b.RData")

library(aster)


#########################
#Estimate Va(W) for 2015#
#########################

#extract bhat - the estimates of random effects
bhat<- rout2015b$b

#extract the sire effects
bhat.sire<- bhat[grep("paternalID", names(bhat))]

#stem plot of bhat estimates
stem(bhat.sire)#the canonical estimates "look" somewhat normal...as they should

hoom <- predict(rout2015b$obj, newcoef = rout2015b$alpha)
hoom<- matrix(hoom, ncol =5)
hoom<- hoom[ , 5]


#add a component of b_hat that we want to map to the mean value
#parameter scale (from canoonical scale) to fit:block1 fixed effect (alpha[6]= block1),
#predict, and take value of fourth node for first indiv.
map <- function(b) {
  stopifnot(length(b) == 1)
  stopifnot(is.finite(b))
  alpha <- rout2015b$alpha
  alpha[6] <- alpha[6] + b #block 1 effects
  hoom <- predict(rout2015b$obj, newcoef = alpha)
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


hoom<- predict(rout2015b$obj, newcoef = rout2015b$alpha, se.fit=TRUE)
goom <- hoom$gradient
moom<- goom[,5]
moom<- matrix(moom, ncol=5)

#this is additive genetic variation for fitness!
CS_Va2015<- 4*moom[1 ,5]^2 * rout2015b$nu[1]/map(0) 
CS_Va2015 #3.204069


#recall fitness in 2015 = 1.80

#predicted change in fitness
CS_Va2015/1.80


sout<- summary(rout2015b)

#standard error for predicted change in fitness
4*moom[1 ,5]^2 * sout$nu["parental", "Std. Error"]/1.80


#########################
#Estimate Va(W) for 2016#
#########################

#extract bhat - the estimates of random effects
bhat<- rout2016b$b

#extract the sire effects
bhat.sire<- bhat[grep("paternalID", names(bhat))]

#stem plot of bhat estimates
stem(bhat.sire)#the canonical estimates "look" somewhat normal...as they should

hoom <- predict(rout2016b$obj, newcoef = rout2016b$alpha)
hoom<- matrix(hoom, ncol =5)
hoom<- hoom[ , 5]


#add a component of b_hat that we want to map to the mean value
#parameter scale (from canoonical scale) to fit:block1 fixed effect (alpha[7]= block2),
#predict, and take value of fourth node for first indiv.
map <- function(b) {
  stopifnot(length(b) == 1)
  stopifnot(is.finite(b))
  alpha <- rout2016b$alpha
  alpha[7] <- alpha[7] + b #block 2 effects
  hoom <- predict(rout2016b$obj, newcoef = alpha)
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



hoom<- predict(rout2016b$obj, newcoef = rout2016b$alpha, se.fit=TRUE)
goom <- hoom$gradient
moom<- goom[,5]
moom<- matrix(moom, ncol=5)

#this is additive genetic variation for fitness!
CS_Va2016<- 4*moom[251 ,5]^2 * rout2016b$nu[1]/map(0)
CS_Va2016 #0.8562831

sout<-summary(rout2016b)

#recall mean fitness in 2016 = 0.73

#predicted change in fitness
CS_Va2016/0.73

#standard error of predicted change in fitness

4*moom[251 ,5]^2 * sout$nu["parental", "Std. Error"]/0.73

#########################
#Estimate Va(W) for 2017#
#########################

#extract bhat - the estimates of random effects
bhat<- rout2017b$b

#extract the sire effects
bhat.sire<- bhat[grep("paternalID", names(bhat))]

#stem plot of bhat estimates
stem(bhat.sire)#the canonical estimates "look" somewhat normal...as they should

hoom <- predict(rout2017b$obj, newcoef = rout2017b$alpha)
hoom<- matrix(hoom, ncol =5)
hoom<- hoom[ , 5]


#add a component of b_hat that we want to map to the mean value
#parameter scale (from canoonical scale) to fit:block1 fixed effect (alpha[8]= block3),
#predict, and take value of fourth node for first indiv.
map <- function(b) {
  stopifnot(length(b) == 1)
  stopifnot(is.finite(b))
  alpha <- rout2017b$alpha
  alpha[8] <- alpha[8] + b #block 3 effects
  hoom <- predict(rout2017b$obj, newcoef = alpha)
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


hoom<- predict(rout2017b$obj, newcoef = rout2017b$alpha, se.fit=TRUE)
goom <- hoom$gradient
moom<- goom[,5]
moom<- matrix(moom, ncol=5)

#this is additive genetic variation for fitness!
CS_Va2017<- 4*moom[971 ,5]^2 * rout2017b$nu[1]/map(0) 
CS_Va2017 #1.572956

sout<-summary(rout2017b)

#recall mean fitness in 2017 = 1.24

#predicted chagne in fitness
CS_Va2017/1.24

#standard error of prediction
4*moom[971 ,5]^2 * sout$nu["parental", "Std. Error"]/1.24
