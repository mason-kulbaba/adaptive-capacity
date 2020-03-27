
# The following code calculates the additive genetic effects on total lifetime fitness
# from the three reaster output files generated in: KW_W_reaster_analyses.R

setwd("C:/Users/Mason Kulbaba/Dropbox/git/adaptive-capacity/VaW_W_analyses/KW")

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
#parameter scale (from canoonical scale) to fit:block1 fixed effect (alpha[9]= block4),
#predict, and take value of fourth node for first indiv.
map <- function(b) {
  stopifnot(length(b) == 1)
  stopifnot(is.finite(b))
  alpha <- rout2015b$alpha
  alpha[9] <- alpha[9] + b #block 4 effects
  hoom <- predict(rout2015b$obj, newcoef = alpha)
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


hoom<- predict(rout2015b$obj, newcoef = rout2015b$alpha, se.fit=TRUE)
goom <- hoom$gradient
moom<- goom[,5]
moom<- matrix(moom, ncol=5)

#this is additive genetic variation for fitness!
KW_Va2015<- 4*moom[1570 ,5]^2 * rout2015b$nu[1]/map(0) 
KW_Va2015 #3.527548

sout<-summary(rout2015b)

#recall mean fitness in 2015 = 3.02

#predicted change in fitness 
KW_Va2015/3.02

#standard error for pred. change in fitness
4*moom[1570 ,5]^2 * sout$nu["parental", "Std. Error"]/3.02


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
#parameter scale (from canoonical scale) to fit:block1 fixed effect (alpha[10]= block5),
#predict, and take value of fourth node for first indiv.
map <- function(b) {
  stopifnot(length(b) == 1)
  stopifnot(is.finite(b))
  alpha <- rout2016b$alpha
  alpha[10] <- alpha[10] + b #block 5 effects
  hoom <- predict(rout2016b$obj, newcoef = alpha)
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



hoom<- predict(rout2016b$obj, newcoef = rout2016b$alpha, se.fit=TRUE)
goom <- hoom$gradient
moom<- goom[,5]
moom<- matrix(moom, ncol=5)

#this is additive genetic variation for fitness!
KW_Va2016<- 4*moom[2251 ,5]^2 * rout2016b$nu[1]/map(0)
KW_Va2016 #1.571878

#recall mean fitness in 2016 = 1.88

#predicted change in fitness

KW_Va2016/1.88


sout<- summary(rout2016b)

#standard error for predicted change in fitness
4*moom[2251 ,5]^2 * sout$nu["parental", "Std. Error"]/1.88


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
#parameter scale (from canoonical scale) to fit:block1 fixed effect (alpha[11]= block6),
#predict, and take value of fourth node for first indiv.
map <- function(b) {
  stopifnot(length(b) == 1)
  stopifnot(is.finite(b))
  alpha <- rout2017b$alpha
  alpha[11] <- alpha[11] + b #block 6 effects
  hoom <- predict(rout2017b$obj, newcoef = alpha)
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


hoom<- predict(rout2017b$obj, newcoef = rout2017b$alpha, se.fit=TRUE)
goom <- hoom$gradient
moom<- goom[,5]
moom<- matrix(moom, ncol=5)

#this is additive genetic variation for fitness!
KW_Va2017<- 4*moom[2303 ,5]^2 * rout2017b$nu[1]/map(0) 
KW_Va2017 #1.409605

sout<-summary(rout2017b)

#recall mean fitness in 2017 = 1.08

#predicted change in fitness
KW_Va2017/1.08

#standard error
4*moom[2303 ,5]^2 * sout$nu["parental", "Std. Error"]/1.08

