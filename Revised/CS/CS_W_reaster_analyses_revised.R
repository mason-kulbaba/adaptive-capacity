
#The following code performs reaster analyses that will be used to calculate the 
#additive genetic variance for fitness (Va.W) and mean expected fitness (W) for 
#the pedigreed "greenhouse" cohort in the CERA site for three years (2015-2017)


#Steps to the analyiss

#   1. Estimate mean expected fitness in "greenhouse" generations, and identify block used
#      to represent the site. This block will then be used in the calcualtion of Va(W).

#   2. Reaster models to produce output using block identified in step 1.

#   3. G x Time (environment) analysis


#Note: the same graphical model is used in all above three steps.

#The varibles occur in the following order:

# 1.  "Germ" -> germination (0=no, 1=yes), Bernoulli
# 2.  "flw" -> survival to flowering (0=no, 1=yes), Bernoulli
# 3.  "total.pods" -> total number of pods produced, Poisson
# 4.  "total.pods.collected" -> number of pods collected (subsampled pods), Poisson
# 5.  "totalseeds" -> total number of seeds counted from collected pods, Poisson



#################

#Begin with data preparation

#Load 3-year data file for CERA

setwd("C:/Users/Mason Kulbaba/Dropbox/git/adaptive-capacity/Revised/CS")

csdat<- read.csv("csdata.csv")

#divide into year-specific files
cs2015<- subset(csdat, year==2015)
cs2016<- subset(csdat, year==2016)
cs2017<- subset(csdat, year==2017)

#drop unused levels
cs2015<- droplevels(cs2015)
cs2016<- droplevels(cs2016)
cs2017<- droplevels(cs2017)


###################################################################################
# 1. Estimate mean expected fitness in greenhouse genration for 2015 and 2016 year#
###################################################################################

#isolate greenhouse generation from 2016 and 2017 data (only greenhouse cohort in 2015)
cs2016<- subset(cs2016, cohort=="greenhouse")
cs2016<- droplevels(cs2016)

cs2017<- subset(cs2017, cohort=="greenhouse")
cs2017<- droplevels(cs2017)


#make sure maternalID and paternalID is a factor
cs2015$maternalID<- as.factor(cs2015$maternalID)
cs2015$paternalID<- as.factor(cs2015$paternalID)

cs2016$maternalID<- as.factor(cs2016$maternalID)
cs2016$paternalID<- as.factor(cs2016$paternalID)

cs2017$maternalID<- as.factor(cs2017$maternalID)
cs2017$paternalID<- as.factor(cs2017$paternalID)

#set response variables -> these represent variables in graphical model
vars<- c("Germ","flw","total.pods", "total.pods.collected", "totalseeds")


#reshape data so that all response variables are located in a single vector in a new data
#set called "redata"
redata2015 <- reshape(cs2015, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")

redata2016 <- reshape(cs2016, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")

redata2017 <- reshape(cs2017, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")

#Designation of fitness variable for 2015 data
fit <- grepl("totalseeds", as.character(redata2015$varb))
fit<- as.numeric(fit)

redata2015$fit <- fit

#check
with(redata2015, sort(unique(as.character(varb)[fit == 0])))
with(redata2015, sort(unique(as.character(varb)[fit == 1])))


#add a variable "root" to redata files, where value is 1
redata2015<- data.frame(redata2015, root=1)
redata2016<- data.frame(redata2016, root=1)


#make sure block, row, and position are factors
redata2015$block<- as.factor(redata2015$block)
redata2015$row<- as.factor(redata2015$row)
redata2015$position<- as.factor(redata2015$position)

redata2016$block<- as.factor(redata2016$block)
redata2016$row<- as.factor(redata2016$row)
redata2016$position<- as.factor(redata2016$position)

#load aster package
library(aster)

#set graphical mode and dist. for fitness nodes (preds)
pred<- c(0,1,2,3,4)
fam<- c(1,1,2,1,2)

#describe dist. of preds.
sapply(fam.default(), as.character)[fam]



#Designate fitness variable for 2016 data
fit <- grepl("totalseeds", as.character(redata2016$varb))
fit<- as.numeric(fit)


#Designation of fitness variable for 2017 data
fit <- grepl("totalseeds", as.character(redata2017$varb))
fit<- as.numeric(fit)

redata2017$fit <- fit

#check
with(redata2017, sort(unique(as.character(varb)[fit == 0])))
with(redata2017, sort(unique(as.character(varb)[fit == 1])))


#add a variable "root" to redata files, where value is 1
redata2017<- data.frame(redata2017, root=1)


#make sure block, row, and position are factors
redata2017$block<- as.factor(redata2017$block)
redata2017$row<- as.factor(redata2017$row)
redata2017$position<- as.factor(redata2017$position)

#check
with(redata2017, sort(unique(as.character(varb)[fit == 0])))
with(redata2017, sort(unique(as.character(varb)[fit == 1])))




#############################################################################

# 2. Reaster analyses for 2015, 2016, and 2017 greenhouse cohort. For clarity, the 
#    actual calculation of Va(W) will be performed in a different script (Vaw_Calculation.R).


####combine maternal and paternal into a single parental effect to be estimated
modmat.sire2015 <- model.matrix(~ 0 + fit:paternalID, redata2015)
modmat.dam2015 <- model.matrix(~ 0 + fit:maternalID, redata2015)
modmat.siredam2015 <- cbind(modmat.sire2015,modmat.dam2015)

modmat.sire2016 <- model.matrix(~ 0 + fit:paternalID, redata2016)
modmat.dam2016 <- model.matrix(~ 0 + fit:maternalID, redata2016)
modmat.siredam2016 <- cbind(modmat.sire2016,modmat.dam2016)

modmat.sire2017 <- model.matrix(~ 0 + fit:paternalID, redata2017)
modmat.dam2017 <- model.matrix(~ 0 + fit:maternalID, redata2017)
modmat.siredam2017 <- cbind(modmat.sire2017,modmat.dam2017)

#save(redata2015, file="redata2015.RData")
#save(redata2016, file="redata2016.RData")
#save(redata2017, file="redata2017.RData")


#save(modmat.siredam2015, file="modmat.siredam2015.RData")
#save(modmat.siredam2016, file="modmat.siredam2016.RData")
#save(modmat.siredam2017, file="modmat.siredam2017.RData")


#Reaster analyses

#NOTE: we will use the same graphical model and statistical distributions that we used
#      in the fixed-effects aster models to estiamte mean expected fitness



#2015 analysis

rout2015b_sub<- reaster(resp~ fit + varb +fit:(block),list(parental=~0 + modmat.siredam2015),
                +pred, fam, varb, id, root, data=redata2015) 


sum.2015<- summary(rout2015b_sub)

sum.2015

save(rout2015b_sub, file="rout2015b_sub.RData")

load("rout2015b_sub.RData")

#2016 analysis


rout2016b_sub<- reaster(resp~ fit + varb +fit:(block),list(parental=~0 + modmat.siredam2016),
                    +pred, fam, varb, id, root, data=redata2016)


sum.2016<- summary(rout2016b_sub) 

sum.2016

save(rout2016b_sub, file="rout2016b_sub.RData")


#2017 analysis

rout2017b_sub<- reaster(resp~ fit + varb +fit:(block),list(parental=~0 + modmat.siredam2017),
                   +pred, fam, varb, id, root, data=redata2017) 


sum.2017<- summary(rout2017b_sub)

sum.2017

save(rout2017b_sub, file="rout2017b_sub.RData")

#######################################################################################
#pulls up the fixed effect aster model needed for predict.aster 
class(rout2015b_sub$obj) 
bhat <- rout2015b_sub$b 
bhat.sire<- bhat[grep("paternalID", names(bhat))]

#21 sires
length(bhat.sire)
#################################################################################AP

#extract bhat - the estimates of random effects
bhat<- rout2015b_sub$b

#extract the sire effects
bhat.sire<- bhat[grep("paternalID", names(bhat))]

#stem plot of bhat estimates
stem(bhat.sire)



# This is the revised mapping function to *properly* account for subsampling
map_factory <- function(rout) {
  stopifnot(inherits(rout, "reaster"))
  alpha <- rout$alpha
  # modmat for one individual
  modmat <- rout$obj$modmat[1, ,  , drop = FALSE] # first indiv. in block1A
  # set root = 1
  root <- array(1, dim = dim(modmat)[1:2])
  is.subsampling <- grepl("total.pods.collected", vars) # remove subsample node?
  function (b) {
    stopifnot(is.numeric(b))
    stopifnot(length(b) == 1)
    if (! is.finite(b)) return(NaN)
    alpha["fit"] <- alpha["fit"] + b # should these also correspond to fit"block1A, too?
    # predict.aster doesn't use argument x when is.always.parameter = TRUE
    # but it still checks x.  So x must be supplied and valid even though
    # it does not affect the result.
    xi <- try(predict(rout$obj, newcoef = alpha, modmat = modmat,
                      root = root, model.type = "conditional",
                      is.always.parameter = TRUE, x = root), silent = TRUE)
    if (inherits(xi, "try-error")) return(NaN)
    return(prod(xi[! is.subsampling]))
  }
}

map <- map_factory(rout2015b_sub)

#Try it

#Be sure to vectorize the map function
vectorized.map <- Vectorize(map)

#produce breeding values on the mean value parameter scale
bhat.sire.mu <- vectorized.map(bhat.sire)

#Plot the graph of the function `map`.

curve(vectorized.map, from = -1/4, to = 1/4,
      xlab="b", ylab=expression(mu(b)))

#Plot the density of estimated breeding values mapped to mean value parameter scale.

plot(density(bhat.sire.mu), main = "", xlab = "total seeds") # x-axis seems to have very large values


cs.2015<- bhat.sire.mu

cs.2015<- cs.2015[2:21]


# Calculate Va(W)
integrand_factory <- function(rout) {
  stopifnot(inherits(rout, "reaster"))
  sigma <- as.numeric(rout$sigma) * 2
  function (b) {
    stopifnot(is.numeric(b))
    foo <- vectorized.map(b) * b * dnorm(b, sd = sigma)
    foo[is.na(foo)] <- 0
    return(foo)
  }
}
integrand <- integrand_factory(rout2015b_sub)

iout <- integrate(integrand, lower = -Inf, upper = Inf)
iout$message
iout$value
# cov(m(Z), Z)^2 / var(Z)
vaw.hat <- iout$value^2 / environment(integrand)$sigma^2
vaw.hat



##############################################################################
#
# 2016 analysis


load("rout2016b_sub.RData")


#pulls up the fixed effect aster model needed for predict.aster 
class(rout2016b_sub$obj) 
bhat <- rout2016b_sub$b 
bhat.sire<- bhat[grep("paternalID", names(bhat))]

#21 sires
length(bhat.sire)
#################################################################################AP

#extract bhat - the estimates of random effects
bhat<- rout2016b_sub$b

#extract the sire effects
bhat.sire<- bhat[grep("paternalID", names(bhat))]

#stem plot of bhat estimates
stem(bhat.sire)



# This is the revised mapping function to *properly* account for subsampling
map_factory <- function(rout) {
  stopifnot(inherits(rout, "reaster"))
  alpha <- rout$alpha
  # modmat for one individual
  modmat <- rout$obj$modmat[1, ,  , drop = FALSE] # 251 first indiv in block2B
  modmat[ , , "fit:block2B"] <- 0
  # set root = 1
  root <- array(1, dim = dim(modmat)[1:2])
  is.subsampling <- grepl("total.pods.collected", vars) # remove subsample node?
  function (b) {
    stopifnot(is.numeric(b))
    stopifnot(length(b) == 1)
    if (! is.finite(b)) return(NaN)
    alpha["fit"] <- alpha["fit"] + b # should these also correspond to fit"block1A, too?
    # predict.aster doesn't use argument x when is.always.parameter = TRUE
    # but it still checks x.  So x must be supplied and valid even though
    # it does not affect the result.
    xi <- try(predict(rout$obj, newcoef = alpha, modmat = modmat,
                      root = root, model.type = "conditional",
                      is.always.parameter = TRUE, x = root), silent = TRUE)
    if (inherits(xi, "try-error")) return(NaN)
    return(prod(xi[! is.subsampling]))
  }
}

map <- map_factory(rout2016b_sub)

#Try it

#Be sure to vectorize the map function
vectorized.map <- Vectorize(map)

#produce breeding values on the mean value parameter scale
bhat.sire.mu <- vectorized.map(bhat.sire)

#Plot the graph of the function `map`.

curve(vectorized.map, from = -1/4, to = 1/4,
      xlab="b", ylab=expression(mu(b)))

#Plot the density of estimated breeding values mapped to mean value parameter scale.

plot(density(bhat.sire.mu), main = "", xlab = "total seeds") # x-axis seems to have very large values

cs.2016<- bhat.sire.mu

cs.2016<- cs.2016[2:21]


# Calculate Va(W)
integrand_factory <- function(rout) {
  stopifnot(inherits(rout, "reaster"))
  sigma <- as.numeric(rout$sigma) * 2
  function (b) {
    stopifnot(is.numeric(b))
    foo <- vectorized.map(b) * b * dnorm(b, sd = sigma)
    foo[is.na(foo)] <- 0
    return(foo)
  }
}
integrand <- integrand_factory(rout2016b_sub)

iout <- integrate(integrand, lower = -Inf, upper = Inf)
iout$message
iout$value
# cov(m(Z), Z)^2 / var(Z)
vaw.hat <- iout$value^2 / environment(integrand)$sigma^2
vaw.hat

##############################################################################
#
# 2017 analysis

#pulls up the fixed effect aster model needed for predict.aster 
class(rout2017b_sub$obj) 
bhat <- rout2017b_sub$b 
bhat.sire<- bhat[grep("paternalID", names(bhat))]

#21 sires
length(bhat.sire)
#################################################################################AP

#extract bhat - the estimates of random effects
bhat<- rout2017b_sub$b

#extract the sire effects
bhat.sire<- bhat[grep("paternalID", names(bhat))]

#stem plot of bhat estimates
stem(bhat.sire)



# This is the revised mapping function to *properly* account for subsampling
map_factory <- function(rout) {
  stopifnot(inherits(rout, "reaster"))
  alpha <- rout$alpha
  # modmat for one individual
  modmat <- rout$obj$modmat[971, ,  , drop = FALSE] # 971 first indiv in block 3C
  # set root = 1
  root <- array(1, dim = dim(modmat)[1:2])
  is.subsampling <- grepl("total.pods.collected", vars) # remove subsample node?
  function (b) {
    stopifnot(is.numeric(b))
    stopifnot(length(b) == 1)
    if (! is.finite(b)) return(NaN)
    alpha["fit"] <- alpha["fit"] + b # should these also correspond to fit"block1A, too?
    # predict.aster doesn't use argument x when is.always.parameter = TRUE
    # but it still checks x.  So x must be supplied and valid even though
    # it does not affect the result.
    xi <- try(predict(rout$obj, newcoef = alpha, modmat = modmat,
                      root = root, model.type = "conditional",
                      is.always.parameter = TRUE, x = root), silent = TRUE)
    if (inherits(xi, "try-error")) return(NaN)
    return(prod(xi[! is.subsampling]))
  }
}

map <- map_factory(rout2017b_sub)

#Try it

#Be sure to vectorize the map function
vectorized.map <- Vectorize(map)

#produce breeding values on the mean value parameter scale
bhat.sire.mu <- vectorized.map(bhat.sire)

#Plot the graph of the function `map`.

curve(vectorized.map, from = -1/4, to = 1/4,
      xlab="b", ylab=expression(mu(b)))

#Plot the density of estimated breeding values mapped to mean value parameter scale.

plot(density(bhat.sire.mu), main = "", xlab = "total seeds") # x-axis seems to have very large values

cs.2017<- bhat.sire.mu


##################################################
#
# combine and write file of breeding values

CSb<- rbind(cs.2015, cs.2016, cs.2017)

write.table(CSb, "CSb.csv", sep=",")

##########################################################
#
# Correlations

cor(cs.2015, cs.2016) # 0.38024
cor(cs.2015, cs.2017) # 0.29467
cor(cs.2016, cs.2017) # 0.45765
