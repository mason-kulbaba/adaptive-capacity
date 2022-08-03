
#The following code troubleshoots the subsampling issue in aster models

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

# Germ -> flw -> total.pods -> total.pods.collected -> totalseeds

pred<- c(0,1,2,3,4)
fam<- c(1,1,2,1,2) # note use Bernoulli for subsampling node

#describe dist. of preds.
sapply(fam.default(), as.character)[fam]




#############################################################################

#  Reaster analyses for 2015, 2016, and 2017 greenhouse cohort. For clarity, the 
#  actual calculation of Va(W) will be performed in a different script (Vaw_Calculation.R).


####combine maternal and paternal into a single parental effect to be estimated
modmat.sire2015 <- model.matrix(~ 0 + fit:paternalID, redata2015)
modmat.dam2015 <- model.matrix(~ 0 + fit:maternalID, redata2015)
modmat.siredam2015 <- cbind(modmat.sire2015,modmat.dam2015)

#NOTE: from AP "For convenience of mapping from canonical to mean value parameter scale, 
# we force `fit` into the model by putting it early in the fixed effects formula."

#   I do not understand this reasoning. How does this make mapping more convenient? 


rout2015b_sub<- reaster(resp~ fit + varb + fit:(block) ,list(parental=~0 + modmat.siredam2015),
                pred, fam, varb, id, root, data=redata2015) 



rout.sum<-summary(rout2015b_sub)

#check model summary

rout.sum

#save(rout2015b_sub, file="rout2015b_sub.RData")

load("rout2015b_sub.RData")

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
stem(bhat.sire)#the canonical estimates appear somewhat normal...as they should.

###############################################
#
# This is the "old" mapping function, where we used block1 fixed effects for mapping alpha[6]
#

#map <- function(b) {
#  stopifnot(length(b) == 1)
#  stopifnot(is.finite(b))
#  alpha <- rout2015b$alpha
#  alpha[6] <- alpha[6] + b #block 1 effects
#  hoom <- predict(rout2015b_sub$obj, newcoef = alpha)
#  hoom <- matrix(hoom, ncol = 6)
#  return(hoom[1, 6])# individual in block 1, fifth node
#    }



# This is the revised mapping function to *properly* account for subsampling
map_factory <- function(rout) {
  stopifnot(inherits(rout, "reaster"))
  alpha <- rout$alpha
  # modmat for one individual
  modmat <- rout$obj$modmat[1, ,  , drop = FALSE]
  # Set Yloc = 0 - presume this would be 'block1A' in our analysis
  modmat[ , , "(Intercept)"] <- 0 # fit:block1A has been aliased into intercept
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

