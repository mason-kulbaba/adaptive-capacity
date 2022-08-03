
#The following code performs reaster analyses that will be used to calculate the 
#additive genetic variance for fitness (Va.W) and mean expected fitness (W) for 
#the pedigreed "greenhouse" cohort in the Grey Cloud Dunes site for three years (2015-2017)


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

#Load 3-year data file for Grey Cloud Dunes

setwd("C:/Users/Mason Kulbaba/Dropbox/git/adaptive-capacity/Revised/GC")

gcdat<- read.csv("gcdata.csv")

#divide into year-specific files
gc2015<- subset(gcdat, year==2015)
gc2016<- subset(gcdat, year==2016)
gc2017<- subset(gcdat, year==2017)

#drop unused levels
gc2015<- droplevels(gc2015)
gc2016<- droplevels(gc2016)
gc2017<- droplevels(gc2017)


###################################################################################
# 1. Estimate mean expected fitness in greenhouse genration for 2015 and 2016 year#
###################################################################################

#isolate greenhouse generation from 2016 and 2017 data (only greenhouse cohort in 2015)
gh2016<- subset(gc2016, cohort=="greenhouse")
gh2016<- droplevels(gh2016)

gh2017<- subset(gc2017, cohort=="greenhouse")
gh2017<- droplevels(gh2017) 


#make sure maternalID and paternalID is a factor
gc2015$maternalID<- as.factor(gc2015$maternalID)
gc2015$paternalID<- as.factor(gc2015$paternalID)

gh2016$maternalID<- as.factor(gh2016$maternalID)
gh2016$paternalID<- as.factor(gh2016$paternalID)

gh2017$maternalID<- as.factor(gh2017$maternalID)
gh2017$paternalID<- as.factor(gh2017$paternalID)

#set response variables -> these represent variables in graphical model
vars<- c("Germ","flw","total.pods", "total.pods.collected", "totalseeds")


#reshape data so that all response variables are located in a single vector in a new data
#set called "redata"
redata2015 <- reshape(gc2015, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")

redata2016 <- reshape(gh2016, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")

redata2017 <- reshape(gh2017, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")

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

#save(redata2015, file="redata2015.RData")
#save(redata2016, file="redata2016.RData")
#save(redata2017, file="redata2017.RData")

#load aster package
library(aster)

#set graphical mode and dist. for fitness nodes (preds)
pred<- c(0,1,2,3,4)
fam<- c(1,1,2,1,2) # Note: subset node is under a Bernoulli distribution. 

#describe dist. of preds.
#sapply(fam.default(), as.character)[fam]

#fixed effect model for 2015
#aout2015a<- aster(resp~varb, pred, fam, varb, id, root, data=redata2015)

#aout2015<- aster(resp~varb + fit:(block), pred, fam, varb, id, root, data=redata2015)

#anova(aout2015a, aout2015)

#summary(aout2015, show.graph = T)


#Designate fitness variable for 2016 data
fit <- grepl("totalseeds", as.character(redata2016$varb))
fit<- as.numeric(fit)

redata2016$fit<- fit

#check
#with(redata2016, sort(unique(as.character(varb)[fit == 0])))
#with(redata2016, sort(unique(as.character(varb)[fit == 1])))

#Fixed-effects aster model for 2016 data
#aout2016a<- aster(resp~varb , pred, fam, varb, id, root, data=redata2016)

#aout2016<- aster(resp~varb + fit:(block), pred, fam, varb, id, root, data=redata2016)

#anova(aout2016a, aout2016)

#summary(aout2016, show.graph = T)


#Designation of fitness variable for 2017 data
fit <- grepl("totalseeds", as.character(redata2017$varb))
fit<- as.numeric(fit)

redata2017$fit <- fit

#check
#with(redata2017, sort(unique(as.character(varb)[fit == 0])))
#with(redata2017, sort(unique(as.character(varb)[fit == 1])))


#add a variable "root" to redata files, where value is 1
redata2017<- data.frame(redata2017, root=1)


#make sure block, row, and position are factors
redata2017$block<- as.factor(redata2017$block)
redata2017$row<- as.factor(redata2017$row)
redata2017$position<- as.factor(redata2017$position)

#check
with(redata2017, sort(unique(as.character(varb)[fit == 0])))
with(redata2017, sort(unique(as.character(varb)[fit == 1])))

#Fixed-effects aster model for 2016 data
#aout2017a<- aster(resp~varb, pred, fam, varb, id, root, data=redata2017)

#aout2017<- aster(resp~varb + fit:(block), pred, fam, varb, id, root, data=redata2017)

#anova(aout2017a, aout2017)

#summary(aout2017, show.graph = T)


##############################################################
#Estimate mean fitness for 2015, 2016 (greenhouse cohort), and 2017 (greenhouse cohort) data

#generate MLE of saturated model mean value parameter vector: mu
#pout2015<- predict.aster(aout2015, se.fit=TRUE)
#pout2016<- predict.aster(aout2016, se.fit=TRUE)
#pout2017<- predict.aster(aout2017, se.fit=TRUE)


#make up covariate data for hypothetical
#indivs that meet "typical" criteria:
#Therefore, "make up" covariate data for hypothetical individuals
#that re comparable and obtain mean values for them
##############

#make data.frame of indivudals for each block (1-8)
#fred2015 <- data.frame(block=levels(redata2015$block),
 #                  Germ=1, flw=1, total.pods=1, total.pods.collected=1, totalseeds=1,root = 1)

#fred2016 <- data.frame(block=levels(redata2016$block),
  #                     Germ=1, flw=1, total.pods=1, total.pods.collected=1, totalseeds=1,root = 1)

#fred2017 <- data.frame(block=levels(redata2017$block),
   #                    Germ=1, flw=1, total.pods=1, total.pods.collected=1, totalseeds=1,root = 1)


#reshape the "made up data" just as the actual data
#renewdata2015 <- reshape(fred2015, varying = list(vars),
    #                 direction = "long", timevar = "varb",
     #                times = as.factor(vars), v.names = "resp")

#renewdata2016 <- reshape(fred2016, varying = list(vars),
     #                    direction = "long", timevar = "varb",
      #                   times = as.factor(vars), v.names = "resp")

#renewdata2017 <- reshape(fred2017, varying = list(vars),
       #                  direction = "long", timevar = "varb",
        #                 times = as.factor(vars), v.names = "resp")

#make character string from "varb" of renewdata,
#without actual values (e.g., the layers of varb in renewdata)
#layer<- gsub("[0-9]", "", as.character(renewdata2015$varb))
#layer<- gsub("[0-9]", "", as.character(renewdata2016$varb))
#layer<- gsub("[0-9]", "", as.character(renewdata2017$varb))

#add layer to renewdata
#renewdata2015<- data.frame(renewdata2015, layer= layer)
#renewdata2016<- data.frame(renewdata2016, layer= layer)
#renewdata2017<- data.frame(renewdata2017, layer= layer)

# add totalseeds in new layer col of renewdata as numeric, called fit
#fit<- as.numeric(layer=="totalseeds")

#add fit to renewdata
#renewdata2015<- data.frame(renewdata2015, fit = fit)
#renewdata2016<- data.frame(renewdata2016, fit = fit)
#renewdata2017<- data.frame(renewdata2017, fit = fit)


#Generate fintess estimates and standard errors for each block
#nblock<- nrow(fred2015)#all data has same number of blocks so any file will do
#nnode<- length(vars)
#amat<- array(0, c(nblock, nnode, nblock))
#dim(amat)# makes an 8 x 5 x 8 matrix

#only want means for k'th individual that contributes to expected
#fitness, and want to add only totalseeds entries

#makes 8 , 4x8 matrices that 

#foo<- grepl("totalseeds", vars)
#for(k in 1:nblock)
#  amat[k, foo, k]<- 1

#generate predicted valuses using aout object, with renewdata, and amat format
#pout.amat2015<- predict(aout2015, newdata= renewdata2015, varvar= varb,
 #                   idvar= id, root = root, se.fit=TRUE, amat = amat)

#pout.amat2016<- predict(aout2016, newdata= renewdata2016, varvar= varb,
 #                       idvar= id, root = root, se.fit=TRUE, amat = amat)

#pout.amat2017<- predict(aout2017, newdata= renewdata2017, varvar= varb,
 #                       idvar= id, root = root, se.fit=TRUE, amat = amat)

#combine std.err with estimates, and then round
#to three decimal places
#gc2015<- cbind(pout.amat2015$fit, pout.amat2015$se.fit)
#gc2016<- cbind(pout.amat2016$fit, pout.amat2016$se.fit)
#gc2017<- cbind(pout.amat2017$fit, pout.amat2017$se.fit)

#rownames(gc2015)<- as.character(fred2015$block)
#rownames(gc2016)<- as.character(fred2016$block)
#rownames(gc2017)<- as.character(fred2017$block)

#colnames(gc2015)<- c("Expected Fitness", "SE")
#colnames(gc2016)<- c("Expected Fitness", "SE")
#colnames(gc2017)<- c("Expected Fitness", "SE")


#####################################################################
#The median expected fitness for each year is used to represent     #
#the site-specific fitness. The block effects representing these    #
# fitness values will be used to calcualte and convert  Va(W)       #
#estimates from the canonical to man-value parameter scale          #
#####################################################################

#round(gc2015, 3) # median: block 6 - 1.112
#round(gc2016, 3) # median: block 6 - 0.640
#round(gc2017, 3) # median: block 7 - 1.062


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

#save(modmat.siredam2015, file="modmat.siredam2015.RData")
#save(modmat.siredam2016, file="modmat.siredam2016.RData")
#save(modmat.siredam2017, file="modmat.siredam2017.RData")

#Reaster analyses

#NOTE: we will use the same graphical model and statistical distributions that we used
#      in the fixed-effects aster models to estiamte mean expected fitness



#2015 analysis
#rout2015a<- reaster(resp~varb,list(parental=~0 + modmat.siredam2015),
 #                  +pred, fam, varb, id, root, data=redata2015) 

#rout2015b_sub<- reaster(resp~ fit + varb +fit:(block),list(parental=~0 + modmat.siredam2015),
 #               +pred, fam, varb, id, root, data=redata2015) 

# anova(rout2015a, rout2015b)

#sum.2015<- summary(rout2015b_sub)

#sum.2015

#save(rout2015b_sub, file="rout2015b_sub.RData")

load("rout2015b_sub.Rdata")

#2016 analysis
#rout2016a<- reaster(resp~varb,list(parental=~0 + modmat.siredam2016),
 #                  +pred, fam, varb, id, root, data=redata2016)

#rout2016b_sub<- reaster(resp~ fit+ varb +fit:(block),list(parental=~0 + modmat.siredam2016),
 #                   +pred, fam, varb, id, root, data=redata2016)
# anova(rout2016a, rout2016b)


#sum.2016<- summary(rout2016b_sub) 

#sum.2016

#save(rout2016b_sub, file="rout2016b_sub.RData")

load("rout2016b_sub.Rdata")

#2017 analysis
#rout2017a<- reaster(resp~varb,list(parental=~0 + modmat.siredam2017),
 #                  +pred, fam, varb, id, root, data=redata2017) 

#rout2017b_sub<- reaster(resp~ fit + varb +fit:(block),list(parental=~0 + modmat.siredam2017),
 #                  +pred, fam, varb, id, root, data=redata2017)


# anova(rout2017a, rout2017b)

#sum.2017<- summary(rout2017b_sub)

#sum.2017

#save(rout2017b_sub, file="rout2017b_sub.RData")

load("rout2017b_sub.Rdata")

#####################################
#
# The following code is to use the new mapping function, that correctly accounts for
# subsetting (thinning), and calculated correct breeding values
#


#2015 data

#pulls up the fixed effect aster model needed for predict.aster 
class(rout2015b_sub$obj) 
bhat <- rout2015b_sub$b 
bhat.sire<- bhat[grep("paternalID", names(bhat))]

# 42 sires
length(bhat.sire)

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
  modmat <- rout$obj$modmat[1880, ,  , drop = FALSE]# 1880 first indiv in block 6A
 # modmat[ , , "fit:block6A"] <- 0
  
   # set root = 1
  root <- array(1, dim = dim(modmat)[1:2])#
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

plot(density(bhat.sire.mu), main = "", xlab = "total seeds")

GC2015b<- bhat.sire.mu

##########################################################################
#
# 2016 breeding values
#

#pulls up the fixed effect aster model needed for predict.aster 
class(rout2016b_sub$obj) 
bhat <- rout2016b_sub$b 
bhat.sire<- bhat[grep("paternalID", names(bhat))]

# 42 sires
length(bhat.sire)

#extract bhat - the estimates of random effects
bhat.2016<- rout2016b_sub$b

#extract the sire effects
bhat.sire2016<- bhat.2016[grep("paternalID", names(bhat))]

#stem plot of bhat estimates
stem(bhat.sire2016)

# This is the revised mapping function to *properly* account for subsampling
map_factory <- function(rout) {
  stopifnot(inherits(rout, "reaster"))
  alpha <- rout$alpha
  # modmat for one individual
  modmat <- rout$obj$modmat[2251, , ,drop = FALSE] # 2251 fist individual in block6B
 # modmat[ , , "fit:block6B"] <- 0
  
   # set root = 1
  root <- array(1, dim = dim(modmat)[1:2])#
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
bhat.sire.mu <- vectorized.map(bhat.sire2016)

#Plot the graph of the function `map`.

curve(vectorized.map, from = -1/4, to = 1/4,
      xlab="b", ylab=expression(mu(b)))

#Plot the density of estimated breeding values mapped to mean value parameter scale.

plot(density(bhat.sire.mu), main = "", xlab = "total seeds")

GC2016b<- bhat.sire.mu


##########################################################################
#
# 2017 breeding values
#

#pulls up the fixed effect aster model needed for predict.aster 
class(rout2017b_sub$obj) 
bhat <- rout2017b_sub$b 
bhat.sire<- bhat[grep("paternalID", names(bhat))]

# 42 sires
length(bhat.sire)

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
  modmat <- rout$obj$modmat[3126, ,  , drop = FALSE] # 3126 first indiv in block7C
 # modmat[ , , "fit:block7C"] <- 0
  
  # set root = 1
  root <- array(1, dim = dim(modmat)[1:2])#
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

plot(density(bhat.sire.mu), main = "", xlab = "total seeds")

GC2017b<- bhat.sire.mu

# Correlations between years

cor(GC2015b, GC2016b) #0.1303375
cor(GC2015b, GC2017b) # 0.5338379
cor(GC2016b, GC2017b)# 0.08058884

####################################
#
# Combine all three years of breeding values, and print csv file

GCb<- rbind(GC2015b, GC2016b, GC2017b)

write.table(GCb, file="GCb.csv", sep=",")
