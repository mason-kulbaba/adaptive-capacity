library(aster)
load("doubleboot_GC2016_runA.RData")
flist <- ls()
va_star <- GC_Va.star
se.va_star <- se.GC_Va.star
se.va_star_sd <- se.GC_Va.star_sd
rm(list = flist)
load("doubleboot_GC2016_runD.RData")
va_star <- c(va_star,GC_Va.star)
se.va_star <- c(se.va_star,se.GC_Va.star)
se.va_star_sd <- c(se.va_star_sd, se.GC_Va.star_sd)
rm(list = flist)
load("doubleboot_GC2016_runFHSB.RData")
va_star <- c(va_star,GC_Va.star)
se.va_star <- c(se.va_star,se.GC_Va.star)
se.va_star_sd <- c(se.va_star_sd, se.GC_Va.star_sd)
rm(list = flist)
load("doubleboot_GC2016_runFHSE.RData")
va_star <- c(va_star,GC_Va.star)
se.va_star <- c(se.va_star,se.GC_Va.star)
se.va_star_sd <- c(se.va_star_sd, se.GC_Va.star_sd)
rm(list = flist)
load("doubleboot_GC2016_runB.RData")
va_star <- c(va_star,GC_Va.star)
se.va_star <- c(se.va_star,se.GC_Va.star)
se.va_star_sd <- c(se.va_star_sd, se.GC_Va.star_sd)
rm(list = flist)
load("doubleboot_GC2016_runE.RData")
va_star <- c(va_star,GC_Va.star)
se.va_star <- c(se.va_star,se.GC_Va.star)
se.va_star_sd <- c(se.va_star_sd, se.GC_Va.star_sd)
rm(list = flist)
load("doubleboot_GC2016_runFHSC.RData")
va_star <- c(va_star,GC_Va.star)
se.va_star <- c(se.va_star,se.GC_Va.star)
se.va_star_sd <- c(se.va_star_sd, se.GC_Va.star_sd)
rm(list = flist)
load("doubleboot_GC2016_runFHSF.RData")
va_star <- c(va_star,GC_Va.star)
se.va_star <- c(se.va_star,se.GC_Va.star)
se.va_star_sd <- c(se.va_star_sd, se.GC_Va.star_sd)
rm(list = flist)
load("doubleboot_GC2016_runC.RData")
va_star <- c(va_star,GC_Va.star)
se.va_star <- c(se.va_star,se.GC_Va.star)
se.va_star_sd <- c(se.va_star_sd, se.GC_Va.star_sd)
rm(list = flist)
load("doubleboot_GC2016_runFHSA.RData")
va_star <- c(va_star,GC_Va.star)
se.va_star <- c(se.va_star,se.GC_Va.star)
se.va_star_sd <- c(se.va_star_sd, se.GC_Va.star_sd)
# rm(list = flist)  Now we have them all and want to keep the various variables

########################################################
# check
########################################################

length(va_star)

length(se.va_star)

length(se.va_star_sd)


# assuming these are all 100, proceed


GC_Va.star <- va_star
se.GC_Va.star <- se.va_star
se.GC_Va.star_sd <- se.va_star_sd

#### Back into the original CI code
GC_Va.hat <- VaW(rout2016) # generates Va(W) from data

# calculate t for IQR-based confidence interval
t <- (GC_Va.star - GC_Va.hat) / se.GC_Va.star

# set confidence level
conf.level <- 0.95

# generate critical values with t for IQR-based confidence interval
crit <- quantile(t, probs = c((1 - conf.level) / 2, (1 + conf.level) / 2))  # note this uses the t just calculated

# generate IQR and standard deviation measures of scale
foo <- GC_Va.star    # the outer loop Va's
se.GC_Va.hat <- diff(sort(foo)[c(25,76)])  # middle 50 percent of 100 or so

se.GC_Va.hat_sd <- sd(foo)  # standard deviation calculation of scale

# Conf. Int. for GC_Va with IQR scale metric
GC_Va.hat + crit * se.GC_Va.hat

# calculate t for standard deviation-based confidence interval
t <- (GC_Va.star - GC_Va.hat) / se.GC_Va.star_sd

# generate critical values with t for standard deviation-based confidence interval
crit <- quantile(t, probs = c((1 - conf.level) / 2, (1 + conf.level) / 2))

# Conf. Int. for GC_Va with standard deviation scale metric
GC_Va.hat + crit * se.GC_Va.hat_sd

proc.time()


