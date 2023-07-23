##This script includes a function to compare simulation to empirical data and sets up a workflow to 
##sweep parameter values and produce a parameter space map.

##load libraries

library(sp)
library(rgdal)
library(maptools)
library(rgeos)
library(colourvalues)
library(scales)
library(readxl)
library(lattice)

##FUNCTION DEFINITIONS

##define diffusion function
diffusion <- function(x, wm, a, cap=1){
  res <- rep(NA,length(x))
  for (i in 1:nrow(wm)){
    res[i] <- (x[i]*a) + sum(wm[i,] * x)
  } 
  res[res > 1] <- cap
  return(res)
}


##define diffusion function with random trait loss where traits are more or less likely to be lost
##loss values are drawn from a uniform distribution
##@param maximum is the value for the maximum loss value allowed

diffusion.nloss <- function(x, wm, a, cap=1, maximum = 0.1){
  res <- rep(NA,length(x))
  for (i in 1:nrow(wm)){
    res[i] <- (x[i]*a) + sum(wm[i,] * x)
    res[i] <- res[i] - runif(n = 1, min = 0, max = maximum)*res[i] #trait loss
  } 
  res[res > 1] <- cap
  return(res)
}


##define distance decay function

distfun.exp <- function(x,b){exp(-b*x)} #Distance decay function



##Function that compares simulated to empirical timeseries using the sum of absolute differences (function):
##The function has as arguments two time series as two vectors and does:
##1.Produce a new vector with the absolute value of the difference between each element of the input vectors.
##2.Sum the elements of that vector to produce a single value.

compare.fun <- function(empir, sim){
  compar <- abs(empir - sim)
  compar <- sum(compar)
  return(compar)
}


##load simulation data/parameter values
nnodes <- 2000 # number of nodes
seeds <- c(1)
phlen <- c(175,258,142,163,49,138,38,72,145,55,35,57,128,133)
dur <- sum(phlen)
tstep <- dur
alpha <- c(1)
kappa <- c(0.4)
loss <- 1 #trait loss allowed when loss = 1, set loss = 0 for no loss
mu <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)

##load spatial data
directory <- getwd()
aegean <- readOGR(dsn=directory, layer="aeg") #vector map
utm35Nwgs84 <- CRS("+init=epsg:32635") #set crs
aegean <- spTransform(aegean, utm35Nwgs84) #project to specified crs

##load innovation data

##load data
wheel_data <- read_excel("wheel_percentages_may2022_tidy.xlsx")
wheel_data0 <- wheel_data[,c(2:15)] #select only percentages per phase
wheel_data1 <- as.matrix(wheel_data0)

##turn "na" into NA values

for (j in 1:ncol(wheel_data1)){
  for (i in 1:nrow(wheel_data1)){
    if (wheel_data1[i,j] == "na"){
      wheel_data1[i,j] <- NA
    }
    else {
      wheel_data1[i,j] <- wheel_data1[i,j]
    }
  }
}

##turn all values numeric

vou <- apply(wheel_data1, c(1, 2), as.numeric)

##produce time series

tseries <- apply(vou, 2, mean, na.rm = TRUE) ##mean
tseriesmed <- apply(vou, 2, median, na.rm = TRUE) ##median
tseries_sd <- apply(vou, 2, sd, na.rm = TRUE) ##standard deviations


##Load simulated origin points: set up origin node matrix

##This sets up 9 origin points

tleft <- SpatialPoints(coords = data.frame(-50009.3, 4690051), proj4string=CRS("+init=epsg:32635"))
mleft <- SpatialPoints(coords = data.frame(-50009.3, 4200000), proj4string=CRS("+init=epsg:32635"))
bleft <- SpatialPoints(coords = data.frame(200009.3, 3920000), proj4string=CRS("+init=epsg:32635"))
bcenter <- SpatialPoints(coords = data.frame(300009.3, 3880000), proj4string=CRS("+init=epsg:32635"))
bright <- SpatialPoints(coords = data.frame(800009.3, 4040000), proj4string=CRS("+init=epsg:32635"))
mright <- SpatialPoints(coords = data.frame(800009.3, 4200000), proj4string=CRS("+init=epsg:32635"))
tright <- SpatialPoints(coords = data.frame(800009.3, 4500000), proj4string=CRS("+init=epsg:32635"))
tcenter <- SpatialPoints(coords = data.frame(300009.3, 4670051), proj4string=CRS("+init=epsg:32635"))
center <- SpatialPoints(coords = data.frame(300009.3, 4200000), proj4string=CRS("+init=epsg:32635"))
all <- rbind(tleft, mleft, bleft, bcenter, bright, mright, tright, tcenter, center)  

##SIMULATION - Double nested loop, outer loop sweeps alpha, inner loop sweeps kappa
set.seed(seeds[1])

maxloss_list <- vector(mode = "list", length = length(mu))
maxloss_results_list <- vector(mode = "list", length = length(mu))
for (m in 1:length(mu)){

meansmatrixlistlist <- vector(mode = "list", length = length(alpha))
mediansmatrixlistlist <- vector(mode = "list", length = length(alpha))
sdmatrixlistlist <- vector(mode = "list", length = length(alpha))
resultsllist <- vector(mode = "list", length = length(alpha))
nranl <- vector(mode = "list", length = length(alpha))
for (a in 1:length(alpha)){
  nran0 <- spsample(aegean, n=nnodes, type="random")
  nran00 <- nran0[-1] #remove one point from nran0 so that when I add the origin I will still have a pop = nnodes
  nran <- rbind(nran00, all[7]) ##add index of origin point from the 'all' vector
  nran$siteid <- row.names(nran) ##add siteid column to nran
  distmat <- gDistance(nran, byid=TRUE) 
  weights <- distmat/1000 #meters to km
  resultsl <- vector(mode = "list", length = length(kappa))
  weightl <- vector(mode = "list", length = length(kappa))
  for (l in 1:length(kappa)){
    weightl[[l]] <- apply(weights, 1, FUN=distfun.exp, b=kappa[l])
    diag(weightl[[l]]) <- 0 #so that nodes don't increase on their own
    
    trait <- rep(c(0),each = nnodes)
    trait[nnodes] <- 0.1 # this sets the origin point trait proportion value
    trait0 <- trait
    
    results <- matrix(NA, nrow = nnodes, ncol = tstep) # results matrix (size is determined by parameter values set at the beginning)
    results[,1] <- trait0
    for (p in 2:ncol(results)){
      if(loss == 0){
        trait <- diffusion(x = trait, wm = weightl[[l]], a=alpha[a])
      } else {
        trait <- diffusion.nloss(x = trait, wm = weightl[[l]], a=alpha[a], maximum = mu[m])
      }
      results[,p] <- trait
    }
    resultsl[[l]] <- results
    print("Done")
  }
  
  
  ##lists of means and medians for different population values
  meansmatrixlist <- vector(mode = "list", length = length(kappa))
  for (b in 1:length(kappa)){
    meansmatrixlist[[b]] <- apply(resultsl[[b]], 2, mean) * 100
  }
  
  sdmatrixlist <- vector(mode = "list", length = length(kappa))
  for (b in 1:length(kappa)){
    sdmatrixlist[[b]] <- apply(resultsl[[b]], 2, sd) * 100
  }  
  
  meansmatrixlistlist[[a]] <- meansmatrixlist
  sdmatrixlistlist[[a]] <- sdmatrixlist
  resultsllist[[a]] <- resultsl
  nranl[[a]] <- nran
}

maxloss_list[[m]] <- meansmatrixlistlist
maxloss_results_list[[m]] <- resultsllist
}




##This section fits a logistic curve to simulated trends to explore the model's sensitivity to the parameters
## It fits logistic model using SSlogis, a method that self starts fitting a logistic function.


##Function that fits SSlogis
## the if statement controls for the case where the self starting model gives an error
## the else case (where there is no error) fits the logistic model
## the function's output is the scaling value used for each case


extract.scaling.fun <- function(trend, tsteps, thresh = 0.1){
  if (trend[length(trend)] < thresh){ ## this bit ensures we won't get an error when the simulation doesn't look s shaped
    scal <- NA
  }
  else {
    logisticModelSS <- nls(trend~SSlogis(tsteps, Asym, xmid, scal), na.action = na.omit,control=nls.control(maxiter=1000))
    scal <- coef(logisticModelSS)[3] 
  }
  return(scal)
}

##apply to raw data

##first, create rawmeans object string the trends for diff values of mu, alpha and kappa without aggregation
rawmeans <- list()
for (m in 1:length(mu)){
    ##apply compare.fun function to the data
    rawmeans[[m]] <- maxloss_list[[m]][[1]][[1]]/100
}


##then, create a matrix storing each scaling param value for the simulation of each mu/kappa combination

scal_mat <- matrix(NA, nrow = length(mu), ncol = 1)
for (m in 1:length(mu)){
    scal_mat[m,1] <- extract.scaling.fun(trend = rawmeans[[m]], tsteps = 1:1588, thresh = 0.1)
  }
 

##plot param space map

##pdf("sweep_scaling_loss.pdf")

colnames(scal_mat) <- as.character(kappa)
rownames(scal_mat) <- as.character(mu)
lattice::levelplot(scal_mat, xlab="m",
                   ylab="kappa",col.regions = terrain.colors(100))

##dev.off()


## This section compareS empirical to simulation via monte carlo simulation of the empirical trend


phend <- c() ##this object has labels for the year (BC) that each phase ends on
phend[1] <- 2475
for (i in 2:length(phlen)){
  phend[i] <- phend[i-1] - phlen[i]
}


vec <- c(0,phlen)
vec2 <- c()
for (i in 1:length(vec)){
  vec2[i] <- sum(vec[1:i])
}
vec2 <- vec2[2:15]



mcsamptemp <- vector(mode = "list", length = 1000)
for (i in 1:length(mcsamptemp)){
  citser <- vector(mode = "list", length = 14)
  for (j in 1:length(citser)){
    citser[[j]] <- c(tseries[j],
                     round(runif(1, min = phend[j], max = phend[j] + phlen[j])))
  }
  mcsamptemp[[i]] <- do.call("rbind",citser)
}



##make sure monte carlo simulated empirical proportions do not exceed 100 % and are not below 0 %

for (i in 1:length(mcsamptemp)){
  for (j in 1:14){
    if (mcsamptemp[[i]][j,1] > 100){
      mcsamptemp[[i]][j,1] = 100
    }
    if (mcsamptemp[[i]][j,1] < 0){
      mcsamptemp[[i]][j,1] = 0
    }
    else {
      mcsamptemp[[i]][j,1] = mcsamptemp[[i]][j,1]
    }
  }
}


##for each of the 1000 simulations, get a scaling parameter value

mc_scalings <- c()
for (i in 1:length(mcsamptemp)){
  mc_scalings[i] <- extract.scaling.fun(trend = mcsamptemp[[i]][,1]/100, tsteps =mcsamptemp[[i]][,2])
} ##scaling values are negative because the BC date inverts the x axis


## This plots a histogram of the scaling values

hist(-mc_scalings) ## mc scalings are normally distributed, the - symbol helps re-invert the time axis

## compare empirical mc scalings to the scaling produced by the one of the simulations-
## each comparison produces a distribution of scaling differences

scal_diff <- matrix(list(), nrow=length(mu), ncol=length(kappa))
for (j in 1:ncol(scal_mat)){
  for (i in 1:nrow(scal_mat)){
    if (is.na(scal_mat[i,j])){
      scal_diff[[i,j]] <- NA
    }
    else {
      scal_diff[[i,j]] <- abs(scal_mat[i,j] - (-mc_scalings))
    }
  }
}

##Plot parameter space map where each cell is the mean of the differences between the 1000 mc sim for the empirical and the simulated

scal_diff_mean <- matrix(NA, nrow=length(mu), ncol=length(kappa))
for (j in 1:ncol(scal_diff)){
  for (i in 1:nrow(scal_diff)){
    if (is.na(scal_diff[[i,j]])){
      scal_diff_mean[i,j] <- NA
    } else {
      scal_diff_mean[i,j] <- mean(scal_diff[[i,j]])
    }
  } 
}


##Plots


pdf("sweep_compare_meanmcempir_simulated.pdf")
lattice::levelplot(scal_diff_mean, xlab="mu",
                   ylab="kappa",col.regions = terrain.colors(100), row.values = mu,
                   column.values = kappa)
dev.off()



pdf("vcompare_meanmcempir_scaling.pdf")

colnames(scal_diff_mean) <- as.character(kappa)
rownames(scal_diff_mean) <- as.character(mu)
lattice::levelplot(scal_diff_mean, xlab="mu",
                   ylab="kappa",col.regions = terrain.colors(100))

dev.off()







