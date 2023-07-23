##This script loads and processes wheel pottery data
##It also includes a function to compute confidence intervals for the mean

##libraries
library("readxl")
library("scales")

##seed for random number generation
set.seed(13)


##load data
wheel_data <- read_excel("wheel_percentages_may2022_tidy.xlsx")
wheel_data0 <- wheel_data[,c(2:15)] #select only percentages per phase
wheel_data1 <- as.matrix(wheel_data0)

phlen <- c(175,258,142,163,49,138,38,72,145,55,35,57,128,133)
phend <- c() ##this object has labels for the year (BC) that each phase ends on
phend[1] <- 2475
for (i in 2:length(phlen)){
  phend[i] <- phend[i-1] - phlen[i]
}
  
  
  
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


##plot
##mean
plot(tseries, type = "n", ylim = c(0,100))
lines(tseries)

##median
plot(tseriesmed, type = "n", ylim = c(0,100))
lines(tseriesmed, lwd=3)

##Calculate confidence intervals using the method suggested by Shennan in Quantifying Archaeology (1997:83):
## 1. Use Student's t-distribution for a given number of degrees of freedom instead of the 
## Z-score associated with a normal distribution, as the sample size per phase is smaller than 30.
## 2. Formula: mean +/-t value for the number of degrees of freedom I want, multiplied by
## the ratio of the standard deviation to the square root of the sample size


conf.int.fun <- function(samplemeans, stdvs, data){
  samplesize0 <- colSums(!is.na(data)) ##estimates the sample size for all timesteps (without NAs)
  confints <- list()
  for (i in 1:ncol(data)){
    timestep <- i
    samplesize <- samplesize0[timestep]
    tvalue <- qt(.95,df = samplesize-1) ## appropriate t value for 95 % confidence level
    conf <- tvalue*(stdvs[timestep]/sqrt(samplesize))
    confints[[i]] <- c(samplemeans[timestep] - conf, samplemeans[timestep] + conf) 
  }
  confints2 <- list(do.call(rbind, confints)[,1], do.call(rbind, confints)[,2])
  return(confints2)  
}

##apply function to wheel data and plot the mean over time with a 95 % confidence envelope

wheel <- conf.int.fun(samplemeans = tseries, stdvs = tseries_sd, data = vou)

pdf("wheel_timeseries_CI.pdf")
plot(tseries, type = "n", xlab = "Time (Phases)", ylab = "Mean percentage of wheel-made pottery", ylim = c(0,100))
polygon(c(1:14, rev(1:14)), c(wheel[[1]], rev(wheel[[2]])), density = 20, col = 'grey', border = NA)
lines(tseries, type = "l", ylim = c(0, 100), lwd = 3)
lines(wheel[[1]], type = "l", col = "grey", lwd = 2)
lines(wheel[[2]], type = "l", col = "grey", lwd = 2)

dev.off()

##create line graph that respects the non-linearity of the time x axis
##How: use as y the tseries or simulated trend values, and as x the following vector:
## this vector takes the phlen vector with the phase lengths and adds a 0 to it:
## 0 175 258 142 163  49 138  38  72 145  55  35  57 128 133
## write a loop that for each element (n = 15) adds an element to a new vector
## such that the value at iteration i is equal to the sum of the previous ones
vec <- c(0,phlen)
vec2 <- c()
for (i in 1:length(vec)){
  vec2[i] <- sum(vec[1:i])
}
vec2 <- vec2[2:15]

##plot with non-linear x axis (time), where the width of each bin corresponds to the length of each phase

pdf("wheel_timeseries_phaselength_CI.pdf")

plot(x = vec2, y = tseries, type = "n",xaxt = 'n', xlab = "Time (Years BC)", ylab = "Mean percentage of wheel-made pottery", ylim = c(0,100))
polygon(c(vec2, rev(vec2)), c(wheel[[1]], rev(wheel[[2]])), density = 20, col = 'grey', border = NA)
lines(x = vec2, y = tseries, type = "l", ylim = c(0, 100), lwd = 3)
lines(x = vec2, y = wheel[[1]], type = "l", col = "grey", lwd = 2)
lines(x = vec2, y = wheel[[2]], type = "l", col = "grey", lwd = 2)
axis(side=1,at=vec2,labels=paste(phend), cex.axis=1)


dev.off()

## Confidence envelope via Monte Carlo simulation
##1. Sample the x coordinate from the time interval keeping the phase mean as the y, temporal uncertainty.
##2. Sample the y coordinate from the confidence interval multiple times and keep x s the phase end year.
## i.e. sample size uncertainty
##3. Sample both the x and the y coordinate as above, make it 1000 * 1000 simulations.

## 1: temporal uncertainty

vec3 <- vec[2:15]
mctemp <- vector(mode = "list", length = 1000)
for (i in 1:length(mctemp)){
  tser <- vector(length = 14)
  for (j in 1:length(tser)){
    tser[j] <- round(runif(1, min = phend[j], max = phend[j] + vec3[j]))
  }
  mctemp[[i]] <- tser
}

##plot 

pdf("mc_temp.pdf")

plot(x = rev(mctemp[[1]]),type = "n", y = tseries, xlab = "Time (Years BC)", 
     ylab = "Mean percentage of wheel-made pottery", ylim = c(0,100),
     main = "1000 simulations reflecting temporal uncertainty", xaxt ="n")
for (i in 1:length(mctemp)){
  lines(x = rev(mctemp[[i]]), y = tseries, type = "l", ylim = c(0, 100), lwd = 0.1, col = alpha("lightblue", 0.2)) 
}
axis(side=1,at=mctemp[[1]],labels=rev(mctemp[[1]]), cex.axis=0.9)


dev.off()


## 2 : sample size uncertainty

mcsamp <- vector(mode = "list", length = 1000)
for (i in 1:length(mcsamp)){
  ciser <- vector(length = 14)
  for (j in 1:length(ciser)){
    mcwheel <- conf.int.fun(samplemeans = tseries, stdvs = tseries_sd, data = vou)
    ciser[j] <- runif(1, min = mcwheel[[1]][j], max = mcwheel[[2]][j])
  }
  mcsamp[[i]] <- ciser
}

##plot 

pdf("mc_samp.pdf")

plot(x = vec2, y = tseries, xlab = "Time (Years BC)", 
     ylab = "Mean percentage of wheel-made pottery", ylim = c(0,100),
     main = "1000 simulations reflecting sample size uncertainty",
     type = "n",xaxt = 'n')
for (i in 1:length(mcsamp)){
  lines(x = vec2, y = mcsamp[[i]], type = "l", ylim = c(0, 100), lwd = 0.1, col = alpha("lightblue", 0.2)) 
}
axis(side=1,at=vec2,labels=paste(phend), cex.axis=1)

dev.off()


## 3 : temporal and sample size uncertainty ##NEEDS DEBUG

mcsamptemp <- vector(mode = "list", length = 1000)
for (i in 1:length(mcsamptemp)){
  citser <- vector(mode = "list", length = 14)
  for (j in 1:length(citser)){
    mcwheel <- conf.int.fun(samplemeans = tseries, stdvs = tseries_sd, data = vou)
    citser[[j]] <- c(runif(1, min = mcwheel[[1]][j], max = mcwheel[[2]][j]),
                   round(runif(1, min = phend[j], max = phend[j] + vec3[j])))
  }
  mcsamptemp[[i]] <- do.call("rbind",citser)
}


##plot 

pdf("mc_samptemp.pdf")

plot(x = rev(mcsamptemp[[1]][,2]), y = mcsamptemp[[1]][,1], xlab = "Time (Years BC)", 
     ylab = "Mean percentage of wheel-made pottery", ylim = c(0,100),
     main = "1000 simulations reflecting temporal and sample size uncertainty",
     type = "n", xaxt = "n", xlim = c(1000, 2600))
for (i in 1:length(mcsamptemp)){
  lines(x = rev(mcsamptemp[[i]][,2]), y = mcsamptemp[[i]][,1], ylim = c(0, 100), lwd = 0.1,
        col = alpha("lightblue", 0.2)) 
}
axis(side=1,at=seq.int(from = 1000, to = 2600, by = 200),labels=rev(seq.int(from = 1000, to = 2600, by = 200)), cex.axis=0.9)

dev.off()
