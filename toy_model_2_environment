##load libraries

library(rgdal)
library(maptools)
library(rgeos)
library(colourvalues)
library(magick)
library(animation)
library(scales)
library(readxl)
library(lattice)
library(raster)

##Function definitions

##define diffusion function
diffusion <- function(x, wm, a, cap=1){
  res <- rep(NA,length(x))
  for (i in 1:nrow(wm)){
    res[i] <- (x[i]*a) + sum(wm[i,] * x)
  } 
  res[res > 1] <- cap
  return(res)
}

##define diffusion function with random trait loss
diffusion.nloss <- function(x, wm, a, cap=1){
  res <- rep(NA,length(x))
  for (i in 1:nrow(wm)){
    res[i] <- (x[i]*a) + sum(wm[i,] * x)
    res[i] <- res[i] - runif(1, min = 0, max = res[i]) #trait loss
  } 
  res[res > 1] <- cap
  return(res)
}


##define distance decay function

distfun.exp <- function(x,b){exp(-b*x)} #Distance decay function

##define function to get coordinates of points along a circle. Arguments: the angle theta and the radius of the circle r

circle.point.fun <- function(theta, radius){
  x <- radius * cos(theta)
  y <- radius * sin(theta)
  coord <- c(x,y)
  return(coord)
}


##LOAD DATA


##load simulation data/parameter values
ttstep <- 14 # number of timesteps
phlen <- c(175,258,142,163,49,138,38,72,145,55,35,57,128,133)
dur <- sum(phlen)
tstep <- dur
alpha <- c(1)
kappa <- c(0.5)
origin <- 0.1
tstep <- 100
rad <- 10 #radius of circle
div <- 10
n <- 20
loss <- 1 #trait loss allowed when loss = 1, set loss = 0 for no loss

##set up the points and network

set.seed(1) ##seed for random number generation

angles <- round(runif(n, min = 0, max = 360)) #n = 20 nodes randomly placed on a circle
pointlist <- vector(mode="list", length = length(angles))
for (i in 1:length(angles)){
  pointlist[[i]] <- circle.point.fun(theta = angles[i], radius =rad)
}
points <- do.call("rbind", pointlist)
plot(points)

##1000 simulations

##simulation_list <- vector(mode = "list", length = 1000)
##for (w in 1:1000){

##SIMULATION - Double nested loop, outer loop sweeps alpha, inner loop sweeps kappa

weightmatl <- vector(mode = "list", length = length(alpha))
meansmatrixlistlist <- vector(mode = "list", length = length(alpha))
mediansmatrixlistlist <- vector(mode = "list", length = length(alpha))
sdmatrixlistlist <- vector(mode = "list", length = length(alpha))
resultsllist <- vector(mode = "list", length = length(alpha))
nranl <- vector(mode = "list", length = length(alpha))
for (a in 1:length(alpha)){
  nran <- points
  distmat <- dist(nran, method = "euclidean", diag = T, upper = T)
  weights <- as.matrix(distmat)
  resultsl <- vector(mode = "list", length = length(kappa))
  weightl <- vector(mode = "list", length = length(kappa))
  for (l in 1:length(kappa)){
    weightl[[l]] <- apply(weights, 1, FUN=distfun.exp, b=kappa[l])
    diag(weightl[[l]]) <- 0 #so that nodes don't increase on their own
    trait <- rep(c(0),each = n)
    trait[1] <- 0.1 # this sets the origin point trait proportion value
    trait0 <- trait
    
    results <- matrix(NA, nrow = n, ncol = tstep) # results matrix (size is determined by parameter values set at the beginning)
    results[,1] <- trait0
    for (p in 2:ncol(results)){
      if(loss == 0){
        trait <- diffusion(x = trait, wm = weightl[[l]], a=alpha[a])
      } else {
        trait <- diffusion.nloss(x = trait, wm = weightl[[l]], a=alpha[a])
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
    sdmatrixlist[[b]] <- apply(resultsl[[b]], 2, sd)
  }  
  
  meansmatrixlistlist[[a]] <- meansmatrixlist
  sdmatrixlistlist[[a]] <- sdmatrixlist
  resultsllist[[a]] <- resultsl
  nranl[[a]] <- nran
  weightmatl[[a]] <- weightl
}

##simulation_list[[w]] <- meansmatrixlistlist
##}

##save results of 1000 simulations
##capture.output(simulation_list, file = "toy2.csv")


##edge list

empt <- vector(mode = "list", length = length(angles))
for (z in 1:length(angles)){
  empt[[z]] <- t(replicate(length(angles), nranl[[1]][z,]))
  empt[[z]] <- cbind(empt[[z]], points)
}
empt1 <- do.call("rbind", empt)

##Plotting function
plot.fun.toy3 <- function(Timestep,Nodes,Results, meantrend, edges, aind, kind,
                          upper, lower){ #aind is index of alpha, kind is index of kappa
  Nodes <- as.data.frame(Nodes)
  colnames(Nodes) <- c("x", "y")
  Nodes$Value <- Results[,Timestep]
  Nodes$Colour  <- colourvalues::colour_values(Nodes$Value)
  par(mfrow=c(1,2))
  par(mar=c(2, 2, 2, 2))
  plot(x = Nodes$x , y = Nodes$y, col = alpha(Nodes$Colour, 0.6), pch = 19,
       asp = 1, cex = 3,xaxt = "n", yaxt = "n",
       main = paste("t = ", Timestep, ", N = ",length(angles),", a = " , alpha[aind], ", k = ",
                    kappa[kind]))
  for (u in 1:nrow(edges)){
    segments(x0 = edges[u,1], y0 = edges[u,2], x1 = edges[u,3], y1 = edges[u,4], col = alpha("grey", 0.1)) 
  }
  par(mar=c(6, 4, 6, 2))
  plot(meantrend, ylim = c(0,100), cex = 0.1, 
       ylab= "Mean trait proportion value", xlab="Time")
  lines(meantrend, col = "purple", lwd = 3)
  lines(upper, col = "grey", lwd = 1)
  lines(lower, col = "grey", lwd = 1)
  points(x = Timestep , y = meantrend[Timestep], col= "purple2", bg = "yellow", pch=21, ylim = c(0,100), cex=1.6)
  par(mfrow=c(1,1)) 
}  


##Data
resultaki <- resultsllist[[1]][[1]]
meanaki <- meansmatrixlistlist[[1]][[1]]
sdaki <- sdmatrixlistlist[[1]][[1]]
upper <- meanaki + sdaki
lower <- meanaki - sdaki


##calculate node degree with the sna package

library(sna)

weightz <- weightmatl[[1]][[1]]
deg <- sna::degree(weightz,ignore.eval=FALSE) ##second argument includes weights in the calculation of betw

##check correlation between strength and trait values 

t <- tstep
correl <- cor.test(resultaki[,t],deg, method = 'pearson')

pdf("degree_correlation.pdf", width = 10, height = 7, paper = "USr")

par(mar=c(6, 6, 4, 4))
plot(x = deg, y = resultaki[,t],
     ylab = "Trait proportion value", xlab = "Node Strength", 
     pch = 19, col = alpha("grey", 0.5), main = paste('r = ',correl$estimate))

dev.off()


##test, apply function to a single timestep

##pdf("timestep.pdf", width = 10, height = 7, paper = "USr")

plot.fun.toy3(Timestep = t,Nodes=nranl[[1]], Results=resultaki, meantrend = meanaki, edges = empt1, aind = 1, kind = 1,
              upper = upper, lower = lower)

##dev.off()

##map and timeseries

saveGIF({
  for (i in 1:tstep) plot.fun.toy3(Timestep = i,Nodes=nranl[[1]], Results=resultaki, meantrend = meanaki, edges = empt1, aind = 2,
                                   kind = 2,
                                   upper = upper, lower = lower)},
  movie.name = "simulation_video.gif",interval = 0.1, ani.height = 500, ani.width = 500
)


##plot 1000 simulations

pdf("toy2_1000simulations.pdf", width = 10, height = 7, paper = "USr")

plot(simulation_list[[1]][[1]][[1]], ylim = c(0,100), cex = 0.1, 
      ylab= "Mean trait proportion value", xlab="Time")
for (i in 1:1000){
  lines(simulation_list[[i]][[1]][[1]], col = alpha("grey", 0.2), lwd = 2)
}

dev.off()





