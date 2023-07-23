##This script includes a function to compare simulation to empirical data and sets up a workflow to 
##sweep parameter values and produce a parameter space map.

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


##LOAD DATA


##load simulation data/parameter values
alpha <- 1
kappa <- 0.4
origin <- 0.01
tstep <- 100
n <- 200

##coordinates for the four points delineating the square 
xmn <- 0
xmx <- 100
ymn <- 0
ymx <- 100


##set up the points and network

set.seed(30)

##define toy environment (a square with a side = 10)
square <- rbind(c(xmn,ymn), c(xmx, ymn),c(xmx, ymx), c(xmn, ymx), c(xmn,ymn))
#plot(square, asp = 1)
#lines(square)
 
##random nodes on toy environment
pointx <- runif(n, min = xmn, max = xmx) #n 
pointy <- runif(n, min = ymn, max = ymx) #n 
points <- matrix(NA, nrow = n, ncol = 2)
points[,1] <- pointx
points[,2] <- pointy


##plot nodes in environment
#plot(points, asp = 1)


##SIMULATION - Double nested loop, outer loop sweeps alpha, inner loop sweeps kappa
weightmatl <- vector(mode = "list", length = length(alpha))
meansmatrixlistlist <- vector(mode = "list", length = length(alpha))
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
    trait[1] <- origin # this sets the origin point trait proportion value
    trait0 <- trait
    
    results <- matrix(NA, nrow = n, ncol = tstep) # results matrix (size is determined by parameter values set at the beginning)
    results[,1] <- trait0
    for (p in 2:ncol(results)){
      trait <- diffusion.nloss(x = trait, wm = weightl[[l]], a=alpha[a])
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
  
  meansmatrixlistlist[[a]] <- meansmatrixlist
  resultsllist[[a]] <- resultsl
  nranl[[a]] <- nran
  weightmatl[[a]] <- weightl
}


##PLOT SYSTEM

##edge list

empt <- vector(mode = "list", length = n)
for (z in 1:n){
  empt[[z]] <- t(replicate(n, nranl[[1]][z,]))
  empt[[z]] <- cbind(empt[[z]], points)
}
empt1 <- do.call("rbind", empt)

##PLOTTING FUNCTION DEFINITION BEGINS HERE
plot.fun.toy4 <- function(Timestep,Nodes,Results, meantrend, edges, aind, kind){ #aind is index of alpha, kind is index of kappa
  Nodes <- as.data.frame(Nodes)
  colnames(Nodes) <- c("x", "y")
  Nodes$Value <- Results[,Timestep]
  Nodes$Colour  <- colourvalues::colour_values(Nodes$Value)
  par(mfrow=c(1,2))
  par(mar=c(2, 2, 2, 2))
  plot(x = Nodes$x , y = Nodes$y, col = alpha(Nodes$Colour, 0.6), pch = 19,
       asp = 1, cex = 2,xaxt = "n", yaxt = "n",
       main = paste("t = ", Timestep, ", N = ",n,", a = " , alpha[aind], ", k = ",
                    kappa[kind]))
  par(mar=c(6, 2, 6, 2))
  plot(meantrend, ylim = c(0,100), cex = 0.1, 
       ylab= "Mean trait proportion value", xlab="Time")
  lines(meantrend, col = "blue", lwd = 3)
  points(x = Timestep , y = meantrend[Timestep], col= "red", pch=19, ylim = c(0,100), cex=1.6)
  par(mfrow=c(1,1)) 
}  


##Data
resultaki <- resultsllist[[1]][[1]]
meanaki <- meansmatrixlistlist[[1]][[1]]


##calculate betweenness centrality

##try with the sna package

library(sna)

weightz <- weightmatl[[1]][[1]]
betw <- sna::betweenness(weightz,ignore.eval=FALSE) ##second argument includes weights in the calculation of betw

##check correlation between values 
t <- 30
plot(x = betw, y = log(resultaki[,t]))

##pdf("toy5n500k0.4a1.2origin0.01tstep10.pdf", width = 10, height = 7, paper = "USr")

plot.fun.toy4(Timestep = t,Nodes=nranl[[1]], Results=resultaki, meantrend = meanaki, edges = empt1, aind = 1, kind = 1)

##dev.off()

##map and timeseries

saveGIF({
  for (i in 1:tstep) plot.fun.toy4(Timestep = i,Nodes=nranl[[1]], Results=resultaki, meantrend = meanaki, edges = empt1, aind = 1,
                                   kind = 1)},
  movie.name = "13012023_toy5_a1k0.8n500.gif",interval = 0.1, ani.height = 500, ani.width = 500
)


##Load results for qualitative sweep for N, kappa, alpha (separate panels for each N (200, 300, 400, 500, 600, 700), seed 30
## variable recorded was presence, somewhat presence, or absence of cluster-driven diffusion dummy coded as 0, 0.5, 1, NA 
## where NA stands for unsuccessful diffusion beyond the origin cluster

N200 <- c(1, NA, NA, NA, NA, 1, 1, NA, NA, NA, 0.5, 1, 1, NA, NA, 0.5, 1, 1, 1, 1, 0.5, 1, 1, 1, 1)
N200 <- matrix(N200, nrow = length(kappa), ncol = length(alpha))
N200r <- raster(N200)


N300 <- c(0.5, 1, 1, 1, NA, 0, 1, 0.5, 1, 1, 0, 1, 1, 1, 1, 0, 0.5, 1, 1, 0.5, 0, 0, 0.5, 0.5, 0.5)
N300 <- matrix(N300, nrow = length(kappa), ncol = length(alpha))
N300r <- raster(N300)

N400 <- c(0,1,1,NA,NA,0, 0.5, 1, 1, 1, 0,0, 1, 1, 1, 0,0,0,0.5,1,0,0,0,0.5,0.5)
N400 <- matrix(N400, nrow = length(kappa), ncol = length(alpha))
N400r <- raster(N400)


N500 <- c(0,0.5,1,1,NA,0,0.5,1,1,1,0,0,1,1,1,0,0,0.5,1,1,0,0,0,0,1)
N500 <- matrix(N500, nrow = length(kappa), ncol = length(alpha))
N500r <- raster(N500)


N600 <- c(0,0.5,1,1,1,0,0,1,1,1,0,0,0.5,1,1,0,0,0,0.5,0.5,0,0,0,0,0.5)
N600 <- matrix(N600, nrow = length(kappa), ncol = length(alpha))
N600r <- raster(N600)


N700 <- c(0,0,1,1,1,0,0,1,1,1,0,0,0.5,1,1,0,0,0,0.5,0.5,0,0,0,0,0)
N700 <- matrix(N700, nrow = length(kappa), ncol = length(alpha))
N700r <- raster(N700)

##PLOTS

pdf("toymodel4_qualitative_seed30.pdf")

plotlist <- c(N200r, N300r, N400r, N500r , N600r, N700r)
Ns <-seq(from = 200, to = 700, by = 100)
par(mfrow=c(2,3))
for (i in 1:length(plotlist)){
  plot(plotlist[[i]], xaxt = "n", yaxt = "n", xlab = "alpha", ylab = "kappa", main = paste("N =",Ns[i]))
  axis(1, at = seq(0.1,0.9, by = 0.2),labels = alpha)
  axis(2, at = rev(seq(0.1,0.9, by = 0.2)),labels = kappa)
}

dev.off()

##plot param space map for scaling param for all N

pdf("sweep14092022_toymodel4_nodes_scaling_seed30.pdf")

N_scal_list <- list(scal_matN200, scal_matN300, scal_matN400, scal_matN500, scal_matN600, scal_matN700)
N_scal_list2 <- lapply(N_scal_list, FUN = raster)
Ns <-seq(from = 200, to = 700, by = 100)
par(mfrow=c(2,3))
for (i in 1:length(N_scal_list2)){
  plot(N_scal_list2[[i]],main = paste("N =",Ns[i]), 
       xlab = "kappa", ylab = "alpha", xaxt = "n", yaxt = "n")
  axis(1, at = seq(0.1,0.9, by = 0.2),labels = kappa)
  axis(2, at = rev(seq(0.1,0.9, by = 0.2)),labels = alpha)
}

dev.off()

