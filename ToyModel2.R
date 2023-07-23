##Toy Model 2 -  Ring networks!

##load libraries
library(rgdal)
library(maptools)
library(rgeos)
library(colourvalues)
library(magick)
library(animation)
library(scales)

##Parameter values

kappa <- 0.5
alpha <- 1.2
origin <- 0.1 
seeds <- c(1)
phlen <- c(175,258,142,163,49,138,38,72,145,55,35,57,128,133)
dur <- sum(phlen)
tstep <- 100
rad <- 10 ##radius of circle
div <- 10
loss <- 1 ##trait loss allowed when loss = 1, set loss = 0 for no loss


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

##define function to get coordinates of points along a circle! Arguments: the angle theta and the radius of the circle r

circle.point.fun <- function(theta, radius){
  x <- radius * cos(theta)
  y <- radius * sin(theta)
  coord <- c(x,y)
  return(coord)
}

##set up the points and network

#angles <- seq(from = 0, to = 1080, by = 120) #n = 10
angles <- seq(from = 0, to = 1140, by = 60) #n = 20
pointlist <- vector(mode="list", length = length(angles))
for (i in 1:length(angles)){
  pointlist[[i]] <- circle.point.fun(theta = angles[i], radius =rad)
}
points <- do.call("rbind", pointlist)
plot(points)

#case with two nodes
#n1b <- c(1,1)
#n2b <- c(2,1)
#points <- rbind(n1b, n2b)
#angles <- c(NA, NA)
#plot(points, asp = 1)

meansmatrixlistlist <- vector(mode = "list", length = length(seeds))
mediansmatrixlistlist <- vector(mode = "list", length = length(seeds))
sdmatrixlistlist <- vector(mode = "list", length = length(seeds))
resultsllist <- vector(mode = "list", length = length(seeds))
nranl <- vector(mode = "list", length = length(seeds))
for (s in 1:length(seeds)){
  set.seed(seeds[s])
  nran <- points
  distmat <- dist(nran, method = "euclidean", diag = T, upper = T)
  weights <- as.matrix(distmat)
  resultsl <- vector(mode = "list", length = length(kappa))
  weightl <- vector(mode = "list", length = length(kappa))
  for (l in 1:length(kappa)){
    weightl[[l]] <- apply(weights, 1, FUN=distfun.exp, b=kappa[l])
    diag(weightl[[l]]) <- 0 #so that nodes don't increase on their own
    trait <- rep(c(0),each = nrow(points))
    trait[1] <- 0.01 # this sets the origin point trait proportion value
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
    sdmatrixlist[[b]] <- apply(resultsl[[b]], 2, sd) * 100
  }  
  
  
  meansmatrixlistlist[[s]] <- meansmatrixlist
  sdmatrixlistlist[[s]] <- sdmatrixlist
  resultsllist[[s]] <- resultsl
  nranl[[s]] <- nran
}

##edge list

empt <- vector(mode = "list", length = length(angles))
for (z in 1:length(angles)){
  empt[[z]] <- t(replicate(length(angles), nranl[[1]][z,]))
  empt[[z]] <- cbind(empt[[z]], points)
}
empt1 <- do.call("rbind", empt)

##plot 

#pdf("toy2environment.pdf", width = 10, height = 7, paper = "USr")


plot(points, asp = 1, pch = 19, col = alpha("lightblue", 0.8), cex = 3,
     xaxt = "n", yaxt = "n", ylab = "", xlab = "")
for (u in 1:nrow(empt1)){
  segments(x0 = empt1[u,1], y0 = empt1[u,2], x1 = empt1[u,3], y1 = empt1[u,4], col = alpha("grey", 0.1)) 
}


#dev.off()


##Plotting function
plot.fun.toy2 <- function(Timestep,Nodes,Results, meantrend, edges, upper, lower){
  Nodes <- as.data.frame(Nodes)
  colnames(Nodes) <- c("x", "y")
  Nodes$Value <- Results[,Timestep]
  Nodes$Colour  <- colourvalues::colour_values(Nodes$Value)
  par(mfrow=c(1,2))
  par(mar=c(2, 2, 2, 2))
  plot(x = Nodes$x , y = Nodes$y, col = alpha(Nodes$Colour, 0.6), pch = 19,
       asp = 1, cex = 3,xaxt = "n", yaxt = "n",
       main = paste("t = ", Timestep, ", N = ",length(angles),", a = " , alpha, ", k = ",
                    kappa[1]))
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



##test, apply function to a single timestep

pdf("toy2.pdf", width = 10, height = 7, paper = "USr")

t <- 5
plot.fun.toy2(Timestep = t,Nodes=nranl[[1]], Results=resultaki, meantrend = meanaki, edges = empt1,
              upper = upper, lower = lower)

dev.off()

##map and timeseries

saveGIF({
  for (i in 1:tstep) plot.fun.toy2(Timestep = i,Nodes=nranl[[1]], Results=resultaki, meantrend = meanaki, edges = empt1)},
  movie.name = "15082022_toy2k1a1origin0.000000001.gif",interval = 0.1, ani.height = 500, ani.width = 500
)

