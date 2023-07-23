##load libraries

library(rgdal)
library(maptools)
library(rgeos)
library(colourvalues)
library(magick)
library(animation)
library(scales)

##Define diffusion function - with localised a
diffusion <- function(x, wm, a, cap=1){
  res <- rep(NA,length(x))
  for (i in 1:nrow(wm)){
    res[i] <- (x[i]*a[i]) + sum(wm[i,] * x)
  } 
  res[res > 1] <- cap
  return(res)
}

##define distance decay function

distfun.exp <- function(x,b){exp(-b*x)} #Distance decay function


##load data/parameter values
ttstep <- 14 # number of timesteps
nnodes <- 2000 # number of nodes
seeds <- c(1)
phlen <- c(175,258,142,163,49,138,38,72,145,55,35,57,128,133)
dur <- sum(phlen)
tstep <- dur
alpha <- 1.2
#kappa <- seq(from=0.5, to=0.7,by = 0.1)
kappa <- 0.6

##localised a is drawn from a lognormal distribution, with parameters mu, sigma,
## their values are set here

mean <- 1
sd <- 0.5 


directory <- getwd()
aegean <- readOGR(dsn=directory, layer="aeg") #vector map
utm35Nwgs84 <- CRS("+init=epsg:32635") #set crs
aegean <- spTransform(aegean, utm35Nwgs84) #project to specified crs


##set up origin node matrix

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


#test  
plot(aegean)
plot(all, add = T, col="blue", pch=19)

##START SIMULATION

meansmatrixlistlist <- vector(mode = "list", length = length(seeds))
mediansmatrixlistlist <- vector(mode = "list", length = length(seeds))
resultsllist <- vector(mode = "list", length = length(seeds))
nranl <- vector(mode = "list", length = length(seeds))
for (s in 1:length(seeds)){
  set.seed(seeds[s])
  nran0 <- spsample(aegean, n=nnodes, type="random")
  nran00 <- nran0[-1] #remove one point from nran0 so that when i add the origin i will still have a pop = nnodes
  nran <- rbind(nran00, all[1]) ##set to index of desired origin point from the 'all' vector
  nran$siteid <- row.names(nran) ##add siteid column to nran
  alpha_local <- rlnorm(nrow(nran), meanlog=mu, sdlog=sigma) ##lognormal
  ##alpha_local <- rnorm(nrow(nran), mean=1, sd=0.5) ##normal with mean = 1 and sd =0.5  
  distmat <- gDistance(nran, byid=TRUE)
  weights <- distmat/1000 #meters to km
  resultsl <- vector(mode = "list", length = length(kappa))
  weightl <- vector(mode = "list", length = length(kappa))
  for (l in 1:length(kappa)){
    weightl[[l]] <- apply(weights, 1, FUN=distfun.exp, b=kappa[l])
    diag(weightl[[l]]) <- 0 #so that nodes don't increase on their own
    
    bg <- which.max(nran$x) #find index of point with max longitude (easternmost)
    trait <- rep(c(0),each = nnodes)
    trait[bg] <- 0.1 # this sets the origin point trait proportion value
    ## trait[nnodes] <- 0.1  #alternatively, use one of the origin points in the 'all' vector 
    trait0 <- trait
    
    results <- matrix(NA, nrow = nnodes, ncol = tstep) # results matrix (size is determined by parameter values set at the beginning)
    results[,1] <- trait0
    for (p in 2:ncol(results)){
      trait <- diffusion(x = trait, wm = weightl[[l]], a=alpha_local)
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

## Plotting function
plot.fun.3.local <- function(Timestep,Nodes,Results, meantrend,upper, lower){ 
  Nodes$Value <- Results[,Timestep]
  Nodes$Colour  <- colourvalues::colour_values(Nodes$Value)
  par(mfrow=c(1,2))
  par(mar=c(1, 2, 1, 2))
  plot(aegean, col="white", border="black", lwd=0.1,
       main = paste("t = ", Timestep, ", N = ",nnodes,", k = ", kappa))
  points(Nodes, pch=19,col = alpha(Nodes$Colour, 0.6), cex=1)
  par(mar=c(6, 4, 6, 2))
  plot(meantrend, ylim = c(0,100), cex = 0.1, 
       ylab= "Mean trait proportion value", xlab="Time")
  lines(meantrend, col = "purple", lwd = 1)
  lines(upper, col = "grey", lwd = 1)
  lines(lower, col = "grey", lwd = 1)
  points(x = Timestep , y = meantrend[Timestep], col= "purple2", bg = "yellow", pch=21, ylim = c(0,100), cex=1.6)
  par(mfrow=c(1,1)) 
}  


##Simplify data names

resultaki <- resultsllist[[1]][[1]]
meanaki <- meansmatrixlistlist[[1]][[1]]
sdaki <- sdmatrixlistlist[[1]][[1]]
upper <- meanaki + sdaki
lower <- meanaki - sdaki

##Apply function to a single timestep

pdf("lognorm.pdf", width = 10, height = 7, paper = "USr")

t <- 400
plot.fun.3.local(Timestep = t,Nodes=nranl[[1]], Results=resultaki, meantrend = meanaki, upper = upper, lower = lower)

dev.off()


##map and timeseries, using plot.fun 3

saveGIF({
  for (i in 500:900) plot.fun.3(Timestep = i,Nodes=nranl[[1]], Results=resultaki, meantrend = meanaki)},
  movie.name = "localised_a_video.gif",interval = 0.1, ani.height = 500, ani.width = 500
)



## OPTIONAL plot 11/03/2022 geo

n7000g <- meansmatrixlistlist[[1]][[1]]

pop <- c(1000, 2000, 3000, 4000, 5000, 6000)
rbPal <- colorRampPalette(c('lightblue','blue'))
popcols <- rbPal(length(pop))[as.numeric(cut(pop,breaks = length(pop)))]

pdf("pop1000to6000geok0.6a1.2.pdf")

plot(n1000g, 
     ylim = c(0,100), cex = 0.1, ylab= "Mean adoption percentage", xlab="Time", sub=paste("kappa = ",kappa,", alpha = " , alpha))
lines(n1000g, col = popcols[1], lwd = 3)
lines(n2000g, col = popcols[2], lwd = 3)
lines(n3000g, col = popcols[3], lwd = 3)
lines(n4000g, col = popcols[4], lwd = 3)
lines(n5000g, col = popcols[5], lwd = 3)
lines(n6000g, col = popcols[6], lwd = 3)
legend("topleft",legend=as.character(pop),fill = popcols, col=popcols, cex=0.6, title="N")


dev.off()
