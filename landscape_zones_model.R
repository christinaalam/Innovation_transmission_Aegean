##load libraries


library(rgdal)
library(maptools)
library(rgeos)
library(colourvalues)
library(magick)
library(animation)
library(scales)
library(raster)
library(raptr)
library(dismo)
library(sf)
library(terra)

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


##load data/parameter values

nnodes <- 2000 # number of nodes
seeds <- c(1)
phlen <- c(175,258,142,163,49,138,38,72,145,55,35,57,128,133)
dur <- sum(phlen)
tstep <- dur
alpha <- 1.2
kappa <- 0.2
#kappa <- seq(from=0.5, to=0.7,by = 0.1) #alternatively, if you want multiple kappa values
loss <- 1 #trait loss allowed when loss = 1, set loss = 0 for no loss


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


#plot origins
plotcolors <- c("darkblue", "cyan4","darkseagreen", "darkred", "darkmagenta", "forestgreen", 
                "deeppink2", "black", "darkgoldenrod3" )


pdf("origins.pdf", width = 10, height = 7, paper = "USr")
plot(aegean)
for (i in 1:length(all)){
  plot(all[i], add = T, col=plotcolors[i], pch=19, cex = 2)
}
dev.off()

##Define thresholds and probabilities for the landscape classification

coastthreshold <- 5000 #the distance threshold defining coastal regions <= 5000 metres from the coast
islandthreshold <- 10000 #the size of islands where a point will be force-plotted if they are not randomly assigned a point
slopethreshold <- 10
coaststeepprob <- 0.05 # the probability that points will be plotted on coastal land
coastflatprob <- 0.70
inlandflatprob <- 0.20 # the probability that points will be plotted on inland land
inlandsteepprob <- 1-(coastflatprob + coaststeepprob + inlandflatprob)


##Load spatial data


proxcoast.new <- raster("distance_to_coast_utm35n.tif") #proximity to coast raster
aegeandem.new <- raster("aegean_dem_utm35n.tif") #aegean dem raster
directory <- getwd()
aegean <- readOGR(dsn=directory, layer="aeg") #vector map
utm35Nwgs84 <- CRS("+init=epsg:32635") #set crs
aegean <- spTransform(aegean, utm35Nwgs84) #project to specified crs
proxcoast.new <- raster::crop(proxcoast.new, aegean) #crop proxcoast.new to the extent of the aegean vector map
aegeandem.new <- raster::crop(aegeandem.new, aegean) #crop proxcoast.new to the extent of the aegean vector map
aegeandem.agg <- raster::aggregate(aegeandem.new, fact=5)


## Reclassify aegeandem.agg to get a second layer where NAs are 1 and all other values are NA

values <- getValues(aegeandem.agg)
values <- na.omit(values)
min0 <- min(values) ##upper limit for elev values
max0 <- max(values) ##upper limit for elev values
rclmat <- matrix(c(min0, max0, NA, NA, NA, 1), ncol=3, byrow=TRUE)
rclmat
masklayer <- reclassify(aegeandem.agg, rclmat)

##Get new distance to coast layer using the new DEM 
coastdist <- terra::distance(masklayer)

##Get slope layer using the new DEM 
aegeanslope <- raster::terrain(aegeandem.agg, opt="slope", unit="degrees")

##Take the distance to coast layer and classify it as a categorical raster (e.g. with a value 1 and 5, where 1 is coast and 5 is inland

values <- getValues(coastdist)
values <- na.omit(values)
maxval <- max(values) ##upper limit for proximity values
rclmat <- matrix(c(0, coastthreshold, 1, coastthreshold, maxval, 5), ncol=3, byrow=TRUE)
rclmat
coastdist1to5 <- reclassify(coastdist, rclmat)

##Take the slope layer and classify it as a categorical raster (e.g. with a value 3 and 4, where 3 is flat and 4 is steep)

values <- getValues(aegeanslope)
values <- na.omit(values)
maxval <- max(values) ##upper limit for slope values
minval <- min(values)
rclmat <- matrix(c(minval-1, slopethreshold, 3, slopethreshold, maxval+1, 4), ncol=3, byrow=TRUE)
rclmat
aegeanslope3to4 <- reclassify(aegeanslope, rclmat)
plot(aegeanslope3to4)

##Add coastdist1to5 (the categorical distance to coast) and aegeanslope3to4 (the categorical slope) using the raster:overlay() function.

#ddrast <- raster::overlay(coastdist1to5,aegeanslope3to4, fun=function(x,y){return(x+y)})

##or just add them

addrast <- coastdist1to5 + aegeanslope3to4
rclmat <- matrix(c(3, 0, 9, 4, 8, 3, 5, 2, 4, 1),
                 ncol=2, byrow=TRUE)
zones <- reclassify(addrast, rclmat) #in the resulting layer, coastal is 1, inland flat is 2, and inland steep is 3

## Then use a final reclassify to assign each zone the probabilities I want

rclmat <- matrix(c(1, coastflatprob, 2, coaststeepprob, 3, inlandflatprob, 4, inlandsteepprob), ncol=2, byrow=TRUE)
rclmat
proxclasses <- reclassify(zones, rclmat)

pdf("landscapezones1.pdf", width = 10, height = 7, paper = "USr")
plot(proxclasses)
dev.off()


##START SIMULATION
set.seed(seeds[1]) ##seed for random number generation


##10 simulations

##nsims <- 10
##simulation_list <- vector(mode = "list", length = 10)
##for (w in 1:nsims){

meansmatrixlistlist <- vector(mode = "list", length = length(seeds))
mediansmatrixlistlist <- vector(mode = "list", length = length(seeds))
sdmatrixlistlist <- vector(mode = "list", length = length(seeds))
resultsllist <- vector(mode = "list", length = length(seeds))
nranl <- vector(mode = "list", length = length(seeds))
for (s in 1:length(seeds)){
  nran0 <- dismo::randomPoints(proxclasses, n=nnodes, prob=TRUE) 
  nran <- nran0[1: nnodes,] # force exact sample size
  nran <- SpatialPoints(coords = nran,proj4string=CRS("+init=epsg:32635"))
  nran00 <- nran[-1] #remove one point from nran0 so that when i add the origin i will still have a pop = nnodes
  nran <- rbind(nran00, all[7]) ##add index of origin point of interest from the 'all' vector
  nran$siteid <- row.names(nran) ##add siteid column to nran
  print("ok")
  distmat <- gDistance(nran, byid=TRUE)
  weights <- distmat/1000 #meters to km
  resultsl <- vector(mode = "list", length = length(kappa))
  weightl <- vector(mode = "list", length = length(kappa))
  for (l in 1:length(kappa)){
    weightl[[l]] <- apply(weights, 1, FUN=distfun.exp, b=kappa[l])
    diag(weightl[[l]]) <- 0 #so that nodes don't increase on their own
    
    bg <- which.max(nran$x) #find index of point with max longitude (easternmost)
    trait <- rep(c(0),each = nnodes)
    trait[nnodes] <- 0.1 # this sets the origin point trait proportion value
    trait0 <- trait
    
    results <- matrix(NA, nrow = nnodes, ncol = tstep) # results matrix (size is determined by parameter values set at the beginning)
    results[,1] <- trait0
    for (p in 2:ncol(results)){
      if(loss == 0){
        trait <- diffusion(x = trait, wm = weightl[[l]], a=alpha)
      } else {
        trait <- diffusion.nloss(x = trait, wm = weightl[[l]], a=alpha)
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

##simulation_list[[w]] <- meansmatrixlistlist
##}

##END SIMULATION


##plot 10 simulations

pdf("loss_7origin_a1.2k0.2n2000seed1_landscape2simulations10.pdf", width = 10, height = 7, paper = "USr")

plot(simulation_list[[1]][[1]][[1]], ylim = c(0,100), cex = 0.1, 
     ylab= "Mean trait proportion value", xlab="Time", col = alpha("grey", 0.2))
for (i in 1:nsims){
  lines(simulation_list[[i]][[1]][[1]], col = alpha("grey", 0.2), lwd = 2)
}

dev.off()

##Simulation Data aggregation


## Plotting function

plot.fun.realistic.loss <- function(Timestep,Nodes,Results, meantrend, upper, lower){
  Nodes$Value <- Results[,Timestep]
  Nodes$Colour  <- colourvalues::colour_values(Nodes$Value)
  par(mfrow=c(1,2))
  par(mar=c(1, 2, 1, 2))
  plot(aegean, col="white", border="black", lwd=0.1,
       main = paste("t = ", Timestep, ", N = ",nnodes,", a = " , alpha, ", k = ", kappa))
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


##Rename data
resultaki <- resultsllist[[1]][[1]]
meanaki <- meansmatrixlistlist[[1]][[1]]
sdaki <- sdmatrixlistlist[[1]][[1]]
upper <- meanaki + sdaki
lower <- meanaki - sdaki



##test, apply function to a single timestep

pdf("t1000.pdf", width = 10, height = 7, paper = "USr")

t <- 1000
plot.fun.realistic.loss(Timestep = t,Nodes=nranl[[1]], Results=resultaki, meantrend = meanaki, upper = upper, lower = lower)

dev.off()



##map and timeseries video, using plot.fun 3

saveGIF({
  for (i in 1:500) plot.fun.realistic.loss(Timestep = i,Nodes=nranl[[1]], Results=resultaki, meantrend = meanaki7,
                                             upper = upper, lower = lower)},
  movie.name = "simulation_video_landscape.gif",interval = 0.1, ani.height = 500, ani.width = 500
)


dev.off()



