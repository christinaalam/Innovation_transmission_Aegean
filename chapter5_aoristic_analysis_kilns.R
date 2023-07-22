##chapter 5, kiln aoristic analysis


## Aoristic Analysis


##load archSeries library

library(archSeries)


##load kiln data

d <- read.csv("KilnPresAoristic_Dec2022.csv")
kilns <- d[,c(1,3,4,5,6,10)]
kilns <- kilns[kilns$remove == 0,]
kilnss <- kilns[,1:3]
kilnss$From <- kilnss$From*(-1)
kilnss$To <- kilnss$To*(-1)


##rename From and To columns to Start and End because otherwise the function won't work

kilnss$Start <- kilnss$From
kilnss$End <- kilnss$To
kilnssss <- data.table(kilnss[,4:5])
aorist_kilns <- archSeries::aorist(kilnssss, bin.width = 1)

aoristic_weights <- aorist_kilns[[3]]
aoristic_dates <- c((min(kilnssss$Start)+1):max(kilnssss$End))
length(aoristic_dates)==length(aoristic_weights)

##plot 

pdf("aoristic_kilns.pdf")

plot(x = aoristic_dates, y = aoristic_weights, type = "n",xaxt = 'n', xlab = "Time (Years BC)", ylab = "aoristic weight sum")
polygon(c(aoristic_dates, rev(aoristic_dates)), c(aoristic_weights,rev(aoristic_weights)), density = 200, col = "lightblue"
        , border = "lightblue", lwd = 3)
axis(side=1,at=aoristic_dates,labels=paste((aoristic_dates)*-1), cex.axis=0.6)


dev.off()

##aoristic weight for each location over time, formula after
## Crema (2012:447)

##function calculating aoristic weight per year for each kiln location

aoristic.weight.fun <- function(timestep=1, start, end){
  weight <- timestep/(end-start)
  return(weight)
}


##apply to kiln data, using a 10 year interval

site_aoristic <- aoristic.weight.fun(timestep = 1, start = kilnss$Start, end = kilnss$End)
kilnss$aoristic_w <- site_aoristic

##now add each time interval value as a column to each kiln
kilnss_mat <- as.matrix(kilnss[,c(-1, -4, -5)])

tm <- vector(mode = "list")
for (i in 1:nrow(kilnss)){
  bef <- rep(0, times = (min(kilnssss$Start) - kilnss_mat[i,1])*-1)
  pres <- rep(site_aoristic[i], times = kilnss_mat[i,2] - kilnss_mat[i,1]) 
  aft <- rep(0, times = (max(kilnssss$End) - kilnss_mat[i,2]))
  tm[[i]] <- c(bef, pres, aft)
}

## Temporal maps with aoristic weights


##load libraries

library(rgdal)
library(maptools)
library(rgeos)
library(colourvalues)
library(magick)
library(animation)
library(scales)
library(sf)


points2 <- st_as_sf(kilns[,4:5], coords = c("Long", "Lat"), crs = 4326)
points3 <- st_transform(points2, 2100) ##project to greek grid
points4 <- as(points3, "Spatial")

tm_points <- do.call("rbind", tm) ## each row is a kiln, each column is the aoristic weight at a specific year
plot(points4)


##load aegean data

directory <- getwd()
aegean <- readOGR(dsn=directory, layer="aeg") #vector map
utm35Nwgs84 <- CRS("+init=epsg:32635") #set crs
aegean <- spTransform(aegean, utm35Nwgs84) #project to specified crs
points5 <- spTransform(points4, utm35Nwgs84) #project to specified crs


##Plotting function
plot.fun <- function(Timestep,Nodes,weights,region, date){ 
  Nodes$Value <- weights[,Timestep]
  Nodes$Colour  <- colourvalues::colour_values(Nodes$Value, palette = "reds")
  par(mar=c(1, 2, 1, 2))
  plot(region, col="white", border="grey", lwd=0.1,
       main = paste("Date : ", -date, "BC"))
  points(Nodes, pch=19,col = alpha(Nodes$Colour, 0.4), cex=3)
}  

##plot

dates <- cbind(aoristic_dates, c(1:length(aoristic_dates)))

##add date you are interested in, using negative numbers for BC dates

mapdate <- -2400
timestep <- which(dates[,1]==mapdate)

pdf("aoristic_kiln_map_2400BC.pdf")
plot.fun(Timestep=timestep, Nodes = points5, weights = tm_points, region = aegean, date = mapdate)
dev.off()
