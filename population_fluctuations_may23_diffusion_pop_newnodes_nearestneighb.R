library(sp)
library(rgdal)
library(maptools)
library(rgeos)


set.seed(1) #seed for reproducible results
#number of timesteps
tstep <- 160
ph <- c(1:tstep) #vector of timestep values
phases <- c("Phase 1","Phase 2","Phase 3","Phase 4","Phase 5","Phase 6","Phase 7","Phase 8","Phase 9","Phase 10","Phase 11", "Phase 12","Phase 13","Phase 14") #labels for plots per phase

##population parameters
## possibly useless, I can just use the amplitude, popmax <- 200 # maximum population
car <- 1000 # carrying capacity-midline
period <- 100 # if the period is set to be  period = tstep (number of timesteps) the model represents population growth only.
# a, set below is amplitude, if set to 0 (meaning that popmax = c) the model is one of stable population, should be  <= c, otherwise we get negative population values which are meaningless
b <- (2 * pi)/period

#Parameters sweeped
d <- c(0,20, 40) ## horizontal shift possibly important, as it determines initial conditions
horizontal <- d[3]
#kappainc <- 0.01
kappa <- 0.4
ampinc <- 100
amp <- seq(from = 0, to = 200, by = ampinc)


#source node parameters
initial <- 0.1

## load empirical data
directory <- getwd()
aegean <- readOGR(dsn=directory, layer="aeg") #vector map
utm35Nwgs84 <- CRS("+init=epsg:32635") #set crs
aegean <- spTransform(aegean, utm35Nwgs84)

##END DATA LOADING


##BEGIN FUNCTION DEFINITIONS

## Defines the sine function
sin.exp <- function(A,B,C,t,D){A*sin(B*(t+D))+C}

## Defines the cost function
distfun.exp <- function(x,k){exp(-k*x)} #Distance decay function


##diffusion function definition
diffusion <- function(x, wm, cap=1){
  res <- rep(NA,length(x))
  for (i in 1:nrow(wm)){
    res[i] <- x[i] + sum(wm[i,] * x)
  } 
  res[res > 1] <- cap
  return(res)
}


# fluct_nodes: A new function that will take the starting nodes and the *fluctuation* vector and do the following:
# Create a list with the starting nodes as the first element.
# For each element of *fluctuation*, if the value is positive, use spsample() and generate as many points as the value, then rbind # # that to the previous element on the list of nodes and store the whole thing as a new element in that list.
# If the value is negative, take the  previous element in the list, sample randomly as many points as the absolute value of that negative value and then subset the previous element of the list minus the sampled points and store that as a new element in the list.


fluct_nodes <- function(start_nodes,fluctuation){
  slist <- vector(mode = "list", length = tstep) #empty list to intialise
  slist[[1]] <- start_nodes
  for (i in 1:length(fluctuation)){
    if (fluctuation[i] > 0){
      add <- spsample(aegean, n=fluctuation[i], type = "random", iter = 10)
      slist[[i + 1]] <- rbind(slist[[i]],add)
    } else if (fluctuation[i] < 0){
      sample <- sample(x = c(1:length(slist[[i]])),abs(fluctuation[i]))
      dtnodes <- as.matrix(as.data.frame(slist[[i]]))
      subsetnodes <- dtnodes[-sample,]
      subsetnodes <- subsetnodes[,1:2]   
      slist[[i + 1]] <- SpatialPoints(coords = subsetnodes,proj4string = utm35Nwgs84) ##to remove the sampled points
    } else {
      slist[[i + 1]] <- slist[[i]]
    }
  }
  return(slist)
}


##new_node: function to set value for new nodes to be equal to that of their nearest neighbour

#1. takes the new distance matrix for every population increase phase
#2. for each new node, locates the index of the node with the largest value (maximum value on the row)
#3. goes to the trait proportion value vector for the previous phase and takes the value corresponding to that index
#4. assigns the new node this value

new_node <- function(listnewnodes,distmatrix,prevtraitlist){
  for (i in 1:length(listnewnodes)){
    x <- which.max(distmatrix[i+length(prevtraitlist),1:length(prevtraitlist)]) #index of nearest neighbour
    listnewnodes[i] <- prevtraitlist[x]
  }
  return(listnewnodes)
}

##END FUNCTION DEFINITIONS


##BEGIN SIMULATION


##Loop to sweep amplitude 

resultsl <- vector(mode = "list", length = length(amp))
for(e in 1:length(amp)){
  site_number <- round(sin.exp(A=amp[e], B=b, C=car, t=ph, D=horizontal))
  ##plot(site_number, ylim=c(0,500), xlab="Time", ylab="Sites",main = "Population size per phase", sub=paste("amplitude = ",a[e], ", period = ",period,", carrying capacity = ",car), font.sub=3)
  ##lines(site_number)
  fluctuation1 <- round(diff(site_number))
  nodes1 <- spsample(aegean, n=site_number[1], type="random",iter = 10)
  
  negfluctuation <- which(fluctuation1 < 0) #indices of negative fluctuations
  fluctuationneg <- fluctuation1[fluctuation1 < 0] #size of negative fluctuation
  
  # applies fluct_nodes and gets cost matrices
  nodeslist <- fluct_nodes(start_nodes = nodes1,fluctuation = fluctuation1)
  
  for (i in 1:tstep){
    nodeslist[[i]]$siteid <- row.names(nodeslist[[i]]) ##add site id column to each element of nodeslist, to plot   network w. sparch
  }
  
  
  ## Euclidean distance matrix list
  d_matrixl <- vector(mode = "list", length = tstep)
  for (i in 1:tstep){
    d_matrixl[[i]] <- gDistance(nodeslist[[i]], byid=TRUE)/1000
  }
  ## Travel Cost weight matrix list
  w_matrixl <- vector(mode = "list", length = tstep)
  for (i in 1:tstep){
    w_matrixl[[i]] <- apply(d_matrixl[[i]], c(1,2), distfun.exp, k=kappa)
    diag(w_matrixl[[i]]) <- 0 #so that nodes don't increase on their own
  }
  
  rm(d_matrixl) ##remove distance matrix to save space
  #get indices of removed values 
  
  if(all(fluctuation1 >= 0)){
    print("All values are non-negatives!")
  } else {
    removed1 <- list()
    for (i in 1:length(negfluctuation)){
      removed1[[i]] <- sample(x = c(1:length(nodeslist[[negfluctuation[i]]])),abs(fluctuationneg[i]))
    }
    names(removed1) <- negfluctuation      ## this gives a name to each element of removed1 equal to its timestep index
    
    ## This object is the one I need to use for the diffusion process
    removed <- vector(mode = "list", length = length(fluctuation1))
    for (i in 1:length(fluctuation1)){
      iter <- as.character(i)
      if (fluctuation1[i] < 0){
        removed[[i]] <- removed1[[which(names(removed1) == iter)]]
      } else if (fluctuation1[i] >= 0){
        removed[[i]] <- c(0)
      }
    }
  }
  
  
  ## Initialise trait values
  bg <- which.max(nodeslist[[1]]$x) #find index of point with max longitude (easternmost)
  trait <- rep(c(0),each=nrow(nodeslist[[1]]))
  trait[bg] <- initial # this sets the origin point trait proportion value
  trait0 <- trait
  
  
  ## This applies the diffusion function for all timesteps.
  
  trait <- trait0
  nresults2 <- vector(mode = "list", length = tstep) # results list (size is determined by parameter values set at the beginning)
  nresults2[[1]] <- trait0
  for (p in 2:tstep){  
    if (fluctuation1[p-1] > 0){
      add1 <- rep(NA, times=fluctuation1[p-1])
      add2 <- new_node(listnewnodes=add1, distmatrix=w_matrixl[[p]], prevtraitlist=nresults2[[p-1]]) ##applies new_node function
      trait <- c(trait, add2)
    } else if (fluctuation1[p-1] < 0){
      rm <- removed[[p-1]]
      trait <- trait[-rm]  #get the correct index from the removed list
    } else if (fluctuation1[p-1] == 0){
      trait <- trait
    }
    trait <- diffusion(x = trait, wm = w_matrixl[[p]]) 
    nresults2[[p]] <- trait
    print("done 1 tstep") ##test
  }
  
  ## remove NAs and replace them with 0 for nresults2 list
  
  for (j in 1:length(nresults2)){
    nresults2[[j]][is.na(nresults2[[j]])] <- 0
  }
  
  
  ##store trend means in the trendmeans matrix/list TODO figure out best way to store multiple outputs
  resultsl[[e]] <- sapply(X = nresults2, FUN = mean)
  print("done") ##test
}

##END SIMULATION


##BEGIN PLOTTING


## plot the line graphs in the trendmeans matrix as a single plot and colour code based on amplitude
rbPal <- colorRampPalette(c('orange','purple2'))
amplitudecols <- rbPal(length(amp))[as.numeric(cut(amp,breaks = length(amp)))]

## plot and save as pdf

pdf("sinmeanadoptionnearestneighbourexperiment11.pdf")

plot(resultsl[[1]],ylim = c(0,1), col="white", xlab="Time", ylab="Mean trait proportion value",
      main = "k = 0.4, a = 1, carrying capacity = 1000", sub = "Horizontal shift =  40")
  for (i in 1:length(resultsl)){
    lines(resultsl[[i]], type="l", ylim = c(0,1),col=amplitudecols[i],lwd=2)
  }
  legend("right",legend=as.character(amp),fill = amplitudecols, col=amplitudecols, cex=1, title="Amplitude")

dev.off()

sitel <- vector(mode = "list", length = length(amp))
for(e in 1:length(amp)){
  sitel[[e]] <- round(sin.exp(A=amp[e], B=b, C=car, t=ph, D=horizontal))
}

pdf("sinmeanadoptionNovertimeExperiment8nearestneighbour.pdf")

plot(sitel[[1]], ylim=c(0,3000), xlab="Time", ylab="N (number of communities)",main = "N over time, period = 100, c = 1000, d = 40",
     col="white")
for (i in 1:length(amp)){
  lines(sitel[[i]],col=amplitudecols[i],lwd=2)
  legend("topleft",legend=as.character(amp),fill = amplitudecols, col=amplitudecols, cex=1, title="Amplitude")
}

dev.off()

##load libraries for plotting
library(scales)
library(colourvalues)

##FUNCTION DEFINITION 3 BEGINS HERE - NEEDS DEBUGGING 07/03
plot.fun.realistic.loss <- function(Timestep,Nodes,Results, meantrend){
  Nodes$Value <- Results[[Timestep]]
  Nodes$Colour  <- colourvalues::colour_values(Nodes$Value)
  par(mfrow=c(1,2))
  par(mar=c(1, 2, 1, 2))
  plot(aegean, col="white", border="black", lwd=0.1)
  points(Nodes, pch=19,col = alpha(Nodes$Colour, 0.6), cex=1)
  par(mar=c(6, 4, 6, 2))
  plot(meantrend, ylim = c(0,1), cex = 0.1, 
       ylab= "Mean trait proportion value", xlab="Time")
  lines(meantrend, col = "purple", lwd = 1)
  points(x = Timestep , y = meantrend[Timestep], col= "purple2", bg = "yellow", pch=21, ylim = c(0,100), cex=1.6)
  par(mfrow=c(1,1)) 
}  



##Data
meanaki <- resultsl[[3]]


##test, apply function to a single timestep

pdf("t1500_loss_anatolian_origin_a1k0.4maxloss_0_n2000seed1.pdf", width = 10, height = 7, paper = "USr")

t <- 10
plot.fun.realistic.loss(Timestep = t,Nodes=nodeslist[[t]], Results=nresults2, meantrend = meanaki)

dev.off()

##load libraries

library(magick)
library(animation)

##map and timeseries, using plot.fun 3

saveGIF({
  for (i in 1:160) plot.fun.realistic.loss(Timestep = i,Nodes=nodeslist[[i]], Results=nresults2, meantrend = meanaki)},
  movie.name = "popfluct052023k0.2experiment10nn.gif",interval = 0.1, ani.height = 500, ani.width = 500
)


