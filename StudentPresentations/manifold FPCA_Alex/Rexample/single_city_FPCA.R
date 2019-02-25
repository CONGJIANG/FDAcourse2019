rm(list=ls())
source('functions.R')
#1. Get data for December
data.templist <- list()
for (index in 1:252){
  filename <- paste('~/Dropbox/Wind Data/File ',index,'.csv',sep='') 
  #filename <- paste('C:/Users/Inori/Dropbox/Wind Data/File ',index,'.csv',sep='') 
  tempdata <- read.csv(filename,skip=15,header=TRUE)
  tempdata <- tempdata[,c('Year','Month','Day','Time','Wind.Dir..10s.deg.')]
  colnames(tempdata) <- c('Year','Month','Day','Time','Wind Dir')
  data.templist[[index]] <- tempdata
}

data <- do.call(rbind,data.templist)



#
dec.data <- data[data$Month==12,]
numofrep <- nrow(dec.data)/24

datamat   <- matrix(NA, nrow = numofrep, ncol = 24)
rownames  <- rep(NA, numofrep)
year      <- seq(1998,2018,1)
day       <- seq(1,31,1)
index     <- 1
indicator <- matrix(NA, nrow = numofrep, ncol = 2)
for (jj in 1:21){
  for (kk in 1:31){
    tempmat           <- dec.data[dec.data$Year == year[jj],]
    tempmat           <- tempmat[tempmat$Day  == day[kk],]
    datamat[index,]   <- tempmat$`Wind Dir`
    indicator[index,1] <- year[jj]
    indicator[index,2] <- day[kk]
    index <- index + 1
  }
}
#--------------------------------------------------------------------------
#2. Wind rose plot
library(ggplot2)
plot.wind(datamat[,1] * 10,city='Vancouver')


#--------------------------------------------------------------------------
#3. log-map data around frechet mean 
logmap_data <- matrix(NA, nrow = nrow(datamat),ncol = ncol(datamat))
rotate.angle <- rep(NA,24)

for (ii in 1:24){
  print(paste('Hour',ii))
  data_vec          <- datamat[,ii] * 10

  # current.theta <- 120
  # nsteps        <- 5000
  # learn.rate    <- 2
  # stepsize      <- .1
  # jj = 1
  # while (jj <= nsteps){
  #   print(jj)
  #   grad.dist     <-  (distance.fun(current.theta + stepsize, data_vec = data_vec) - distance.fun(current.theta - stepsize, data_vec = data_vec))/(2*stepsize)
  #   current.theta <- current.theta - learn.rate * grad.dist
  #   jj            <- jj + 1
  # }

 
  #mean_vec[ii]   <- current.theta
  mean_coord     <- c(cos(mean_vec[ii] * pi / 180), sin(mean_vec[ii] * pi / 180))
  # plot(0,0,xlim=c(-1,1),ylim=c(-1,1),pch=' ')
  # draw.circle(x = 0,y = 0,radius = 1)
  # points(x = mean_coord[1],y = mean_coord[2],pch=16,col='red')
  # 
  # 
  nonna   <- which(!is.na(data_vec))
  tempmat <- matrix(NA,nrow = length(data_vec),ncol =2 )
  for (index in 1:length(nonna)){
    x_vec <- c(cos(data_vec[nonna[index]] * pi / 180),sin(data_vec[nonna[index]] * pi / 180))
    vec_u <- x_vec - (mean_coord %*% x_vec) * mean_coord
    v_vec <- (vec_u / sqrt(vec_u %*% vec_u)) * acos(mean_coord %*% x_vec)
    tempmat[nonna[index],] <- v_vec
  }
  # points(x = tempmat[2,1], y = tempmat[2,2],col='blue',pch = '1')
  # points(x = tempmat[3,1], y = tempmat[3,2],col='blue',pch = '2')
  # points(x = cos(data_vec[2] * pi / 180), y = sin(data_vec[2] * pi / 180),pch='1')
  # points(x = cos(data_vec[3] * pi / 180), y = sin(data_vec[3] * pi / 180),pch='2')
  # 
  
  tempmean <- apply(tempmat,2,function(x) return(mean(x,na.rm=TRUE)))

  plot(x = tempmat[,1],y=tempmat[,2],xlim=c(-3,3),ylim=c(-3,3),main=ii,ylab='')  
  points(x = tempmean[1],y=tempmean[2],xlim=c(-1,1),ylim=c(-1,1),pch = 16)
  
  # 
  if (coef(lm(tempmat[,2] ~ tempmat[,1]))[2] > 0){
    angle <- 2*pi - (atan(coef(lm(tempmat[,2] ~ tempmat[,1]))[2]))
  }else{
    angle <- - (atan(coef(lm(tempmat[,2] ~ tempmat[,1]))[2]))
  }
  
  rotate.angle[ii] <- angle
  rotation.mat <- matrix(c(cos(angle),-sin(angle),sin(angle),cos(angle)), byrow = TRUE,nrow = 2)
  
  rotated.points <- matrix(NA, nrow = nrow(tempmat), ncol = ncol(tempmat))
  for (index in 1:length(nonna)){
    prep_vec <- tempmat[nonna[index],]
    prep_vec <- prep_vec - tempmean
    rotated.points[nonna[index],] <-  rotation.mat   %*% prep_vec
  }
  lines(rotated.points[,1],rotated.points[,2])
  mean(rotated.points[,1],na.rm=TRUE)
  # 
  # anti.rotate.mat <- matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)), byrow = TRUE,nrow = 2)
  # anti.rotate <-matrix(NA, nrow = nrow(tempmat), ncol = ncol(tempmat))
  # for (index in 1:length(nonna)){
  #   prep_vec    <- rotated.points[nonna[index],]
  #   anti.rotate[nonna[index],] <- anti.rotate.mat %*% prep_vec
  # }
  # lines(anti.rotate[,1],anti.rotate[,2],col='blue')
  
  tempvec <- rep(NA, length(data_vec))
  for (jj in 1:length(nonna)){
    tempvec[nonna[jj]] <- rotated.points[nonna[jj],1]
  }
  logmap_data[,ii] <- tempvec
}
#--------------------------------------------------------------------------
#3. FPCA on logmapped data
par(mfrow=c(2,1))
plot(seq(0,23,1),logmap_data[1,],type='l')
plot(seq(0,23,1),datamat[1,],type='l')
par(mfrow=c(1,1))

#fdapace
library(fdapace)
L3 <- MakeFPCAInputs(IDs = rep(1:651, each=24), tVec=rep(seq(0,23,1),651), t(logmap_data))
FPCAdense <- FPCA(L3$Ly, L3$Lt, opt=list(FVEthreshold = 0.99))
#FPCAdense <- FPCA(L3$Ly, L3$Lt, opt=list(maxK =3 ))


# Plot the FPCA object
plot(FPCAdense)
sqrt(FPCAdense$lambda)


est.mat  <- FPCAdense$xiEst %*% t(FPCAdense$phi)

#4.Exponential map submanifold
library(plotrix)
expmap.data <- list()
expmap.ang  <- matrix(NA, nrow = 651, ncol = 24)
for (ii in 1:24){
  data_vec       <- datamat[,ii] * 10
  mean_coord     <- c(cos(mean_vec[ii] * pi / 180), sin(mean_vec[ii] * pi / 180))
 
  nonna   <- which(!is.na(data_vec))
  tempmat <- matrix(NA,nrow = length(data_vec),ncol =2 )
  for (index in 1:length(nonna)){
    x_vec <- c(cos(data_vec[nonna[index]] * pi / 180),sin(data_vec[nonna[index]] * pi / 180))
    vec_u <- x_vec - (mean_coord %*% x_vec) * mean_coord
    v_vec <- (vec_u / sqrt(vec_u %*% vec_u)) * acos(mean_coord %*% x_vec)
    tempmat[nonna[index],] <- v_vec
  }

   tempmean <- apply(tempmat,2,function(x) return(mean(x,na.rm=TRUE)))
  
   line.vec    <- est.mat[,ii]
   line.points <- cbind(line.vec, rep(0,length(line.vec)))
   
   plot(line.points,main=ii)
   angle <- rotate.angle[ii]
   anti.rotate.mat <- matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)), byrow = TRUE,nrow = 2)
   anti.rotate <-matrix(NA, nrow = nrow(tempmat), ncol = ncol(tempmat))
   for (index in 1:length(nonna)){
     prep_vec    <- line.points[nonna[index],]
     anti.rotate[nonna[index],] <- anti.rotate.mat %*% prep_vec + tempmean
   }
   points(anti.rotate[,1],anti.rotate[,2],col='blue')
   points(tempmat[,1],tempmat[,2],col='red')
  
  
  
  shifted.points <- matrix(NA, nrow = nrow(anti.rotate), ncol = ncol(anti.rotate))
  expmap.data[[ii]]    <- matrix(NA, nrow = 651, ncol = 2)
  for (jj in 1:length(nonna)){
    shifted.points[nonna[jj],] <-anti.rotate[nonna[jj],] - tempmean
    tang.vec <- shifted.points[nonna[jj],] 
    expmap.data[[ii]][nonna[jj],] <- cos(sqrt(tang.vec %*% tang.vec)) * mean_coord  + 
    sin(sqrt(tang.vec %*% tang.vec))/(sqrt(tang.vec %*% tang.vec)) * tang.vec
    x_coord <- expmap.data[[ii]][nonna[jj],1]
    y_coord <- expmap.data[[ii]][nonna[jj],2]
    
    if(x_coord > 0 ){
      if(y_coord >0){
        expmap.ang[nonna[jj],ii]  <- atan(expmap.data[[ii]][nonna[jj],2]/expmap.data[[ii]][nonna[jj],1]) 
      }else{
        expmap.ang[nonna[jj],ii]  <- atan(expmap.data[[ii]][nonna[jj],2]/expmap.data[[ii]][nonna[jj],1]) + 2*pi
      }
    }else{
      if(y_coord > 0){
        expmap.ang[nonna[jj],ii]  <- atan(expmap.data[[ii]][nonna[jj],2]/expmap.data[[ii]][nonna[jj],1]) + pi
      }else{
        expmap.ang[nonna[jj],ii]  <- atan(expmap.data[[ii]][nonna[jj],2]/expmap.data[[ii]][nonna[jj],1]) + pi
      }
    }
      
  }

}

plot.index = 5
plot(seq(0,23,1), datamat[plot.index,] * 10,col='blue')
points(seq(0,23,1),expmap.ang[plot.index,]*180 / pi,type='l',lwd=2)

