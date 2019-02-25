presmooth <- function(vec){

  library(MASS)
  library(fda)

  tobs     <- seq(0,23,1)
  knots    <- seq(0,23,1)

  nknots   <- length(knots)
  norder   <- 5
  nbasis   <- length(knots) + norder - 2
  basis    <- create.bspline.basis(c(min(tobs),max(tobs)),nbasis,norder,knots)


  basismat <- eval.basis(tobs, basis)


  # Use quadrature to get integral - Composite Simpson's Rule

  delta <- 0.05
  quadpts <- seq(0,23,delta)
  nquadpts <- length(quadpts)
  quadwts <- as.vector(c(1,rep(c(4,2),(nquadpts-2)/2),4,1),mode="any")
  quadwts <- c(1,rep(c(4,2),(nquadpts-1)/2))
  quadwts[nquadpts] <- 1
  quadwts <- quadwts*delta/3


  # Second derivative of basis functions at quadrature points

  Q2basismat   = eval.basis(quadpts, basis,2);

  # estimates for basis coefficients
  Rmat = t(Q2basismat)%*%(Q2basismat*(quadwts%*%t(rep(1,nbasis))))
  basismat2 = t(basismat)%*%basismat;
  lambda = 1  # smoothing parameter
  Bmat                      = basismat2 + lambda*Rmat;
  chat = ginv(Bmat)%*%t(basismat)%*%vec;

  # fitted value
  yhat = basismat%*%chat;
  return(yhat)
}


sphere.distance <- function(p1,p2){
  return(acos(p1 %*% p2)*180/pi)
}

circmean <- function(circdata){
  x_vec <- cos(circdata * pi / 180)
  y_vec <- sin(circdata * pi / 180)

  x_mean <- sum(x_vec) / length(x_vec)
  y_mean <- sum(y_vec) / length(y_vec)

  cos_mean <- x_mean
  sin_mean <- y_mean

  theta_mean <- atan(sin_mean / cos_mean) * 180 / pi

  if(sin_mean > 0){
      if( cos_mean > 0){
          theta_mean <- theta_mean
      }else{
          theta_mean <- 180 - theta_mean
      }
  }else{
     if( cos_mean > 0){
          theta_mean <- 360 - theta_mean
     }else{
          theta_mean <- 180 + theta_mean
     }
  }

  return(theta_mean)
}



distance.fun <- function(mean.angle,data_vec){
   x_vec <- cos(data_vec * pi / 180)
   y_vec <- sin(data_vec * pi / 180)

   x_mean <- cos(mean.angle * pi / 180)
   y_mean <- sin(mean.angle * pi / 180)
   meanvec <- c(x_mean, y_mean)

   datamat <- cbind(x_vec, y_vec)

   dotprod <- rep(NA, length(data_vec))

   loc <- which(!is.na(data_vec))

   for (jj in 1:length(loc)){
     index = loc[jj]
     dotprod[index] <- datamat[index,] %*% meanvec
   }
   sqr.dist <- mean((acos(dotprod))^2,na.rm=TRUE)
   return(sqr.dist)
}

plot.wind <- function(data_vec,city){
  plot_data <- as.data.frame(data_vec[!is.na(data_vec)])
  colnames(plot_data) <- 'direction'
  # choose bin size (degrees/bin)
  deg <- 30
  # define the range of each bin
  dir.breaks <- seq(0-(deg/2), 360+(deg/2), deg)
  #Now we generate a factor variable, exchanging the directions with the ranges. We’ll also generate some pretty labels and assign them as levels of the new object. Finally we’ll attach the new variable to the main dataset.

  # assign each direction to a bin range
  dir.binned <- cut(plot_data$direction,
                    breaks = dir.breaks,
                    ordered_result = TRUE)
  # generate pretty lables
  dir.labels <- as.character(c(seq(0, 360-deg, by = deg), 0))
  # replace ranges with pretty bin lables
  levels(dir.binned) <- dir.labels
  # Assign bin names to the original data set
  summary(dir.binned)
  plot_data$dir.binned <- dir.binned

  thm <- theme_bw() +
    theme(axis.text.x = element_text(size=8, face = "plain"),
          axis.text.y = element_text(size=8, face = "plain"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=8, face = "plain", hjust = 0.9, vjust = 1.3),
          panel.border = element_blank(),
          panel.grid  = element_blank())


  # initialise the plot
  plt.dirrose <- ggplot() +
    # since the plot background is blank we'll add a series of horizontal lines, at 5000 count intervals, up to 25000.
    geom_hline(yintercept = seq(0, 300, by = 20), colour = "grey60", size = 0.3) +
    # Now we add a darker horizontal line as the top border at 30000.
    geom_hline(yintercept = 300, colour = "black", size = 0.3) +
    # We want 12 vertical lines representing the centers of the 30° ranges.
    geom_vline(xintercept = c(seq(1,12,1)), colour = "grey60", size = 0.3) +
    # On top of everything we place the histogram bars.
    geom_bar(data = plot_data, aes(x = dir.binned), width = 1, colour="black", size = 0.3, alpha=0.5) +
    # Add the x-axis labels
    scale_x_discrete( drop = FALSE, labels = c(0, "", "", 90, "", "", 180, "", "", 270, "", "")) +
    # Add the y-axis labels
    scale_y_continuous(limits = c(0, 300), expand = c(0, 0),
                       breaks = seq(0, 300, by = 20),
                       labels = seq(0, 300, by = 20)) +
    # Add the axis titles
    labs(x = 'Outward step bearing (°)', y = 'Count of outward steps ',
        title=paste('Wind rose plot of ',city,sep='')) +
    # If you only use the plot code up till here you will get a histogram.
    # the next line wraps the histogram into a windrose
    coord_polar(start = -(deg/2)*(pi/180)) +
    # apply theme
    thm

  return(plt.dirrose)

}
get_s_0 <- function(t,h_bw){
  s_0_mat <- matrix(NA, nrow = 651, ncol = 24)
  for (ii in 1:651){
    for (jj in 1:24){
      s_0_mat[ii,jj] <- (1/h_bw)*epane((timemat[ii,jj] - t)/h_bw)
    }
  }
  return(sum(s_0_mat,na.rm=TRUE)/nrow(s_0_mat)/ncol(s_0_mat))
}
get_s_1 <- function(t,h_bw){
  s_1_mat <- matrix(NA, nrow = 651, ncol = 24)
  for (ii in 1:651){
    for (jj in 1:24){
      s_1_mat[ii,jj] <- (1/h_bw)*epane((timemat[ii,jj] - t)/h_bw)*((timemat[ii,jj] - t)/h_bw)
    }
  }
  return(sum(s_1_mat,na.rm=TRUE)/nrow(s_1_mat)/ncol(s_1_mat))
}

get_s_2 <- function(t,h_bw){
  s_2_mat <- matrix(NA, nrow = 651, ncol = 24)
  for (ii in 1:651){
    for (jj in 1:24){
      s_2_mat[ii,jj] <- (1/h_bw)*epane((timemat[ii,jj] - t)/h_bw)*(((timemat[ii,jj] - t)/h_bw)^2)
    }
  }
  return(sum(s_2_mat,na.rm=TRUE)/nrow(s_2_mat)/ncol(s_2_mat))
}

get_r_0 <- function(t,obs,h_bw){
  r_0_mat <- matrix(NA, nrow = 651, ncol = 24)
  for (ii in 1:651){
    for (jj in 1:24){
      r_0_mat[ii,jj] <- (1/h_bw)*epane((timemat[ii,jj] - t)/h_bw)*obs[ii,jj]
    }
  }
  return(sum(r_0_mat,na.rm=TRUE)/nrow(r_0_mat)/ncol(r_0_mat))
}

get_r_1 <- function(t,obs,h_bw){
  r_1_mat <- matrix(NA, nrow = 651, ncol = 24)
  for (ii in 1:651){
    for (jj in 1:24){
      r_1_mat[ii,jj] <- (1/h_bw)*epane((timemat[ii,jj] - t)/h_bw)*((timemat[ii,jj] - t)/h_bw)*obs[ii,jj]
    }
  }
  return(sum(r_1_mat,na.rm=TRUE)/nrow(r_1_mat)/ncol(r_1_mat))
}

epane  <- function(x){
  if (abs(x) <= 1){
    return(3/4*(1-x^2))
  }else{
    return(0)
  }

}


llr.mean <- function(time,data,bw){
  s0_vec <- rep(NA,length(time))
  for (ii in 1:length(time)){
    s0_vec[ii] <- get_s_0(ii,h_bw =bw)
  }

  s1_vec <- rep(NA,length(time))
  for (ii in 1:length(time)){
    s1_vec[ii] <- get_s_1(ii,h_bw =bw)
  }

  s2_vec <- rep(NA,length(time))
  for (ii in 1:length(time)){
    s2_vec[ii] <- get_s_2(ii,h_bw =bw)
  }

  r0_vec <- rep(NA,length(time))
  for (ii in 1:length(time)){
    r0_vec[ii] <- get_r_0(ii,obs =data,h_bw =bw)
  }

  r1_vec <- rep(NA,length(time))
  for (ii in 1:length(time)){
    r1_vec[ii] <- get_r_1(ii,obs = data,h_bw =bw)
  }

  return((r0_vec * s2_vec - r1_vec * s1_vec)/(s0_vec*s2_vec - (s1_vec)^2))
}

S_pq <- function(p,q,s,t,h_R,timemat){
  tempmat <- matrix(NA, nrow = nrow(timemat), ncol = ncol(timemat))
  sum_vec <- rep(NA, nrow = nrow(timemat))
  for (ii in 1:nrow(timemat)){
    index = 1
    tempvec <- rep(NA, length(timemat[ii,])*( length(timemat[ii,])-1))
    for (kk in 1:ncol(timemat)){
      for (jj in 1:ncol(timemat)){
        if( jj != kk){
          tempvec[index] <- (((timemat[ii,jj] - s)/h_R)^p)*(((timemat[ii,kk] - t)/h_R)^q)*
            (1/h_R)*epane((timemat[ii,jj] - s)/h_R)*
            (1/h_R)*epane((timemat[ii,kk] - t)/h_R)
          index = index + 1
        }else{

        }
      }
    }
    sum_vec[ii] <- sum(tempvec)
  }
  return(sum(sum_vec) / length(sum_vec) / nrow(timemat))
}


R_pq <- function(p,q,s,t,h_R,timemat,obs){
  tempmat <- matrix(NA, nrow = nrow(timemat), ncol = ncol(timemat))
  sum_vec <- rep(NA, nrow = nrow(timemat))
  for (ii in 1:nrow(timemat)){
    index = 1
    tempvec <- rep(NA, length(timemat[ii,])*( length(timemat[ii,])-1))
    for (kk in 1:ncol(timemat)){
      for (jj in 1:ncol(timemat)){
        if( jj != kk){
          tempvec[index] <- obs[ii,jj] * obs[ii,kk]*(((timemat[ii,jj] - s)/h_R)^p)*(((timemat[ii,kk] - t)/h_R)^q)*
            (1/h_R)*epane((timemat[ii,jj] - s)/h_R)*
            (1/h_R)*epane((timemat[ii,kk] - t)/h_R)
          index = index + 1
        }else{

        }
      }
    }
    sum_vec[ii] <- sum(tempvec,na.rm= TRUE)
  }
  return(sum(sum_vec) / length(sum_vec) / nrow(timemat))
}


getcov <- function(s,t,mean_vec,h_R,obs, timemat){
    S_00  <- S_pq(p = 0 , q = 0 , s = s , t = t, h_R = h_R, timemat = timemat)
    S_01  <- S_pq(p = 0 , q = 1 , s = s , t = t, h_R = h_R, timemat = timemat)
    S_10  <- S_pq(p = 1 , q = 0 , s = s , t = t, h_R = h_R, timemat = timemat)
    S_11  <- S_pq(p = 1 , q = 1 , s = s , t = t, h_R = h_R, timemat = timemat)
    S_02  <- S_pq(p = 0 , q = 2 , s = s , t = t, h_R = h_R, timemat = timemat)
    S_20  <- S_pq(p = 2 , q = 0 , s = s , t = t, h_R = h_R, timemat = timemat)
    R_00  <- R_pq(p = 0 , q = 0 , s = s , t = t, h_R = h_R, timemat = timemat, obs = obs)
    R_01  <- R_pq(p = 0 , q = 1 , s = s , t = t, h_R = h_R, timemat = timemat, obs = obs)
    R_10  <- R_pq(p = 1 , q = 0 , s = s , t = t, h_R = h_R, timemat = timemat, obs = obs)

    A_1   <- S_20 * S_02 - (S_11)^2
    A_2   <- S_10 * S_02 - S_01 * S_11
    A_3   <- S_01 * S_20 - S_10 * S_11
    B     <- A_1 * S_00 - A_2 * S_10 - A_3 * S_01

    cov.val <- (1/B) * (A_1 * R_00 - A_2 * R_10 - A_3 * R_01) - mean_vec[s+1] * mean_vec[t+1]
    return(cov.val)
}

getcrosscov <- function(s, t, h_R,obs_1,obs_2, timemat){
  S_00  <- S_pq(p = 0 , q = 0 , s = s , t = t, h_R = h_R, timemat = timemat)
  S_01  <- S_pq(p = 0 , q = 1 , s = s , t = t, h_R = h_R, timemat = timemat)
  S_10  <- S_pq(p = 1 , q = 0 , s = s , t = t, h_R = h_R, timemat = timemat)
  S_11  <- S_pq(p = 1 , q = 1 , s = s , t = t, h_R = h_R, timemat = timemat)
  S_02  <- S_pq(p = 0 , q = 2 , s = s , t = t, h_R = h_R, timemat = timemat)
  S_20  <- S_pq(p = 2 , q = 0 , s = s , t = t, h_R = h_R, timemat = timemat)
  Rc_00  <- Rc_pq(p = 0 , q = 0 , s = s , t = t, h_R = h_R, timemat = timemat, obs_1 = obs_1, obs_2 = obs_2)
  Rc_01  <- Rc_pq(p = 0 , q = 1 , s = s , t = t, h_R = h_R, timemat = timemat, obs_1 = obs_1, obs_2 = obs_2)
  Rc_10  <- Rc_pq(p = 1 , q = 0 , s = s , t = t, h_R = h_R, timemat = timemat,  obs_1 = obs_1, obs_2 = obs_2)

  A_1   <- S_20 * S_02 - (S_11)^2
  A_2   <- S_10 * S_02 - S_01 * S_11
  A_3   <- S_01 * S_20 - S_10 * S_11
  B     <- A_1 * S_00 - A_2 * S_10 - A_3 * S_01

  cov.val <- (1/B) * (A_1 * Rc_00 - A_2 * Rc_10 - A_3 * Rc_01)
  return(cov.val)
}


Rc_pq <- function(p,q,s,t,h_R,timemat,obs_1,obs_2){
  tempmat <- matrix(NA, nrow = nrow(timemat), ncol = ncol(timemat))
  sum_vec <- rep(NA, nrow = nrow(timemat))
  for (ii in 1:nrow(timemat)){
    index = 1
    tempvec <- rep(NA, length(timemat[ii,])*( length(timemat[ii,])-1))
    for (kk in 1:ncol(timemat)){
      for (jj in 1:ncol(timemat)){
        if( jj != kk){
          tempvec[index] <- obs_1[ii,jj] * obs_2[ii,kk]*(((timemat[ii,jj] - s)/h_R)^p)*(((timemat[ii,kk] - t)/h_R)^q)*
            (1/h_R)*epane((timemat[ii,jj] - s)/h_R)*
            (1/h_R)*epane((timemat[ii,kk] - t)/h_R)
          index = index + 1
        }else{

        }
      }
    }
    sum_vec[ii] <- sum(tempvec,na.rm= TRUE)
  }
  return(sum(sum_vec) / length(sum_vec) / nrow(timemat))
}



getdata <- function(x){

  #----------------------------------------------------
  data.templist <- list()
  for (index in 1:252){
    filename <- paste('~/Dropbox/Wind Data/File ',index,'.csv',sep='')
    #filename <- paste('C:/Users/Inori/Dropbox/Wind Data/File ',index,'.csv',sep='')
    tempdata <- read.csv(filename,skip=15,header=TRUE)
    tempdata <- tempdata[,c('Year','Month','Day','Time','Wind.Dir..10s.deg.','Temp...C.','Wind.Spd..km.h.','Stn.Press..kPa.')]
    colnames(tempdata) <- c('Year','Month','Day','Time','Wind Dir','Temperature','Wind Spd','AtmPress')
    data.templist[[index]] <- tempdata
  }

  data <- do.call(rbind,data.templist)

  dec.data <- data[data$Month==12,]
  numofrep <- nrow(dec.data)/24
  spd.data <- matrix(NA, nrow = numofrep, ncol = 24)
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
      spd.data[index,]  <- tempmat$'Wind Spd'
      datamat[index,]   <- tempmat$`Wind Dir`
      indicator[index,1] <- year[jj]
      indicator[index,2] <- day[kk]
      index <- index + 1
    }
  }

  van.data <- list(spd.data = spd.data, ang.data = datamat)
  #----------------------------------------------------
  #----------------------------------------------------
  data.templist <- list()
  for (index in 1:21){
    filename <- paste('~/Dropbox/toronto wind data/Toronto File',index,'.csv',sep='')
    #filename <- paste('C:/Users/Inori/Dropbox/Wind Data/File ',index,'.csv',sep='')
    tempdata <- read.csv(filename,skip=15,header=TRUE)
    tempdata <- tempdata[,c('Year','Month','Day','Time','Wind.Dir..10s.deg.','Temp...C.','Wind.Spd..km.h.','Stn.Press..kPa.')]
    colnames(tempdata) <- c('Year','Month','Day','Time','Wind Dir','Temperature','Wind Spd','AtmPress')
    data.templist[[index]] <- tempdata
  }

  data <- do.call(rbind,data.templist)

  dec.data <- data[data$Month==12,]
  numofrep <- nrow(dec.data)/24
  spd.data <- matrix(NA, nrow = numofrep, ncol = 24)
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
      spd.data[index,]  <- tempmat$'Wind Spd'
      datamat[index,]   <- tempmat$`Wind Dir`
      indicator[index,1] <- year[jj]
      indicator[index,2] <- day[kk]
      index <- index + 1
    }
  }
  ton.data <- list(spd.data = spd.data, ang.data = datamat)

  #----------------------------------------------------
  data.templist <- list()
  for (index in 1:21){
    filename <- paste('~/Dropbox/halifax wind data/Halifax File',index,'.csv',sep='')
    #filename <- paste('C:/Users/Inori/Dropbox/Wind Data/File ',index,'.csv',sep='')
    tempdata <- read.csv(filename,skip=15,header=TRUE)
    tempdata <- tempdata[,c('Year','Month','Day','Time','Wind.Dir..10s.deg.','Temp...C.','Wind.Spd..km.h.','Stn.Press..kPa.')]
    colnames(tempdata) <- c('Year','Month','Day','Time','Wind Dir','Temperature','Wind Spd','AtmPress')
    data.templist[[index]] <- tempdata
  }

  data <- do.call(rbind,data.templist)

  dec.data <- data[data$Month==12,]
  numofrep <- nrow(dec.data)/24
  spd.data <- matrix(NA, nrow = numofrep, ncol = 24)
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
      spd.data[index,]  <- tempmat$'Wind Spd'
      datamat[index,]   <- tempmat$`Wind Dir`
      indicator[index,1] <- year[jj]
      indicator[index,2] <- day[kk]
      index <- index + 1
    }
  }
  hal.data <- list(spd.data = spd.data, ang.data = datamat)
  #-------------------------------------------------------------

  return(list(van.data = van.data,
              ton.data = ton.data,
              hal.data = hal.data))
}


logmap <- function(datamat,nsteps = 5000,initial.var = 120){
  logmap_data <- matrix(NA, nrow = nrow(datamat),ncol = ncol(datamat))
  mean_vec    <- rep(NA, 24)
  rotate.angle <- rep(NA,24)

  for (ii in 1:24){
    print(paste('Hour',ii))
    data_vec          <- datamat[,ii] * 10

     current.theta <- initial.var
     nsteps        <- nsteps
     learn.rate    <- 2
     stepsize      <- .1
     jj = 1
     while (jj <= nsteps){
       print(jj)
       grad.dist     <-  (distance.fun(current.theta + stepsize, data_vec = data_vec) - distance.fun(current.theta - stepsize, data_vec = data_vec))/(2*stepsize)
       current.theta <- current.theta - learn.rate * grad.dist
       jj            <- jj + 1
     }


    mean_vec[ii]   <- current.theta
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


    tempvec <- rep(NA, length(data_vec))
    for (jj in 1:length(nonna)){
      tempvec[nonna[jj]] <- rotated.points[nonna[jj],1]
    }
    logmap_data[,ii] <- tempvec
  }
  return(list(logmap_data = logmap_data, mean_vec = mean_vec, rotate.angle = rotate.angle))
}

exp_map <- function(datamat,logmap_data, mean_vec, rotate.angle){
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

    line.vec    <- logmap_data[,ii]
    line.points <- cbind(line.vec, rep(0,length(line.vec)))

    plot(line.points,main=ii)
    angle <- rotate.angle[ii]
    anti.rotate.mat <- matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)), byrow = TRUE,nrow = 2)
    anti.rotate <-matrix(NA, nrow = nrow(tempmat), ncol = ncol(tempmat))
    for (index in 1:length(logmap_data[,ii])){
      prep_vec    <- line.points[index,]
      anti.rotate[index,] <- anti.rotate.mat %*% prep_vec + tempmean
    }
    points(anti.rotate[,1],anti.rotate[,2],col='blue')
    points(tempmat[,1],tempmat[,2],col='red')



    shifted.points <- matrix(NA, nrow = nrow(anti.rotate), ncol = ncol(anti.rotate))
    expmap.data[[ii]]    <- matrix(NA, nrow = 651, ncol = 2)
    for (jj in 1:nrow(datamat)){
      shifted.points[jj,] <- anti.rotate[jj,] - tempmean
      tang.vec <- shifted.points[jj,]
      expmap.data[[ii]][jj,] <- cos(sqrt(tang.vec %*% tang.vec)) * mean_coord  +
        sin(sqrt(tang.vec %*% tang.vec))/(sqrt(tang.vec %*% tang.vec)) * tang.vec
      x_coord <- expmap.data[[ii]][jj,1]
      y_coord <- expmap.data[[ii]][jj,2]

      if(x_coord > 0 ){
        if(y_coord >0){
          expmap.ang[jj,ii]  <- atan(expmap.data[[ii]][jj,2]/expmap.data[[ii]][jj,1])
        }else{
          expmap.ang[jj,ii]  <- atan(expmap.data[[ii]][jj,2]/expmap.data[[ii]][jj,1]) + 2*pi
        }
      }else{
        if(y_coord > 0){
          expmap.ang[jj,ii]  <- atan(expmap.data[[ii]][jj,2]/expmap.data[[ii]][jj,1]) + pi
        }else{
          expmap.ang[jj,ii]  <- atan(expmap.data[[ii]][jj,2]/expmap.data[[ii]][jj,1]) + pi
        }
      }

    }

  }

  return(expmap = expmap.ang)

}




