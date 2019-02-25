rm(list=ls())
source('functions.R')

combined.data <- getdata()


############
#wind rose plots 
library(ggplot2)
library(cowplot)
windplot.van <- plot.wind(combined.data$van.data$ang.data[,1] * 10,city = 'Vancouver')
windplot.ton <- plot.wind(combined.data$ton.data$ang.data[,1] * 10,city = 'Tonroto'  )
windplot.hal <- plot.wind(combined.data$hal.data$ang.data[,1] * 10,city = 'Halifax'  )
quartz()
plot_grid(windplot.van,
          windplot.ton,
          windplot.hal,ncol = 3)

############
#log transform data 
van.logmapped.data <- logmap(combined.data$van.data$ang.data,initial.var = 100)
ton.logmapped.data <- logmap(combined.data$ton.data$ang.data,initial.var = 265)
hal.logmapped.data <- logmap(combined.data$hal.data$ang.data,initial.var = 300)
#-------------
#3. FPCA on logmapped data
par(mfrow=c(2,1))
plot(seq(0,23,1),ton.logmapped.data$logmap_data[1,],type='l')
plot(seq(0,23,1),combined.data$ton.data$ang.data[1,],type='l')
par(mfrow=c(1,1))

#fdapace
load('all_logdata.Rdata')
library(fdapace)
L1      <- MakeFPCAInputs(IDs = rep(1:651, each=24), tVec=rep(seq(0,23,1),651), t(van.logmapped.data$logmap_data)) 
L2      <- MakeFPCAInputs(IDs = rep(1:651, each=24), tVec=rep(seq(0,23,1),651), t(ton.logmapped.data$logmap_data)) 
L3      <- MakeFPCAInputs(IDs = rep(1:651, each=24), tVec=rep(seq(0,23,1),651), t(hal.logmapped.data$logmap_data)) 

(van.pca <- FPCA(L1$Ly, L1$Lt, opt=list(FVEthreshold = 0.99,kernel='epan',methodMuCovEst = 'smooth')))
(ton.pca <- FPCA(L2$Ly, L2$Lt, opt=list(FVEthreshold = 0.99,kernel='epan',methodMuCovEst = 'smooth')))
(hal.pca <- FPCA(L3$Ly, L3$Lt, opt=list(FVEthreshold = 0.99,kernel='epan',methodMuCovEst = 'smooth')))


plot(van.pca)
plot(ton.pca)
plot(hal.pca)

#estimated covariance surface
library(plotly)
library(shiny)
plot_ly(z = ~van.pca$fittedCov) %>% add_surface( contours = list(
  z = list(
    show=TRUE,
    usecolormap=TRUE,
    highlightcolor="#ff0000",
    project=list(z=TRUE)
  )
)
) %>%
  layout(
    scene = list(
      camera=list(
        eye = list(x=1.87, y=0.88, z=-0.64)
      )
    )
  )
plot_ly(z = ~ton.pca$fittedCov) %>% add_surface( contours = list(
  z = list(
    show=TRUE,
    usecolormap=TRUE,
    highlightcolor="#ff0000",
    project=list(z=TRUE)
  )
)
) %>%
  layout(
    scene = list(
      camera=list(
        eye = list(x=1.87, y=0.88, z=-0.64)
      )
    )
  )
plot_ly(z = ~hal.pca$fittedCov) %>% add_surface( contours = list(
  z = list(
    show=TRUE,
    usecolormap=TRUE,
    highlightcolor="#ff0000",
    project=list(z=TRUE)
  )
)
) %>%
  layout(
    scene = list(
      camera=list(
        eye = list(x=1.87, y=0.88, z=-0.64)
      )
    )
  )

#plot some eigenfunctions 

quartz()
par(mfrow=c(3,5))
plot(-van.pca$phi[,1],type='l')
plot(van.pca$phi[,2],type='l')
plot(van.pca$phi[,3],type='l')
plot(van.pca$phi[,4],type='l')
plot(van.pca$phi[,5],type='l')

plot(ton.pca$phi[,1],type='l')
plot(-ton.pca$phi[,2],type='l')
plot(ton.pca$phi[,3],type='l')
plot(ton.pca$phi[,4],type='l')
plot(-ton.pca$phi[,5],type='l')

plot(hal.pca$phi[,1],type='l')
plot(hal.pca$phi[,2],type='l')
plot(hal.pca$phi[,3],type='l')
plot(hal.pca$phi[,4],type='l')
plot(hal.pca$phi[,5],type='l')


#score * eigen functions 
van.est.mat  <- van.pca$xiEst %*% t(van.pca$phi)
ton.est.mat  <- ton.pca$xiEst %*% t(ton.pca$phi)
hal.est.mat  <- hal.pca$xiEst %*% t(hal.pca$phi)


#plot estimated logmap data 
par(mfrow=c(1,3))
plot(van.logmapped.data$logmap_data[1,])
lines(van.est.mat[1,])

plot(ton.logmapped.data$logmap_data[1,])
lines(ton.est.mat[1,])

plot(hal.logmapped.data$logmap_data[1,])
lines(hal.est.mat[1,])
par(mfrow=c(1,1))

#-----------------------------------------------------------------------
#exponential map 
van.expmap <- exp_map(datamat = combined.data$van.data$ang.data,
                      logmap_data =  van.est.mat, mean_vec = van.logmapped.data$mean_vec,
                      rotate.angle = van.logmapped.data$rotate.angle)

ton.expmap <- exp_map(datamat = combined.data$ton.data$ang.data,
                  logmap_data =  ton.est.mat, mean_vec = ton.logmapped.data$mean_vec,
                  rotate.angle = ton.logmapped.data$rotate.angle)

hal.expmap <- exp_map(datamat = combined.data$hal.data$ang.data,
                      logmap_data =  hal.est.mat, mean_vec = hal.logmapped.data$mean_vec,
                      rotate.angle = hal.logmapped.data$rotate.angle)
#---------------------------------------------------------------------------
#check estimation 
par(mfrow=c(2,3))
plot.index = 1
plot(van.logmapped.data$logmap_data[plot.index,])
lines(van.est.mat[plot.index,])

plot(ton.logmapped.data$logmap_data[plot.index,])
lines(ton.est.mat[plot.index,])

plot(hal.logmapped.data$logmap_data[plot.index,])
lines(hal.est.mat[plot.index,])

plot(seq(0,23,1), combined.data$van.data$ang.data[plot.index,] * 10,col='blue')
points(seq(0,23,1),van.expmap[plot.index,]*180 / pi,type='l',lwd=2)

plot(seq(0,23,1), combined.data$ton.data$ang.data[plot.index,] * 10,col='blue')
points(seq(0,23,1),ton.expmap[plot.index,]*180 / pi,type='l',lwd=2)

plot(seq(0,23,1), combined.data$hal.data$ang.data[plot.index,] * 10,col='blue')
points(seq(0,23,1),hal.expmap[plot.index,]*180 / pi,type='l',lwd=2)
par(mfrow=c(1,1))
