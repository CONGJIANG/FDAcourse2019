install.packages("fdadensity")

library(fdadensity)

x <- seq(-5,5,length.out =512)
y <- dnorm(x)
q = pnorm(x)

plot(x,y,type="l",main="probability density function")

#dens2lqd is a function for converting Densities to Quantile Densities 
y.lqd <- dens2lqd( dens=y, dSup = x) 

## dens2quantile() is a function for converting Densities to Quantile Functions

y.quantile <- dens2quantile(dens=y, dSup = x) 
plot(q,y.quantile,type="l",main="quantile function")

plot(q,y.lqd,type="l", ylim=c(0,5),main = "log quantile density function")





data(Top50BabyNames)

# This data set consist of two list element: Top50BabyNames$x, and Top50BabyNames$dens
# *Param* Top50BabyNames$x:  grid of years between 1950 and 2016, of length 67.

# *Param* Top50BabyNames$dens: list of length two, corresponding to male (dens$male) and female(dens$female) names.
#         Each is a 67-by-50 matrix of density estimates, where each column corresponds to a unique
#         baby name given by the corresponding column name.
# It's a density function of the distribution of the popularity index of a certain name. 

densitysample = t(Top50BabyNames$dens$female)
densitygrid =  Top50BabyNames$x

# each row of the density sample sum up to 1
sum(densitysample[1,])
sum(densitysample[2,])
sum(densitysample[3,])

plot(densitygrid,densitysample[1,],type="l")
lines(densitygrid,densitysample[2,],type="l")
lines(densitygrid,densitysample[3,],type="l")

# Perform Transformation FPCA for male baby name densities


fpcaObj = FPCAdens(dmatrix = densitysample, dSup = densitygrid, useAlpha = FALSE,
optns = list(dataType = 'Dense', error = FALSE, methodSelectK = 2))


#It fails because one density value is zero
densitysample[33,]

#useAlpha set it the zero or negative density value to 0.01
fpcaObj = FPCAdens(dmatrix = densitysample, dSup = densitygrid, useAlpha = TRUE,
optns = list(dataType = 'Dense', error = FALSE, methodSelectK = 3))
plot(fpcaObj)


#we may look at the source code for better understanding
FPCAdens

# Note that the resulting FPCA object is based on the log quantile density space, 
# we need to map it back using function lqd2dens()
 
 # suppose I want to map the mean function back
 nu = fpcaObj$mu
 denMu <- lqd2dens(nu, dSup = densitygrid, useSplines = TRUE)
 plot(densitygrid, denMu)
 
 # mapping the mode of variation back
 k = 1
 alpha = c(-0.2,-0.1,0,0.1,0.2)
 obsGrid = fpcaObj$obsGrid
 s = fpcaObj$workGrid
 nu = fpcaObj$mu
 sigma = sqrt(fpcaObj$lambda[k])
 rho = fpcaObj$phi[, k]
 Qmatrix = alpha %*% t(rho * sigma) + matrix(rep(nu, length(alpha)), nrow = length(alpha), byrow = TRUE)
 DENreg <- t(apply(Qmatrix, 1, function(u) lqd2dens(u, dSup = densitygrid, useSplines = TRUE)))
 
 plot(densitygrid,DENreg[1,])
 plot(densitygrid,DENreg[2,])
 plot(densitygrid,DENreg[3,])
 plot(densitygrid,DENreg[4,])

 ## You may also use the function CreateModeOfVarPlotLQ2D() to do it 
 # the variation in FPC1 is horizontal
 CreateModeOfVarPlotLQ2D(fpcaObj, k = 1, dSup = x, Qvec = alpha, main = 'First Mode, Density Space')
 # the variation in FPC2 is vertical
 CreateModeOfVarPlotLQ2D(fpcaObj, k = 2, dSup = x, Qvec = alpha, main = 'Second Mode, Density Space')
 # the variation in FPC3 relates to the heavy-tailedness
 CreateModeOfVarPlotLQ2D(fpcaObj, k = 3, dSup = x, Qvec = alpha, main = 'Third Mode, Density Space')
 
 

 

			