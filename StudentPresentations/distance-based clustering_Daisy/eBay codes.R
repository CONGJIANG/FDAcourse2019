####04-13-08
####This file contains the code used for the following paper 
####Title: "Distance-based clustering of sparsely observed stochastic processes, with applications to online auctions"
####Authors: Jie Peng and Hans Mueller
####Institute: University of California, Davis

########################### The following code is executed in R
source("functions.txt")    ### Load the R functions                      
library(sm)                ### Load R package "sm"

########## PART I: data processing
####read WTP values: data downloaded from the url given in the paper is a .xls file; It is then saved
####as a .txt file for easy input of R 
wtp<-read.csv("WTP.csv", header=TRUE) 

####get relative bidding time (in sections) to the auction starting time 
##date
year.bid<-matrix(0,nrow(wtp),2)
month.bid<-matrix(0,nrow(wtp),2)  
day.bid<-matrix(0,nrow(wtp),2)
day.diff<-numeric(nrow(wtp))   ##day difference
month.len<-c(31,28,31,30,31,30,31,31,30,31,30,31)

for(i in 1:nrow(wtp)){
  bid.time<-strptime(wtp[i,3],format="%b-%d-%Y")
  bid.start<-strptime(wtp[i,5],format="%b-%d-%Y")
  year.bid[i,1]<-strftime(bid.time,format="%Y")
  year.bid[i,2]<-strftime(bid.start,format="%Y")
  
  month.bid[i,1]<-strftime(bid.time,format="%m")
  month.bid[i,2]<-strftime(bid.start,format="%m")
  
  day.bid[i,1]<-strftime(bid.time,format="%d")
  day.bid[i,2]<-strftime(bid.start,format="%d")
  
  if (month.bid[i,1]==month.bid[i,2]){
    day.diff[i]<-as.numeric(day.bid[i,1])-as.numeric(day.bid[i,2])
  }else{
    
    cur<-as.numeric(month.bid[i,2])
    day.diff[i]<-month.len[cur]-as.numeric(day.bid[i,2])+as.numeric(day.bid[i,1])   
  }
  
}

all(year.bid[,1]==year.bid[,2])   ##all same year? yes
all((as.numeric(month.bid[,1])-as.numeric(month.bid[,2]))<=1)   ##all same month? no, but month difference is at most one


##time
time.bid<-matrix(0,nrow(wtp),2)  ##seconds to the 00:00:00 of the corresponding day

for (i in 1:nrow(wtp)){
  bid.time<-strptime(wtp[i,4],format="%H:%M:%S")
  bid.start<-strptime(wtp[i,6],format="%H:%M:%S")
  
  temp1<-as.numeric(strftime(bid.time,format="%H"))*3600
  temp2<-as.numeric(strftime(bid.time,format="%M"))*60
  temp3<-as.numeric(strftime(bid.time,format="%S"))
  time.bid[i,1]<-temp1+temp2+temp3
  
  
  temp1<-as.numeric(strftime(bid.start,format="%H"))*3600
  temp2<-as.numeric(strftime(bid.start,format="%M"))*60
  temp3<-as.numeric(strftime(bid.start,format="%S"))
  time.bid[i,2]<-temp1+temp2+temp3
  
}

##relative bidding time in seconds
rel.BidT<-numeric(nrow(wtp))
for(i in 1:nrow(wtp)){
  rel.BidT[i]<-day.diff[i]*3600*24+time.bid[i,1]-time.bid[i,2]
}



####win bid or not: the first bidding with the closing price is the winner
##
auc.ID<-unique(wtp[,1])  ##auction ID 
win<-numeric(nrow(wtp))
for(i in 1:length(auc.ID)){
  index.cur<-wtp[,1]==auc.ID[i]
  cur<-wtp[index.cur,2]
  cur.max<-cur[length(cur)]
  temp<-numeric(sum(index.cur))
  ind<-1
  count<-1
  while(ind){
    if(cur[count]==cur.max){
      temp[count]<-1
      ind<-0
    }else{
      count<-count+1
    }
  }
  
  win[index.cur]<-temp
}

##manully fix for the 95th auction where the opening and closing prices are the same 
i<-95
index.cur<-wtp[,1]==auc.ID[i]
temp<-numeric(sum(index.cur))
temp[2]<-1
win[index.cur]<-temp

##checks
sum(win)==length(auc.ID)  ##TRUE

##check whether the bidding time is monotonic increasing

time<-numeric(length(auc.ID))

for(i in 1:length(auc.ID)){
  index.cur<-wtp[,1]==auc.ID[i]
  cur<-rel.BidT[index.cur]
  time[i]<-all(diff(cur)>=0)
}
all(time)   ##TRUE


####combine the processed information with wtp data and save as  R dataset
wtp.new<-cbind(wtp, rel.BidT,win)
save(wtp.new, file="wtp.Rdata")



############# PART II: FPCA model fitting
load("wtp.Rdata")   ##load the wtp.Rdata got from PART I

####give new ID's for each bidder in each auction
ID.all<-numeric(nrow(wtp.new))
auc.ID<-unique(wtp.new[,1])

count<-1

for (i in 1:length(auc.ID)){
  index.cur<-wtp.new[,1]==auc.ID[i]
  bid.cur<-wtp.new[index.cur,7]
  bid.ID.cur<-unique(bid.cur)
  id.cur<-numeric(length(bid.cur))
  
  for (j in 1:length(bid.ID.cur)){
    id.cur[bid.cur==bid.ID.cur[j]]<-count
    count<-count+1
  }
  ID.all[index.cur]<-id.cur
  
}

####get rid of the opening bidding set by the seller, that is the one with a bidder ID "*"
index<-wtp.new[,7]!="*"
ebay.sub<-wtp.new[index,]
ID.all<-ID.all[index]
auc.ID.sub<-ebay.sub[,1]
ebay.sub[,1]<-ID.all
auc.ID<-unique(auc.ID.sub)


#### data summary 
length(unique(ebay.sub[,1]))   ##number of bidder trajectories: 1818
length(unique(ebay.sub[,7]))   ##number of bidders: 1122
length(auc.ID)                 ##number of auctions: 157; Note that one auction with all bidder ID being "*" is removed
nrow(ebay.sub)                 ##total number of bids: 3643


####format as Obs--measurement (bids) matrices, T--measurement (bidding) time, 
####and N--number of measurements (bids), each row for one subject

##
ID<-as.numeric(unique(ebay.sub[,1]))   #unique ID
num.auc<-length(ID)
t.max<-max(table(ebay.sub[,1]))   ##maximum number of measurements for one subject (a bidder trajectory)
Obs<-matrix(0,num.auc,t.max)
T<-matrix(0,num.auc,t.max)
N<-numeric(num.auc)
bidder.ID<-numeric(num.auc)     ##for each bidder trajectory
auc.ID.new<-numeric(num.auc)    ##for each bidder trajectory
ebay.sub<-as.matrix(ebay.sub) ##


##
for (i in 1:num.auc){
  print(i)
  ebay.cur<-matrix(ebay.sub[as.numeric(ebay.sub[,1])==ID[i],],ncol=ncol(ebay.sub))
  
  N[i]<-nrow(ebay.cur)
  Obs[i,1:N[i]]<-as.numeric(ebay.cur[,2])
  T[i,1:N[i]]<-as.numeric(ebay.cur[,9])/3600
  bidder.ID[i]<-ebay.cur[,7][1]
  auc.ID.new[i]<-auc.ID.sub[as.numeric(ebay.sub[,1])==ID[i]][1]
}


##summary: number bids of  each bidder trajectory
table(N)
#N
#   1    2    3    4    5    6    7    8    9   10   11   13   14   16   20   22 
#1046  361  179   91   51   28   23   18    9    2    4    1    1    2    1    1 

##check!
sum(N)==nrow(ebay.sub)      ##TRUE
length(unique(auc.ID.new))==length(auc.ID)     
length(unique(bidder.ID))==length(unique(ebay.sub[,7]))
num.auc==length(unique(ebay.sub[,1]))

#### FPCA model fitting
##(a)remove outliers
##(i) cross validation
data<-TranMtoV(Obs,T,N)
y<-data[,2]
t<-data[,3]
hmu.cv<-h.select(t,y,method="cv")
fitmu<-sm.regression(t,y,h=hmu.cv,poly.index=1,eval.points=t)$estimate 

##(ii)remove outliers
removal<-RemOut(Obs,T,N,fitmu,thre=3)   ###threshold for outliers: three standard deviations
Obs.r<-removal[[1]]
T.r<-removal[[2]]
N.r<-removal[[3]]
id.r<-ID[removal[[4]]]
bidder.ID.r<-bidder.ID[removal[[4]]]
auc.ID.r<-auc.ID.new[removal[[4]]]

## (iii) summary after removing outliers
table(N.r)
#N.r
#   1    2    3    4    5    6    7    8    9   10   11   13   14   16   20   22 
#1041  357  177   87   50   28   22   18    9    2    4    1    1    2    1    1 

length(unique(id.r))         ##number of remaining bidder trajectories: 1801
length(unique(bidder.ID.r))  ##number of remaining bidders:1112
length(unique(auc.ID.r))     ##number of remaining aucitons: 157
sum(N.r)                     ##number of remaining bids: 3596


##(iv) winning information after removing outliers
n.r<-nrow(Obs.r)
win.r<-numeric(n.r)
price.r<-numeric(n.r)

for (i in 1:n.r){
  index.w<-as.numeric(ebay.sub[,1])==id.r[i]
  l<-sum(index.w)
  win.r[i]<-any(as.numeric(ebay.sub[index.w,10])==1)
  price.r[i]<-as.numeric(ebay.sub[index.w,2][l])  ##close price
}

##check!
sum(as.numeric(ebay.sub[,10]))  ## true
sum(win.r)==length(unique(auc.ID.r)) ##false. 155 vs. 157, two auctions do not have a winner after the removal of outliers
sum(as.numeric(ebay.sub[,10]))  ##true



##(v) fit FPCA model for the data after removing outliers
##re-format
data.r<-TranMtoV(Obs.r,T.r,N.r)
y.r<-data.r[,2]
t.r<-data.r[,3]

##range of the bidding times and the evaluation grid
L1<-min(t.r)
L2<-max(t.r)
nump<-100
grid.age<-seq(L1,L2,by=(L2-L1)/nump)


##bandwidth selection: use h.select and use method="cv"
selec<-1
if(selec==1){
  #(i) mu
  hmu.cv<-h.select(t.r,y.r,method="cv")
  fitmu<-sm.regression(t.r,y.r,h=hmu.cv,poly.index=1,eval.points=t.r)$estimate 
  
  #G(s,t)I(s!=t)
  conver<-HatGOff(fitmu,data.r,N.r)  ##convert to off diagonal pair forma
  Cova<-conver[[2]]
  CovT<-conver[[3]]
  hcov.cv<-h.select(t(CovT),Cova,method="cv") 
  
  #G(t,t)+\sigma^2
  hsig.cv<-h.select(t.r,(y.r-fitmu)^2,method="cv")
}

##FPCA model fitting with the selected bandwidth
##
lin.result<-LocLinEst.new(Obs.r,T.r,N.r,grid.age,hmu.cv,hcov.cv,hsig.cv)
muhat<-lin.result[[1]]            ##estimated mean curve
covmatrix<-lin.result[[2]]        ##estimated covariance kernel
sig2hat<-lin.result[[4]]/(L2-L1)  ##estimated error variance

##estimate eigenfunctions and eigenvalues by projection: use linear interpolation: G(s,t)I(s!=t) 
eigenl<-EigenC(covmatrix,grid.age)
eigenfl<-eigenl[[1]]
eigenvl<-eigenl[[2]]

K<-5         ##dimension of the process: explain 96%
eigenv.est<-eigenvl[1:K]
eigenv.percent<-eigenvl[1:K]/sum(eigenvl)

eigenv.est
eigenv.percent
sum(eigenv.percent)  ##96%
sig2hat     ##697.4015

## plots
#(i)mu curve and and line connected individual trajectory: (x)
subject<-"wtp_whole"
y.lim<-c(min(as.vector(Obs.r)[as.vector(Obs.r)>0]),max(as.vector(Obs.r)[as.vector(Obs.r)>0]))
x.lim<-c(min(as.vector(T.r)[as.vector(T.r)>0]),max(as.vector(T.r)[as.vector(T.r)>0]))
plot(T.r[1,1:N.r[1]],Obs.r[1,1:N.r[1]],ylim=y.lim,xlim=x.lim,type='n',xlab="Time after auction started (in hours)",
     ylab="Bid price (in dollars)")
for(i in 1:n.r){
  cur<-i
  points(T.r[cur,1:N.r[cur]],Obs.r[cur,1:N.r[cur]],type='o',col=i)
}
points(grid.age,muhat,col='2',type='l', lwd=3)        
dev.print(postscript, file=paste(subject,'.trajectory_mu.ps',sep=""),height=8, width=10)


################# PART III: calculate conditional L_2 distance
####(a)conditional PC score 
PcScr<-FPcScore(Obs.r,T.r,N.r,grid.age,muhat,covmatrix,sig2hat,K) ##estimated conditional score 

####(b) pairwise conditional squared L2 distance
cdist<-CDist.mtr(Obs.r,T.r,N.r,grid.age,muhat,covmatrix,sig2hat,K)  ##estimated cdist
hist(cdist)
write(t(cdist),file="cdist.txt",ncolumns=ncol(cdist))  ##export for matlab analysis


#####################################The following code is executed in matlab
#########################PART IV: nonclassical multidimension scaling by matlab
%%use function mdscale in matlab to do multidimension scaling on the conditional distances cdist
%%use default path setting
%%read in the squared conditional distances
load cdist.txt;
p=2;  % dimension

%non-classic multidimension scaling by mdscale: only metricsstress converges
[Ycdstress,stress]=mdscale(sqrt(cdist),p,'criterion','metricstress');
[Ycdsstress,stress]=mdscale(sqrt(cdist),p,'criterion','metricsstress');
[Ycdsammon,stress]=mdscale(sqrt(cdist),p,'criterion','sammon');

%output the result for further  R analysis
%save ncmd.cdstress.txt Ycdstress -ascii;   %%not converge
save ncmd.cdsstress.txt Ycdsstress -ascii; 
%save ncmd.cdsammon.txt Ycdsammon -ascii;   %%not converge

########################################The following code is executed in R

##########################PART FIVE: kmean clustering of the mdscale result
md.sstress<-read.table("ncmd.cdsstress.txt")   ##load the MDS result by matlab
plot(md.sstress[,1],md.sstress[,2],xlab="dim 1", ylab="dim 2")

####clustering by kmeans
set.seed(1)
kmean.result<-kmeans(md.sstress, 6, iter.max = 5000, nstart =100)$cluster
plot(md.sstress[,1],md.sstress[,2],type='n',xlab="dim 1", ylab="dim 2")
text(md.sstress[,1],md.sstress[,2],kmean.result,cex=0.8,col=1)
kmean.result.ad<-kmean.result


##color and label and plot
index<-1:n.r
index4<-(kmean.result.ad==1)
grp4<-index[index4]
index3<-(kmean.result.ad==6)
grp3<-index[index3]
index1<-(kmean.result.ad==3)
grp1<-index[index1]
index5<-(kmean.result.ad==4)
grp5<-index[index5]
index2<-(kmean.result.ad==5)
grp2<-index[index2]
index6<-(kmean.result.ad==2)
grp6<-index[index6]

label<-kmean.result.ad
label[grp1]<-"E"
label[grp2]<-"L"
label[grp3]<-"H"
label[grp4]<-"A"
label[grp5]<-"S"
label[grp6]<-"F"

par(mfrow=c(1,1))
plot(md.sstress[,1],md.sstress[,2],type='n',xlab="dim 1", ylab="dim 2")
text(md.sstress[grp1,1],md.sstress[grp1,2],label[grp1],cex=1,col=5)
text(md.sstress[grp2,1],md.sstress[grp2,2],label[grp2],cex=1,col=2)
text(md.sstress[grp3,1],md.sstress[grp3,2],label[grp3],cex=1,col=3)
text(md.sstress[grp4,1],md.sstress[grp4,2],label[grp4],cex=1,col=4)
text(md.sstress[grp5,1],md.sstress[grp5,2],label[grp5],cex=1,col=6)
text(md.sstress[grp6,1],md.sstress[grp6,2],label[grp6],cex=1,col=rgb(1,0.5,0))
dev.print(postscript,file=paste(subject,"sstress.kmeans.ps",sep=""),height=8,width=10)


####fit mean curves for each cluster
#
Obs1<-Obs.r[grp1,]
T1<-T.r[grp1,]
N1<-N.r[grp1]
data<-TranMtoV(Obs1,T1,N1)
y1<-data[,2]
t1<-data[,3]
grid.age1<-seq(min(t1),max(t1),by=(max(t1)-min(t1))/nump)
hmu.cv<-h.select(t1,y1,method="cv")
hmu.cv<-10   ##manually adjust the smoothness a little bit for better visualization
fitmu1<-sm.regression(t1,y1,h=hmu.cv,poly.index=1,eval.points=grid.age1)$estimate 

#
Obs2<-Obs.r[grp2,]
T2<-T.r[grp2,]
N2<-N.r[grp2]
data<-TranMtoV(Obs2,T2,N2)
y2<-data[,2]
t2<-data[,3]
grid.age2<-seq(min(t2),max(t2),by=(max(t2)-min(t2))/nump)
hmu.cv<-h.select(t2,y2,method="cv")
hmu.cv<-12
fitmu2<-sm.regression(t2,y2,h=hmu.cv,poly.index=1,eval.points=grid.age2)$estimate 


#
Obs3<-Obs.r[grp3,]
T3<-T.r[grp3,]
N3<-N.r[grp3]
data<-TranMtoV(Obs3,T3,N3)
y3<-data[,2]
t3<-data[,3]
grid.age3<-seq(min(t3),max(t3),by=(max(t3)-min(t3))/nump)
hmu.cv<-h.select(t3,y3,method="cv")
fitmu3<-sm.regression(t3,y3,h=hmu.cv,poly.index=1,eval.points=grid.age3)$estimate 

#
Obs4<-Obs.r[grp4,]
T4<-T.r[grp4,]
N4<-N.r[grp4]
data<-TranMtoV(Obs4,T4,N4)
y4<-data[,2]
t4<-data[,3]
grid.age4<-seq(min(t4),max(t4),by=(max(t4)-min(t4))/nump)
hmu.cv<-h.select(t4,y4,method="cv")
hmu.cv<-12
fitmu4<-sm.regression(t4,y4,h=hmu.cv,poly.index=1,eval.points=grid.age4)$estimate 


##
Obs5<-Obs.r[grp5,]
T5<-T.r[grp5,]
N5<-N.r[grp5]
data<-TranMtoV(Obs5,T5,N5)
y5<-data[,2]
t5<-data[,3]
grid.age5<-seq(min(t5),max(t5),by=(max(t5)-min(t5))/nump)
hmu.cv<-h.select(t5,y5,method="cv")
fitmu5<-sm.regression(t5,y5,h=hmu.cv,poly.index=1,eval.points=grid.age5)$estimate 


##
Obs6<-Obs.r[grp6,]
T6<-T.r[grp6,]
N6<-N.r[grp6]
data<-TranMtoV(Obs6,T6,N6)
y6<-data[,2]
t6<-data[,3]
grid.age6<-seq(min(t6),max(t6),by=(max(t6)-min(t6))/nump)
hmu.cv<-h.select(t6,y6,method="cv")
fitmu6<-sm.regression(t6,y6,h=hmu.cv,poly.index=1,eval.points=grid.age5)$estimate 




##mean curves for each cluster
muhat.a<-c(fitmu1,fitmu2,fitmu3,fitmu4,fitmu5,fitmu6)
par(mfrow=c(3,2))
#
plot(grid.age2,fitmu2,type='l',xlab='Time after auction started (in hours)',ylab='Bid price (in dollars)',ylim=c(min(y.r,muhat.a),max(y.r,muhat.a)),xlim=range (t),main="Mean bid trajectory of cluster L",cex=1,lwd=2)        
points(t2,y2,type='p')

#
plot(grid.age3,fitmu3,type='l',xlab='Time after auction started (in hours)',ylab='Bid price (in dollars)',ylim=c(min(y.r,muhat.a),max(y.r,muhat.a)),xlim=range (t),main="Mean bid trajectory of cluster H",cex=1,lwd=2)        
points(t3,y3,type='p')

#
plot(grid.age5,fitmu5,type='l',xlab='Time after auction started (in hours)',ylab='Bid price (in dollars)' ,ylim=c(min(y.r,muhat.a),max(y.r,muhat.a)),main="Mean bid trajectory of cluster S",cex=1,lwd=2)        
points(t5,y5,type='p')

#
plot(grid.age6,fitmu6,type='l',xlab='Time after auction started (in hours)',ylab='Bid price (in dollars)',ylim=c(min(y.r,muhat.a),max(y.r,muhat.a)),xlim=range (t),main="Mean bid trajectory of cluster F",cex=1,lwd=2)        
points(t6,y6,type='p')

#
plot(grid.age4,fitmu4,type='l',xlab='Time after auction started (in hours)',ylab='Bid price (in dollars)',ylim=c(min(y.r,muhat.a),max(y.r,muhat.a)),xlim=range (t),main="Mean bid trajectory of cluster A",cex=1,lwd=2)        
points(t4,y4,type='p')

#
plot(grid.age1,fitmu1,type='l',xlab='Time after auction started (in hours)',ylab='Bid price (in dollars)',ylim=c(min(y.r,muhat.a),max(y.r,muhat.a)),xlim=range(t),main="Mean bid trajectory of cluster E", cex=1,lwd=2)        
points(t1,y1,type='p')
par(mfrow=c(1,1))
dev.print(postscript,file=paste(subject,".mucurves.groups.ps",sep=""),height=10,width=10)




##individual trajectory for each cluster
y.lim<-c(min(as.vector(Obs.r)),max(as.vector(Obs.r)))
x.lim<-c(min(as.vector(T.r)[as.vector(T.r)>0]),max(as.vector(T.r)[as.vector(T.r)>0]))
par(mfrow=c(3,2))
#
plot(T2[1,1:N2[1]],Obs2[1,1:N2[1]],ylim=y.lim,xlim=x.lim,type='n',xlab="Time after auction started (in hours)",
     ylab="Bid price (in dollars)",main="Individual trajectories of cluster L")

for(i in 1:length(grp2)){
  cur<-i
  points(T2[cur,1:N2[cur]],Obs2[cur,1:N2[cur]],type='o',col=i,pch=22)
  
}

#
plot(T3[1,1:N3[1]],Obs3[1,1:N3[1]],ylim=y.lim,xlim=x.lim,type='n',xlab="Time after auction started (in hours)",
     ylab="Bid price (in dollars)",main="Individual trajectories of cluster  H")
for(i in 1:length(grp3)){
  cur<-i
  points(T3[cur,1:N3[cur]],Obs3[cur,1:N3[cur]],type='o',col=i,pch=22)
  
}




#
plot(T5[1,1:N5[1]],Obs5[1,1:N5[1]],ylim=y.lim,xlim=x.lim,type='n',xlab="Time after auction started (in hours)",
     ylab="Bid price (in dollars)",main="Individual trajectories of cluster  S")
for(i in 1:length(grp5)){
  cur<-i 
  points(T5[cur,1:N5[cur]],Obs5[cur,1:N5[cur]],type='o',col=i,pch=22)
  
}



#
plot(T6[1,1:N4[1]],Obs6[1,1:N4[1]],ylim=y.lim,xlim=x.lim,type='n',xlab="Time after auction started (in hours)",
     ylab="Bid price (in dollars)",main="Individual trajectories of cluster  F")

for(i in 1:length(grp6)){
  cur<-i 
  points(T6[cur,1:N6[cur]],Obs6[cur,1:N6[cur]],type='o',col=i,pch=22) 
}


#
plot(T4[1,1:N4[1]],Obs4[1,1:N4[1]],ylim=y.lim,xlim=x.lim,type='n',xlab="Time after auction started (in hours)",
     ylab="Bid price (in dollars)",main="Individual trajectories of cluster  A")

for(i in 1:length(grp4)){
  cur<-i 
  points(T4[cur,1:N4[cur]],Obs4[cur,1:N4[cur]],type='o',col=i,pch=22) 
}


##
plot(T1[1,1:N1[1]],Obs1[1,1:N1[1]],ylim=y.lim,xlim=x.lim,type='n',xlab="Time after auction started (in hours)",
     ylab="Bid price (in dollars)",main="Individual trajectories  of cluster E")
for(i in 1:length(grp1)){
  cur<-i
  points(T1[cur,1:N1[cur]],Obs1[cur,1:N1[cur]],type='o',col=i,pch=22) 
}


par(mfrow=c(1,1))
dev.print(postscript,file=paste(subject,".indtraject.groups.ps",sep=""),height=10,width=10)




####range of bidding time for each bidder
bid.range<-numeric(n.r)
for (i in 1:n.r){
  bid.range[i]<-T.r[i,N.r[i]]-T.r[i,1]
}

range1<-bid.range[grp1]
range2<-bid.range[grp2]
range3<-bid.range[grp3]
range4<-bid.range[grp4]
range5<-bid.range[grp5]
range6<-bid.range[grp6]


###bidding summary
c(length(grp1),range(t1), mean(N1), range(y1))
c(length(grp2),range(t2), mean(N2), range(y2))
c(length(grp3),range(t3), mean(N3), range(y3))
c(length(grp4),range(t4), mean(N4), range(y4))
c(length(grp5),range(t5), mean(N5), range(y5))
c(length(grp6),range(t6), mean(N6), range(y6))



####winning rate and mean price paid on winning
c(mean(win.r[grp1]),sum(win.r[grp1]),mean(price.r[index1&win.r==1]),mean(price.r[grp1]))
c(mean(win.r[grp2]),sum(win.r[grp2]),mean(price.r[index2&win.r==1]),mean(price.r[grp2]))
c(mean(win.r[grp3]),sum(win.r[grp3]),mean(price.r[index3&win.r==1]),mean(price.r[grp3]))
c(mean(win.r[grp4]),sum(win.r[grp4]),mean(price.r[index4&win.r==1]),mean(price.r[grp4]))
c(mean(win.r[grp5]),sum(win.r[grp5]),mean(price.r[index5&win.r==1]),mean(price.r[grp5]))
c(mean(win.r[grp6]),sum(win.r[grp6]),mean(price.r[index6&win.r==1]),mean(price.r[grp6]))


####histogram of bidding times

par(mfrow=c(3,2))

hist(t2,xlim=range(t),xlab="Time after auction started (in hours)",main="Cluster L",breaks=50,cex.lab=1.5,cex.main=1.5)
hist(t3,xlim=range(t),xlab="Time after auction started (in hours)",main="Cluster H",breaks=50,cex.lab=1.5,cex.main=1.5)
hist(t5,xlim=range(t),xlab="Time after auction started (in hours)",main="Cluster S",breaks=50,cex.lab=1.5,cex.main=1.5)
hist(t6,xlim=range(t),xlab="Time after auction started (in hours)",main="Cluster F",breaks=50,cex.lab=1.5,cex.main=1.5)
hist(t4,xlim=range(t),xlab="Time after auction started (in hours)",main="Cluster A",breaks=50,cex.lab=1.5,cex.main=1.5)
hist(t1,xlim=range(t),xlab="Time after auction started (in hours)",main="Cluster E",breaks=50,cex.lab=1.5,cex.main=1.5)
par(mfrow=c(1,1))
dev.print(postscript,file=paste(subject,".histtime.groups.ps",sep=""),height=8,width=12)



#### histogram of the amount of bids

par(mfrow=c(3,2))
hist(y2,xlim=range(y.r),xlab="Bid amount (in dollars)",main="Cluster L",breaks=50,cex.lab=1.5,cex.main=1.5)
hist(y3,xlim=range(y.r),xlab="Bid amount (in dollars)",main="Cluster H",breaks=50,cex.lab=1.5,cex.main=1.5) 

hist(y5,xlim=range(y.r),xlab="Bid amount (in dollars)",main="Cluster S",breaks=50,cex.lab=1.5,cex.main=1.5)
hist(y6,xlim=range(y.r),xlab="Bid amount (in dollars)",main="Cluster F",breaks=50,cex.lab=1.5,cex.main=1.5)

hist(y4,xlim=range(y.r),xlab="Bid amount (in dollars)",main="Cluster A",breaks=50,cex.lab=1.5,cex.main=1.5)
hist(y1,xlim=range(y.r),xlab="Bid amount (in dollars)",main="Cluster E",breaks=50,cex.lab=1.5,cex.main=1.5)
par(mfrow=c(1,1))
dev.print(postscript,file=paste(subject,".histprice.groups.ps",sep=""),height=8,width=12)


####result
> 
  > ###bidding summary
  > c(length(grp1),range(t1), mean(N1), range(y1))
[1] 221.0000000   0.2969444 167.9950000   1.8461538   0.0100000 235.0000000
> c(length(grp2),range(t2), mean(N2), range(y2))
[1] 386.0000000   0.6041667 167.9997222   1.8031088   5.0000000 244.5000000
> c(length(grp3),range(t3), mean(N3), range(y3))
[1] 399.000000   1.367500 167.999167   1.759398  20.000000 280.000000
> c(length(grp4),range(t4), mean(N4), range(y4))
[1] 285.0000000   0.1813889 167.9994444   2.1859649  12.5000000 283.5000000
> c(length(grp5),range(t5), mean(N5), range(y5))
[1] 298.0000000   0.2472222 167.9988889   2.5134228   0.0100000 218.0000000
> c(length(grp6),range(t6), mean(N6), range(y6))
[1] 212.0000000   0.2363889 167.9786111   1.9716981  10.0000000 230.0000000
> 
  > 
  > 
  > ##winning and mean price paid on winning
  > c(mean(win.r[grp1]),sum(win.r[grp1]),mean(price.r[index1&win.r==1]),mean(price.r[grp1]))
[1] 4.524887e-03 1.000000e+00 2.300000e+02 3.972905e+01
> c(mean(win.r[grp2]),sum(win.r[grp2]),mean(price.r[index2&win.r==1]),mean(price.r[grp2]))
[1]   0.07512953  29.00000000 208.12344828 163.92927461
> c(mean(win.r[grp3]),sum(win.r[grp3]),mean(price.r[index3&win.r==1]),mean(price.r[grp3]))
[1]   0.2431078  97.0000000 237.0779381 215.7275689
> c(mean(win.r[grp4]),sum(win.r[grp4]),mean(price.r[index4&win.r==1]),mean(price.r[grp4]))
[1]   0.09122807  26.00000000 245.28807692 194.14056140
> c(mean(win.r[grp5]),sum(win.r[grp5]),mean(price.r[index5&win.r==1]),mean(price.r[grp5]))
[1]   0.00671141   2.00000000 177.25000000  88.55167785
> c(mean(win.r[grp6]),sum(win.r[grp6]),mean(price.r[index6&win.r==1]),mean(price.r[grp6]))
[1]   0.0000   0.0000      NaN 111.3988