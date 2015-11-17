rm(list=ls(all=T))
library(SIS)
nob<-1000
dm<-10000
xx<-matrix(rnorm(dm*nob),nr=nob,nc=dm)
mu<-matrix(rnorm(nob),nr=nob)
beta<-c(seq(1,20),rep(0,dm-20))
my<-xx%*%beta+mu
fit1<-SIS(xx,my,penalty="lasso")
#fit0<-SIS(xx,my,penalty="lasso",alpha=0.5)

library(glmnet)
library(grpreg)
#fit2<-glmnet(xx,my,alpha=0)
mx<-colMeans(xx)
vx<-apply(xx,2,var)
sx<-sweep(xx,2,apply(xx,2,mean))
sx<-xx%*%diag(1/sqrt(vx))
colnames(sx)<-seq(1,ncol(sx))
sxy<-abs(t(sx)%*%my)
summary(sxy)
id0<-which(sxy>quantile(sxy,0.95))
r1<-grpreg(xx[,id0],my,alpha=0.5)
sr1<-select(r1,"BIC")
res1<-my-cbind(1,xx[,id0])%*%c(sr1$beta)


sxy<-abs(t(sx[,-id0])%*%res1)
summary(sxy)
id1<-as.numeric(colnames(sx[,-id0])[which(sxy>quantile(sxy,0.95))])
r2<-grpreg(xx[,id1],my,alpha=0.5)
sr2<-select(r2,"BIC")

scid<-sort(c(id0,id1))

fit0<-grpreg(xx[,scid],my,alpha=0.5)
sr0<-select(fit0,"BIC")

bfit0<-grpreg(xx,my,alpha=0.5)
b0<-select(bfit0,"BIC")

