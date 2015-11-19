rm(list=ls(all=T))
library(SIS)
set.seed(1000)
nob<-1000
dm<-10000
xx<-matrix(rnorm(dm*nob),nr=nob,nc=dm)
mu<-matrix(rnorm(nob),nr=nob)
beta<-c(seq(1,20),rep(0,dm-20))
my<-xx%*%beta+mu
fit1<-SIS(xx,my,penalty="lasso")
#fit0<-SIS(xx,my,penalty="lasso",alpha=0.5)

aa<-1
library(glmnet)
library(grpreg)
#fit2<-glmnet(xx,my,alpha=0)
mx<-colMeans(xx)
vx<-apply(xx,2,var)
sx<-sweep(xx,2,apply(xx,2,mean))
sx<-xx%*%diag(1/sqrt(vx))
colnames(sx)<-seq(1,ncol(sx))
colnames(xx)<-seq(1,ncol(xx))

sxy<-abs(t(sx)%*%my)
summary(sxy)
id0<-which(sxy>sort(sxy,decreasing=T)[50])
r1<-grpreg(xx[,id0],my,alpha=aa)
sr1<-select(r1,"BIC")
res1<-my-cbind(1,xx[,id0])%*%c(sr1$beta)
bb<-sr1$beta
nid0<-(names(bb)[which(bb!=0)])
nid0<-as.numeric(nid0[-1])

sxy<-abs(t(sx[,-nid0])%*%res1)
summary(sxy)
id1<-as.numeric(colnames(sx[,-nid0])[which(sxy>sort(sxy,decreasing=T)[50])])
r2<-grpreg(xx[,id1],res1,alpha=aa)
sr2<-select(r2,"BIC")
bb<-sr2$beta
nid1<-(names(bb)[which(bb!=0)])
nid1<-as.numeric(nid1[-1])
scid<-sort(c(nid0,nid1))

nr2<-grpreg(xx[,scid],my,alpha=aa)
nsr2<-select(nr2,"BIC")
bb<-nsr2$beta
nid1<-(names(bb)[which(bb!=0)])
nid1<-as.numeric(nid1[-1])
res1<-my-cbind(1,xx[,scid])%*%c(bb)
scid<-nid1

sxy<-abs(t(sx[,-scid])%*%res1)
summary(sxy)
id1<-as.numeric(colnames(sx[,-scid])[which(sxy>sort(sxy,decreasing=T)[50])])
r2<-grpreg(xx[,id1],res1,alpha=aa)
#sr2<-select(r2,"BIC",df.method="active")
sr2<-select(r2,"BIC")
bb<-sr2$beta
nid1<-(names(bb)[which(bb!=0)])
nid1<-as.numeric(nid1[-1])
scid<-sort(c(scid,nid1))




