rm(list=ls(all=T))
library(SIS)
#set.seed(1000)
nob<-1000
dm<-10000
xx<-matrix(rnorm(dm*nob),nr=nob,nc=dm)
mu<-matrix(rnorm(nob),nr=nob)
beta<-c(seq(1,20),rep(0,dm-20))
my<-xx%*%beta+mu
#fit1<-SIS(xx,my,penalty="lasso")

aa<-1
ng<-200
library(glmnet)
library(grpreg)
#fit2<-glmnet(xx,my,alpha=0)
mx<-colMeans(xx)
vx<-apply(xx,2,var)
sx<-sweep(xx,2,apply(xx,2,mean))
sx<-xx%*%diag(1/sqrt(vx))
colnames(sx)<-seq(1,ncol(sx))
colnames(xx)<-seq(1,ncol(xx))

sxy<-abs(t(sx)%*%my)/as.numeric(sqrt(var(my)))
#sxy<-abs(cor(sx,my))
summary(sxy)
id0<-which(sxy>=sort(sxy,decreasing=T)[ng])
#r1<-grpreg(xx[,id0],my,alpha=aa,group.multiplier=1/sxy[id0])
r1<-grpreg(xx[,id0],my,alpha=aa)
sr1<-select(r1,"BIC",df.method="active")
res1<-my-cbind(1,xx[,id0])%*%c(sr1$beta)
bb<-sr1$beta
nid0<-(names(bb)[which(bb!=0)])
nid0<-as.numeric(nid0[-1])

sxy<-abs(t(sx[,-nid0])%*%res1)/as.numeric(sqrt(var(res1)))
#sxy<-abs(cor(sx,my))
summary(sxy)
id1<-as.numeric(colnames(sx[,-nid0])[which(sxy>sort(sxy,decreasing=T)[ng])])
#r2<-grpreg(xx[,id1],res1,alpha=aa,group.multiplier=1/sxy[id1])
r2<-grpreg(xx[,id1],res1,alpha=aa)
sr2<-select(r2,"BIC",df.method="active")
bb<-sr2$beta
nid1<-(names(bb)[which(bb!=0)])
nid1<-as.numeric(nid1[-1])
scid<-sort(c(nid0,nid1))

sxy<-abs(t(sx[,scid])%*%my)/as.numeric(sqrt(var(my)))
summary(sxy)
id0<-which(sxy>=sort(sxy,decreasing=T)[30])
#nr2<-grpreg(xx[,id0],my,alpha=aa,group.multiplier=1/sxy[id0])
nr2<-grpreg(xx[,id0],my,alpha=aa)
nsr2<-select(nr2,"BIC",df.method="active")
bb<-nsr2$beta
nid1<-(names(bb)[which(bb!=0)])
nid1<-as.numeric(nid1[-1])
res1<-my-cbind(1,xx[,id0])%*%c(bb)
scid<-nid1

sxy<-abs(t(sx[,-scid])%*%res1)/as.numeric(sqrt(var(res1)))
summary(sxy)
id1<-as.numeric(colnames(sx[,-scid])[which(sxy>sort(sxy,decreasing=T)[ng])])
#r2<-grpreg(xx[,id1],res1,alpha=aa,group.multiplier=1/sxy[id1])
r2<-grpreg(xx[,id1],res1,alpha=aa)
sr2<-select(r2,"BIC",df.method="active")
bb<-sr2$beta
nid1<-(names(bb)[which(bb!=0)])
nid1<-as.numeric(nid1[-1])
scid<-sort(c(scid,nid1))
#####

sxy<-abs(t(sx[,scid])%*%my)/as.numeric(sqrt(var(my)))
summary(sxy)
id0<-which(sxy>=sort(sxy,decreasing=T)[length(sxy)])
nr2<-grpreg(xx[,id0],my,alpha=aa)
nsr2<-select(nr2,"BIC",df.method="active")
bb<-nsr2$beta
nid1<-(names(bb)[which(bb!=0)])
nid1<-as.numeric(nid1[-1])
res1<-my-cbind(1,xx[,id0])%*%c(bb)
print(scid<-nid1)


