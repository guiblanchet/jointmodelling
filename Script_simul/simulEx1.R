rm(list=ls())
library(HMSC)
#==========================
### Load the simulated data
#==========================
data("simulEx1")
data("simulParamEx1")

#==============
### Build model
#==============
model<-hmsc(simulEx1,family="probit",niter=10000,nburn=1000,thin=10)

#===============
### Plot results
#===============
#---------------------
## Plot results paramX
#---------------------
### True values
truth<-as.vector(simulParamEx1$param$paramX)

### Average
average<-apply(model$results$estimation$paramX,1:2,mean)

### 95% confidence intervals
CI.025<-apply(model$results$estimation$paramX,1:2,quantile, probs=0.025)
CI.975<-apply(model$results$estimation$paramX,1:2,quantile, probs=0.975)
CI<-cbind(as.vector(CI.025),as.vector(CI.975))

### Draw confidence interval plots
plot(0,0,xlim=c(1,nrow(CI)),ylim=range(CI,truth),type="n",xlab="",ylab="",main="paramX")
abline(h=0,col="grey")
arrows(x0=1:nrow(CI),x1=1:nrow(CI),y0=CI[,1],y1=CI[,2],code=3,angle=90,length=0.05)
points(1:nrow(CI),average,pch=15,cex=1.5)
points(1:nrow(CI),truth,col="red",pch=19)

### Mixing object
mixing<-as.mcmc(model,parameters="paramX")

### Draw trace and density plots for all combination of paramters
plot(mixing)

### Convert the mixing object to a matrix
mixingDF<-as.data.frame(mixing)

### Draw boxplot for each parameters
par(mar=c(7,4,4,2))
boxplot(mixingDF,las=2)

### Draw beanplots
library(beanplot)
par(mar=c(7,4,4,2))
beanplot(mixingDF,las=2)

#-----------------------------------------------------------
### Plot random effect estimation through correlation matrix
#-----------------------------------------------------------
corMat<-corRandomEff(model,cor=FALSE)

#________________________________
### First set of latent variables
#________________________________
### Isolate the values of interest
ltri<-lower.tri(apply(corMat[,,,1],1:2,quantile,probs=0.025),diag=TRUE)

### True values
truth<-as.vector(tcrossprod(simulParamEx1$param$paramLatent[[1]])[ltri])

### Average
average<-as.vector(apply(corMat[,,,1],1:2,mean)[ltri])

### 95% confidence intervals
corMat.025<-as.vector(apply(corMat[,,,1],1:2,quantile,probs=0.025)[ltri])
corMat.975<-as.vector(apply(corMat[,,,1],1:2,quantile,probs=0.975)[ltri])
CI<-cbind(corMat.025,corMat.975)

### Plot the results
plot(0,0,xlim=c(1,nrow(CI)),ylim=range(CI,truth),type="n",xlab="",,main="cov(paramLatent[[1,1]])")
abline(h=0,col="grey")
arrows(x0=1:nrow(CI),x1=1:nrow(CI),y0=CI[,1],y1=CI[,2],code=3,angle=90,length=0.05)
points(1:nrow(CI),average,pch=15,cex=1.5)
points(1:nrow(CI),truth,col="red",pch=19)

### Mixing object
mixing<-as.mcmc(model,parameters="paramLatent")

### Draw trace and density plots for all combination of paramters
plot(mixing[[1]])

### Convert the mixing object to a matrix
mixingDF<-as.data.frame(mixing[[1]])

### Draw boxplot for each parameters
par(mar=c(7,4,4,2))
boxplot(mixingDF,las=2)

### Draw beanplots
library(beanplot)
par(mar=c(7,4,4,2))
beanplot(mixingDF,las=2)

### Draw estimated correlation matrix
library(corrplot)
corMat<-corRandomEff(model,cor=TRUE)
averageCor<-apply(corMat[,,,1],1:2,mean)
corrplot(averageCor,method="color",col=colorRampPalette(c("blue","white","red"))(200))

### Draw chord diagram
library(circlize)
corMat<-corRandomEff(model,cor=TRUE)
averageCor<-apply(corMat[,,,1],1:2,mean)
colMat<-matrix(NA,nrow=nrow(averageCor),ncol=ncol(averageCor))
colMat[which(averageCor>0.4,arr.ind=TRUE)]<-"red"
colMat[which(averageCor< -0.4,arr.ind=TRUE)]<-"blue"
chordDiagram(averageCor,symmetric=TRUE,annotationTrack=c("name","grid"),grid.col="grey",col=colMat)

#_________________________________
### Second set of latent variables
#_________________________________
corMat<-corRandomEff(model,cor=FALSE)
### Isolate the values of interest
ltri<-lower.tri(apply(corMat[,,,2],1:2,quantile,probs=0.025),diag=TRUE)

### True values
truth<-as.vector(tcrossprod(simulParamEx1$param$paramLatent[[2]])[ltri])

### Average
average<-as.vector(apply(corMat[,,,2],1:2,mean)[ltri])

### 95% confidence intervals
corMat.025<-as.vector(apply(corMat[,,,2],1:2,quantile,probs=0.025)[ltri])
corMat.975<-as.vector(apply(corMat[,,,2],1:2,quantile,probs=0.975)[ltri])
CI<-cbind(corMat.025,corMat.975)

### Plot the results
plot(0,0,xlim=c(1,nrow(CI)),ylim=range(CI,truth),type="n",xlab="",main="cov(paramLatent[[1,2]])")
abline(h=0,col="grey")
arrows(x0=1:nrow(CI),x1=1:nrow(CI),y0=CI[,1],y1=CI[,2],code=3,angle=90,length=0.05)
points(1:nrow(CI),average,pch=15,cex=1.5)
points(1:nrow(CI),truth,col="red",pch=19)

### Mixing object
mixing<-as.mcmc(model,parameters="paramLatent")

### Draw trace and density plots for all combination of paramters
plot(mixing[[2]])

### Convert the mixing object to a matrix
mixingDF<-as.data.frame(mixing[[2]])

### Draw boxplot for each parameters
par(mar=c(7,4,4,2))
boxplot(mixingDF,las=2)

### Draw beanplots
library(beanplot)
par(mar=c(7,4,4,2))
beanplot(mixingDF,las=2)


### Draw estimated correlation matrix
library(corrplot)
corMat<-corRandomEff(model,cor=TRUE)
averageCor<-apply(corMat[,,,2],1:2,mean)
corrplot(averageCor,method="color",col=colorRampPalette(c("blue","white","red"))(200))

### Draw chord diagram
library(circlize)
corMat<-corRandomEff(model,cor=TRUE)
averageCor<-apply(corMat[,,,2],1:2,mean)
colMat<-matrix(NA,nrow=nrow(averageCor),ncol=ncol(averageCor))
colMat[which(averageCor>0.4,arr.ind=TRUE)]<-"red"
colMat[which(averageCor< -0.4,arr.ind=TRUE)]<-"blue"
chordDiagram(averageCor,symmetric=TRUE,annotationTrack=c("name","grid"),grid.col="grey",col=colMat)

