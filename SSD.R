#fits an SSD to tox data
#plots with CIs generated using a bootsrap function

#read data: "val" = value used in SSD, "Species" =Species

a<-read.csv(file.choose(),header=TRUE)

#create a column with the ranked proportion by val
a <- a[order(a$val),]
a$prop <- ppoints(a$val,0.5)

#quick labelled plot of the datapoints - check for doubles etc
plot(a$val,a$prop,xlab="Concentration (mg/L)",
     ylab="Proportion Species Affected",log="x",xlim=c(1,10000))
text(a$val,a$prop,a$Species,cex=0.6,pos=2)

#fit a log normal distribution
library(MASS)
fit <- fitdistr(a$val,"lognormal")
fit

#derive the 90th 95th and 99th % protection level
hcs <- qlnorm(c(0.1,0.05,0.01),meanlog=fit$estimate[1],
              sdlog=fit$estimate[2])
hcs

#We should also get some CIs for that
#This uses a bootstrapping of the fitted distribution

#create a function that fits random values from the fitted
#distribution to the distribution and then estimate an HCp
#from the new distribution

myboot<-function(fit,p) {
  #resample from fitted distribution
  xr <- rlnorm(fit$n,meanlog=fit$estimate[1],
               sdlog=fit$estimate[2])
  #fit distribution to new data
  fitr<-fitdistr(xr,"lognormal")
  #return HCp
  hc5r<-qlnorm(p,meanlog=fitr$estimate[1],
               sdlog=fitr$estimate[2])
  return(hc5r)
}

#Repeat this function a large number of times (1000) and
#get the quantiles of the bootstrapped values

set.seed(1234)
hc1_booted<-replicate(1000,myboot(fit,p=0.01))
hc1_booted_CIS<-quantile(hc1_booted,probs=c(0.025,0.5,0.975))

hc5_booted<-replicate(1000,myboot(fit,p=0.05))
hc5_booted_CIS<-quantile(hc5_booted,probs=c(0.025,0.5,0.975))

hc10_booted<-replicate(1000,myboot(fit,p=0.1))
hc10_booted_CIS<-quantile(hc10_booted,probs=c(0.025,0.5,0.975))

#Tabulate as protection levels
protValTab<-as.data.frame(rbind(hc1_booted_CIS,hc5_booted_CIS,hc10_booted_CIS))
row.names(protValTab)<-c("99%", "95%", "90%")


#Return table of Protection Levels
protValTab

#Now to plot this

#first we need the data for the curve

myboot2<-function(fit,newxs) {
  #first resample again
  xr<-rlnorm(fit$n,meanlog=fit$estimate[1],
           sdlog=fit$estimate[2])
  #fit
  fitr<-fitdistr(xr,"lognormal")
  #predict for new data
  pyr<-plnorm(newxs,meanlog=fitr$estimate[1],
              sdlog=fitr$estimate[2])
  return(pyr)
}

#new data to predict
newxs <- 10^(seq(log10(0.01),log10(max(a$val)+1000000),
                 length.out=10000))
boots<-replicate(10000,myboot2(fit,newxs))

#now create the data frame for the plot
library(reshape2)
bootdat<-data.frame(boots)
bootdat$newxs<-newxs
bootdat<-melt(bootdat,id="newxs")

#extract the CIs to plot
cis<-apply(boots,1,quantile,c(0.025,0.975))
rownames(cis)<-c("lwr","upr")

#add the fitted values to a new dataframe
pdat<-data.frame(newxs,
                 py=plnorm(newxs,
                           meanlog=fit$estimate[1],
                                 sdlog=fit$estimate[2]))

#add CIs to the fitted values dataframe
pdat<-cbind(pdat,t(cis))

#Break into factors for simple plotting of groups
eggs<-a[which(a$Life.stage == "eggs"),]
larvae<-a[which(a$Life.stage == "larvae"),]
Juvenile<-a[which(a$Life.stage == "Juvenile"),]
adult<-a[which(a$Life.stage == "adult"),]

# x coordinates for species names
a$Fitted <- 10^(log10(qlnorm(a$prop, meanlog = fit$estimate[1], sdlog = fit$estimate[2])) -0.4)

#plot
plot(eggs$val,eggs$prop,
     xlab="Concentration (mg/L)",
     ylab="Proportion Species Affected",
     log="x",ylim=c(0,1),xlim=c(0.1,1000000),pch=0)
text(a$Fitted,a$prop,a$Species,pos=2,cex=0.5)
points(larvae$val,larvae$prop,pch=1)
points(Juvenile$val,Juvenile$prop,pch=2)
points(adult$val,adult$prop,pch=6)

lines(pdat$newxs,pdat$py,type="l")
lines(pdat$newxs,pdat$lwr,type="l",lty=2)
lines(pdat$newxs,pdat$upr,type="l",lty=2)

legend(0.1,1,c("Eggs","Larvae","Juveniles","Adults"),pch=c(0,1,2,6),bty="n")

