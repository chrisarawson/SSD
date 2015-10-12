#fits an SSD to tox data
#plots with CIs generated using a bootsrap function
#Note that the protection values are calculated on the curve fitted to the existing data
#Bootstrapping the curve will alter the mean values so is only for providing
#indicative CIs

#Currently the final SSD plot has the fitted data (black) and bootstrapped (red)
#These converge as the CIs contract - generally as more data is added
#SSDs with low n tend to be further from bootstrapped solutions so are 
#more suspect.

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
  #generate random values in the distribution 
  #(same number as originally used)
  xr <- rlnorm(fit$n,meanlog=fit$estimate[1],
               sdlog=fit$estimate[2])
  #fit distribution to new data
  fitr<-fitdistr(xr,"lognormal")
  #find the quantiles from the new fitted distrbution
  #and return the value at probability p
  hc5r<-qlnorm(p,meanlog=fitr$estimate[1],
               sdlog=fitr$estimate[2])
  return(hc5r)
}

#Repeat this function a large number of times (1000)
#Each repeat gives a value for the value for the pth percentile
#So can generate at list of the values of that percentile
#get the quantiles of the bootstrapped values
#This gives the average (50th percentile), and CIS (2.5th and 97.5th percentiles)

set.seed(-123)
hc1_booted<-replicate(1000,myboot(fit,p=0.01))
hc1_booted_CIS<-quantile(hc1_booted,probs=c(0.025,0.5,0.975))

hc5_booted<-replicate(1000,myboot(fit,p=0.05))
hc5_booted_CIS<-quantile(hc5_booted,probs=c(0.025,0.5,0.975))

hc10_booted<-replicate(1000,myboot(fit,p=0.1))
hc10_booted_CIS<-quantile(hc10_booted,probs=c(0.025,0.5,0.975))

#Tabulate as protection levels
protValTab<-as.data.frame(rbind(hc10_booted_CIS,hc5_booted_CIS,hc1_booted_CIS))
colnames(protValTab)<-c("2.5% CI","Central","97.5% CI")
row.names(protValTab)<-c("99%", "95%", "90%")

#Now to plot this
#first we need the data for the curve

myboot2<-function(fit2,newxs) {
  #first resample again
  #generate a random set of points along the curve
  xr<-rlnorm(fit$n,meanlog=fit2$estimate[1],
           sdlog=fit2$estimate[2])
  #fit a log normal distribution to those random points
  fitr<-fitdistr(xr,"lognormal")
  #calculate the y-value for each x in newxs based on this fit
  pyr<-plnorm(newxs,meanlog=fitr$estimate[1],
              sdlog=fitr$estimate[2])
  return(pyr)
}

#create a list of x values along which to calculate the curve
newxs <- 10^(seq(log10(0.01),log10(max(a$val)+1000000),
                 length.out=1000))
#Run the above function 1000 times. Each time will add a set of y values for each
#x in newxs - so 1000 sets of data.
boots<-replicate(1000,myboot2(fit,newxs))

#extract the mean and CIs to plot
#this asks for the 2.5 and 97.5 percentile for each x across all 1000 datasets in BOOTS
cis<-apply(boots,1,quantile,c(0.025,0.5,0.975))
rownames(cis)<-c("lwr","mean","upr")

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
a$Fitted <- 10^(log10(qlnorm(a$prop, meanlog = fit$estimate[1], 
                             sdlog = fit$estimate[2])) -0.4)

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
lines(pdat$newxs,pdat$mean,type="l",col="red")

legend(0.1,1,c("Eggs","Larvae","Juveniles","Adults","Fitted","Bootstrapped"),
       pch=c(0,1,2,6, NA, NA),lty=c(NA,NA,NA,NA,1,1),
       col=c("black", "black","black","black","black","red"),bty="n")

#create a list of outputs and return
mySumm<-list(DataFit=fit,FittedPVs=hcs,BootstrappedPVs=protValTab)

#Create a generc function to print some sensible results
class(mySumm)<-"SSDSummary"
print.SSDSummary<-function(x,...) {
  cat("SSD Summary","\n\n")
  cat("Log-normal parameters","\n")
  print(x$DataFit)
  cat("\n\n","Fitted Protection Values","\n")
  cat("90%", x$FittedPVs[1],"\n")
  cat("95%", x$FittedPVs[2],"\n")
  cat("99%", x$FittedPVs[3],"\n\n")
  cat("Bootstrapped protection Values","\n")
  print(x$BootstrappedPVs,...)
}

#Return results
mySumm