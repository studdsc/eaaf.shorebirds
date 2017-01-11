
# packages
library(lattice)
library(coda)
library(R2WinBUGS)
library(R2jags)

# data
eacu.nmix.cont <- read.csv("/Users/studdsc/Box Sync/shorebird models/eastern.curlew/data/eacu.nmix.csv", header = TRUE)
eacu.nmix.cont <- subset(eacu.nmix.cont, select=c(site,year,count1,count2),period > 1992)
eacu.nmix.cont <- within(eacu.nmix.cont,{
  count1 <- as.integer(count1)
  count2 <- as.integer(count2)
  })

y <- array(NA, dim = c(18, 2, 20))	

for(k in 1:20){
  sel.rows <- with(eacu.nmix.cont,year == k)
  y[,,k] <- as.matrix(eacu.nmix.cont)[sel.rows, 3:4]
  }
y; str(y)

eacu.cov <- read.csv("/Users/studdsc/Box Sync/shorebird models/eastern.curlew/data/eacu.covariates.csv", header = TRUE)
site <- with(eacu.cov[1:18,],as.vector(site))
year <- with(eacu.cov[1:20,],as.vector(scale(year,scale=FALSE)))
nsites=nrow(y) 
nyears=20
ncounts=ncol(y)

# model
ptm=proc.time()
eacu.nmix.cont<-function(){ 
  
#nested indexing
for (s in 1:nsites){
  ste[s] ~ dnorm(mu.site, tau.site)                   
  }
for (k in 1:nyears){
  beta[k] ~ dnorm(0, 0.001)
  }
  
# likelihood
for (i in 1:nsites){                                 
for (k in 1:nyears){
  eps[i,k] ~ dnorm(0, tau.eps)                            
  N[i,k] ~ dpois(lambda[i,k])               
  log(lambda[i,k]) <-  ste[site[i]] + trend * year[k] + eps[i,k]

# detection probability
for (j in 1:ncounts){
  y[i,j,k] ~ dbin(p[i,j,k], N[i,k])      
  p[i,j,k] <- 1 / (1 + exp(-logit.p[i,j,k])) 
  logit.p[i,j,k] ~ dnorm(beta[k], tau.p)      
	
# goodness of fit 	
  expect[i,j,k] <- p[i,j,k] * N[i,k]
  chi.sq.actual[i,j,k] <- pow((y[i,j,k] - expect[i,j,k]),2) / (expect[i,j,k]+0.5)
  sim[i,j,k] ~ dbin(p[i,j,k], N[i,k])
  chi.sq.sim[i,j,k] <- pow((sim[i,j,k] - expect[i,j,k]),2) / (expect[i,j,k]+0.5)
  }
#	sites and year detection probability 
  site.year.p[i,k] <- mean(p[i,,k])
# residuals
  p.res[i,k] <- (N[i,k]-lambda[i,k])/sqrt(lambda[i,k]) 
  }
}
  
# dervied parameters
for (k in 1:nyears){
  totalN[k] <- sum(N[,k])
  totalL[k] <- sum(lambda[,k])
  detection[k] <- mean(site.year.p[,k])
  p.res.m[k] <- mean(p.res[,k])
  }
  
mean.detection <- mean(detection[])
fit.actual <- sum(chi.sq.actual[,,])
fit.sim <- sum(chi.sq.sim[,,])
c.hat <- fit.actual/fit.sim
bpv <- step(fit.sim-fit.actual)  
  
# priors
tau.site ~ dgamma(0.001,0.001)
tau.trend ~ dgamma(0.001,0.001)
tau.eps ~ dgamma(0.001,0.001)
tau.p ~ dgamma(0.001,0.001)
sd.site <- 1/sqrt(tau.site)
sd.trend <- 1/sqrt(tau.trend) 
sd.eps <- 1/sqrt(tau.eps) 
sd.p <- 1/sqrt(tau.p) 
mu.site ~ dnorm(0,0.001)
trend ~ dnorm(mu.trend, tau.trend)
mu.trend ~ dnorm(0,0.001)
}

# bundle data
jags.data <- list("y", "nsites", "nyears", "ncounts", "site", "year")
  
# initial values
  N<-function(x){
  Nst <- apply(x, c(1, 3), max) + 1
  Nst[is.na(Nst)] <- 4000
  return(Nst)
  }
  
inits <- function()list(N = N(y), tau.site=1, tau.trend=1, tau.eps=1, tau.p=1)

# run chains in parallel
source("/Users/studdsc/Box Sync/shorebird models/bar.tailed.godwit/scripts/jags.parallel2.r")
library(abind)

inits.parallel<-list()
for(i in 1:n.chains){
  inits.parallel[[i]]<-inits()
  }
str(inits.parallel)

# Parameters monitored
params <- c("totalN","totalL","trend","detection","mean.detection","fit.actual","fit.sim",
"c.hat","bpv","p.res.m")

# set random seed
set.seed(52)

# call JAGS from R
eacu.nmix.cont <- 
jags.parallel2(data=jags.data,model=eacu.nmix.cont,inits=inits.parallel,parameters=params, 
n.chains = 3, n.thin = 18, n.iter = 800000, n.burnin = 200000, 
                                 
working.directory = getwd())

# summarize posteriors
print(eacu.nmix.cont, dig = 3)
proc.time()-ptm



traceplot(eacu.nmix.cont)

eacu.nmix.cont.mc <- as.mcmc(eacu.nmix.cont)
plot(eacu.nmix.cont.mc, ask=TRUE)
summary(eacu.nmix.cont.mc, dig=3)
autocorr.plot(eacu.nmix.cont.mc, ask=TRUE)

# evaluate model fit
fit <- eacu.nmix.cont$BUGSoutput$sims.list$fit
fit.new <- eacu.nmix.cont$BUGSoutput$sims.list$fit.new
plot(eacu.nmix.cont$BUGSoutput$sims.list$fit, eacu.nmix.cont$BUGSoutput$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE)
abline(0, 1, lwd = 2, col = "black")
mean(eacu.nmix.cont$BUGSoutput$sims.list$fit.new > eacu.nmix.cont$BUGSoutput$sims.list$fit)
mean(eacu.nmix.cont$BUGSoutput$mean$fit) / mean(eacu.nmix.cont$BUGSoutput$mean$fit.new)

# total abundance counts
max.day.count <- apply(y, c(1, 3), max, na.rm = TRUE)
max.day.count[max.day.count == "-Inf"] <- NA
sum.max.count <- apply(max.day.count, 2, sum, na.rm = TRUE)
sum.max.count[1:20]

#total N
year<-1993:2012
p.res.m.fitted<-apply(eacu.nmix.cont$BUGSoutput$sims.list$p.res.m,2,mean)
n.fitted<-apply(eacu.nmix.cont$BUGSoutput$sims.list$totalN,2,mean)
n.lower.CI<-apply(eacu.nmix.cont$BUGSoutput$sims.list$totalN,2,quantile,probs=0.025)
n.upper.CI <-apply(eacu.nmix.cont$BUGSoutput$sims.list$totalN,2,quantile,probs=0.975)
min.N <- min(c(n.fitted, n.lower.CI), na.rm = TRUE)
max.N <- max(c(n.fitted, n.upper.CI), na.rm = TRUE)
plot(year,n.fitted,type="l",xlim=c(1992,2013),ylim=c(min.N,max.N),
xlab="",ylab="",cex=1.2,col="black",pch=21)
lines(year,n.lower.CI,col="black",lty=2)
lines(year,n.upper.CI,col="black",lty=2)

plot(n.fitted,p.res.m.fitted)

#posterior plot data frame
eacu.nmix.plot.cont2 <- data.frame(cbind(sum.max.count,n.fitted,n.lower.CI,n.upper.CI,p.fitted,
p.lower.CI,p.upper.CI))
write.csv(eacu.nmix.plot.cont2,file="eacu.nmix.plot.cont2.csv",append=FALSE,row.names=TRUE)

#posterior density
trend.den <- as.vector(eacu.nmix.cont$BUGSoutput$sims.list$trend)
totalN.den <- as.vector(eacu.nmix.cont$BUGSoutput$sims.list$totalN[,20])
hist(trend.den)
hist(totalN.den)

par(mfrow=c(1,2))
plot(density(trend.den),col="red")
abline(b=6,v=mean(trend.den))
plot(density(totalN.den))
abline(b=6,v=mean(totalN.den),col="red")

#posterior density data frame
eacu.nmix.density.cont2 <- data.frame(cbind(trend.den,totalN.den,fit,fit.new))
write.csv(eacu.nmix.density.cont2,file="eacu.nmix.density.cont2.csv",append=FALSE,row.names=TRUE)

#detection probability
year<-1993:2012
p.fitted <- apply(eacu.nmix.cont$BUGSoutput$sims.list$detection,2,mean)
p.lower.CI<-apply(eacu.nmix.cont$BUGSoutput$sims.list$detection,2,quantile,probs=0.025)
p.upper.CI <-apply(eacu.nmix.cont$BUGSoutput$sims.list$detection,2,quantile,probs=0.975)
plot(year,p.fitted,ylab="Detection probabilty",xlab="",type="b",ylim=c(0,1),pch=16,frame.plot = FALSE)
segments(year,p.lower.CI,year,p.upper.CI)
