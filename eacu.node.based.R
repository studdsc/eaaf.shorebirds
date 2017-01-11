
# packages
library(lattice)
library(coda)
library(R2WinBUGS)
library(R2jags)

# data
eacu.nmix.node <- read.csv("/Users/studdsc/Box Sync/shorebird models/eastern.curlew/data/eacu.nmix.csv", header = TRUE)
eacu.nmix.node <- subset(eacu.nmix.node, select=c(site,node,year,count1,count2),period > 1992)
eacu.nmix.node <- within(eacu.nmix.node,{
  count1 <- as.integer(count1)
  count2 <- as.integer(count2)
  })

y <- array(NA, dim = c(18, 2, 20))	

for(k in 1:20){
   sel.rows <- with(eacu.nmix.node,year == k)
   y[,,k] <- as.matrix(eacu.nmix.node)[sel.rows, 4:5]
   }
y; str(y)

eacu.cov <- read.csv("/Users/studdsc/Box Sync/shorebird models/eastern.curlew/data/eacu.covariates.csv", header = TRUE)
site <- with(eacu.cov[1:18,],as.vector(site))
node<- with(eacu.cov[1:18,],as.vector(node))
year <- with(eacu.cov[1:20,],as.vector(scale(year,scale=FALSE)))
nsites=nrow(y) 
nyears=20
nnodes=3 
ncounts=ncol(y)

# model
ptm=proc.time()
eacu.nmix.node<-function(){

#nested indexing
for (s in 1:nsites){
	ste[s] ~ dnorm(mu.site, tau.site)                   
	}
for (b in 1:nnodes){
	alpha[b] ~ dnorm(mu.alpha, tau.alpha) 
	trend[b] ~ dnorm (mu.trend, tau.trend)
	}
for (k in 1:nyears){
  beta[k] ~ dnorm(0, 0.001)
	}

# likelihood
for (i in 1:nsites){                                 
for (k in 1:nyears){
  eps[i,k] ~ dnorm(0, tau.eps)                            
	N[i,k] ~ dpois(lambda[i,k])               
	log(lambda[i,k]) <-  ste[site[i]] + alpha[node[i]] + trend[node[i]] * year[k] + eps[i,k]

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
   node1.N[k] <- sum(N[1:3,k]) 
   node2.N[k] <- sum(N[4:11,k])
   node3.N[k] <- sum(N[12:18,k]) 
   detection[k] <- mean(site.year.p[,k])
   p.res.1[k] <- mean(p.res[1:3,k]) 
   p.res.2[k] <- mean(p.res[4:11,k])
   p.res.3[k] <- mean(p.res[12:18,k])
   }
  
trendif12 <- trend[1]-trend[2]
trendif13 <- trend[1]-trend[3]
trendif23 <- trend[2]-trend[3]
or1v2 <- (trend[1]/(1-(trend[1]))+1) / (trend[2]/(1-(trend[2]))+1)
or1v3 <- (trend[1]/(1-(trend[1]))+1) / (trend[3]/(1-(trend[3]))+1)
or2v3 <- (trend[2]/(1-(trend[2]))+1) / (trend[3]/(1-(trend[3]))+1)
mean.detection <- mean(detection[])
fit.actual <- sum(chi.sq.actual[,,])
fit.sim <- sum(chi.sq.sim[,,])
c.hat <- fit.actual/fit.sim
bpv <- step(fit.sim-fit.actual)

# priors  
tau.site ~ dgamma(0.001,0.001)
tau.alpha ~ dgamma(0.001,0.001)
tau.trend ~ dgamma(0.001,0.001)
tau.eps ~ dgamma(0.001,0.001)
tau.p ~ dgamma(0.001,0.001)
sd.site <- 1 / sqrt(tau.site)
sd.alpha <- 1/sqrt(tau.alpha)
sd.trend <- 1/sqrt(tau.trend) 
sd.eps <- 1/sqrt(tau.eps) 
sd.p <- 1/sqrt(tau.p) 
mu.site ~ dnorm(0,0.001) 
mu.alpha ~ dnorm(0,0.001)
mu.trend ~ dnorm(0,0.001)
}

# bundle data
jags.data <- list("y", "nsites", "nyears", "nnodes", "ncounts", "site", "node", "year")

# initial values
N<-function(x){
  Nst <- apply(x, c(1, 3), max) + 1
  Nst[is.na(Nst)] <- 4000
  return(Nst)
  }

inits <- function()list(N = N(y), tau.site=1, tau.alpha=1, tau.trend=1, tau.eps=1, tau.p=1)

# run chains in parallel
source("/Users/studdsc/Box Sync/shorebird models/bar.tailed.godwit/scripts/jags.parallel2.r")
library(abind)

inits.parallel<-list()
for(i in 1:n.chains){
  inits.parallel[[i]]<-inits()
  }
str(inits.parallel)

# parameters monitored
params <- c("node1.N","node2.N","node3.N","trend","trendif12","trendif13","trendif23","or1v2",
"or1v3","or2v3","detection","mean.detection","fit.actual","fit.sim","c.hat","bpv","p.res.1",
"p.res.2","p.res.3")

# set random seed
set.seed(52)

# call JAGS from R
eacu.nmix.node <- jags.parallel2(data=jags.data,model=eacu.nmix.node,inits=inits.parallel,parameters=params, 
n.chains = 3, n.thin = 18, n.iter = 800000, n.burnin = 200000, 
                                 
working.directory = getwd())

# summarize posteriors
print(eacu.nmix.node, dig = 3)
proc.time()-ptm