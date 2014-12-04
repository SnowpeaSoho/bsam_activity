
require(rjags)
library(coda)
library(markovchain)
library(sp)

load("data/AFS_project.RData")

source("fn/check.divetime.tstep.r") # function deals with the HOs
source("fn/dat4jags_reg.r") # adapted to regular start times
source("fn/plot.output.r")
source("fn/ssm4.r")

x.mode <- function(x) {mod.x <- as.numeric(names(abs(sort(-table(x)))[1]))} # mode only
x.mode.pr <- function(x) { # mode and also proportion of state 2 estimates across samples
  x.out <- array(0,dim=c(1,2))
  tab <- abs(sort(-table(x)))
  lab <- as.numeric(names(tab))
  pr.2 <- as.numeric(tab[which(lab==2)])  
  x.out[1,1] <- lab[1]; if (length(pr.2)>0) { x.out[1,2] <- pr.2}
  x.out
}

tod <- TRUE
tstep_hr <- 6
tstep <- tstep_hr/24
n.chains <- 2
n.iter <- 2000
n.thin <- 1
n.burn <- 5000

# prepare data for jags
use <- dat4jags(smru[,1:5],tstep=tstep,tod=tod)
# prepare 'ho' component
dive.check <- check.divetime.tstep(data=use,HO=summ,use.binary=TRUE,dives=dives)

for (itag in 2:length(use)) {
dat1 <- dive.check[[1]][[itag]] # start with FM11-S (short) ; then try long FM07-S
dat1$id <- factor(dat1$id,levels=unique(dat1$id))
dat1 <- list(dat1)

st = proc.time()
fit <- ssm(loc.list=dat1, model="haulout4", LM=FALSE, HO=TRUE,
           adapt=n.burn, samples=n.iter, thin=n.thin, chains=n.chains)

cat("Elapsed time: ", round((proc.time() - st)[3]/60,2), "min \n")

tmp <- t(apply(fit[[1]]$mcmc$b,1,x.mode.pr)) # not very quick, but ok
fit[[1]]$summary$b.mode <- tmp[,1] # record state mode
fit[[1]]$summary$b2.prop <- tmp[,2]/prod(dim(fit[[1]]$mcmc$b)[2:3]) # recordproportion of samples that are s2

plot.output(fit,dat1,haulout4=TRUE,map_range=NULL); #c(75,86,-67,-63))

table(fit[[1]]$summary$b.5,fit[[1]]$summary$b.mode)

length(which(fit[[1]]$summary$b>2))

# new way of looking at state summaries, given multinomial output
b.new <- rep(2, length=nrow(fit[[1]]$summary))
b.new[fit[[1]]$summary$b>2] <- 3 # anything above 2 gets allocated as 3
b.new[fit[[1]]$summary$b.mode==1] <- 1 # highest support for s1
b.new[b.new==2 & fit[[1]]$summary$b2.prop < 0.5] <- 3 # poor support for s2 (i.e result of averaging s1 & s3 samples)

table(dat1[[1]]$ho==1,b.new) # TRUE & 2 should be low
table(dat1[[1]]$ho>=0.5,b.new) # FALSE & 3 should be 0
fit[[1]]$summary$b.new <- b.new
save(fit,dat1,file=paste("output/",dat1[[1]]$id[1],"haulout4.RData",sep=""))

}
