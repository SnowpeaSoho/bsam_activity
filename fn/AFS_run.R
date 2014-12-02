
require(rjags)
library(coda)
library(markovchain)
library(sp)

load("data/AFS_project.RData")

source("fn/check.divetime.tstep.r") # function deals with the HOs
source("fn/dat4jags_reg.r") # adapted to regular start times
source("fn/plot.output.r")
source("fn/ssm4.r")

tod <- TRUE
tstep_hr <- 6
tstep <- tstep_hr/24

# prepare data for jags
use <- dat4jags(smru[,1:5],tstep=tstep,tod=tod)
# prepare 'ho' component
dive.check <- check.divetime.tstep(data=use,HO=summ,use.binary=TRUE,dives=dives)
dat1 <- dive.check[[1]][[1]] # start with FM11-S (short) ; then try long FM07-S
dat1$id <- factor(dat1$id,levels=unique(dat1$id))
dat1 <- list(dat1)

n.chains <- 2
n.iter <- 2000#0#0
n.thin <- 10
n.burn <- 5000#0#0

st = proc.time()
fit <- ssm(loc.list=dat1, model="haulout4", LM=FALSE, HO=TRUE,
           adapt=n.burn, samples=n.iter, thin=n.thin, chains=n.chains)

cat("Elapsed time: ", round((proc.time() - st)[3]/60,2), "min \n")

plot.output(fit,dat1,haulout4=TRUE)
