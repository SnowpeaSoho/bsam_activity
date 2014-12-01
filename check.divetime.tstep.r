# IDs SSM location dates encompassed entirely within a HO period
# NB usually very low for SES (and often at very tail)
# check.haulouts: inputs are 'data' from dat4bugCOV and 'HO'=direct from SMRU table
check.divetime.tstep <- function(data,HO,use.binary=TRUE, dives=NULL)  {
HO_cnt <- array(NA,dim=c(length(data),4))
keep.tally <- list()

for (i in 1:length(data)) {

tstep <- data[[i]]$tstep
dates <- data[[i]]$first.date + 86400*tstep*c(0:c(data[[i]]$RegN-1))

data[[i]]$ho <- matrix(0,nrow=data[[i]]$RegN,ncol=1)
ho <- HO[HO$ref==names(data)[i] ,]
dv <- dives[dives$ref==names(data)[i] ,]

tmp <-merge(data.frame(dates=dates,ref=names(data)[i]), ho[,c("gmt","DIVE_TM")],
  by.x="dates",by.y="gmt",all.x=TRUE,all.y=FALSE)

# want ANY diving to result in the p(non-forage) being < 0.5
# scale so p(non-forage) becomes lowest when highest diving recorded
tmp$index <-  0.5-0.5*tmp$DIVE_TM/max(tmp$DIVE_TM,na.rm=TRUE)
tmp$index[tmp$DIVE_TM==0] <- 1 # no diving == 1
tmp$index2 <- ifelse(tmp$DIVE_TM==0,1,0)   # simple binary answer

if (!is.null(dives))     {
# have passed in individual dive data, may have info missed by summary file
    no_info <- which(is.na(tmp$DIVE_TM) | tmp$DIVE_TM==0) # actually summaries can also miss diving seen in individual obs
    # assign individual dives to fall within a timestep
    dv$tstep <- findInterval(dv$gmt,tmp$dates,right=TRUE,all.inside=FALSE)
    updated <-  no_info[no_info %in% dv$tstep]
    tmp$index[updated] <- runif(length(updated),0,0.5)  # random uniform number distn b/w 0-0.5
    tmp$index2[updated] <- 0 } # assign as 'active'
    
# really have no information here; just set at 0.5 (coin toss)
tmp$index[is.na(tmp$index)] <- 0.5
tmp$index2[is.na(tmp$index2)] <- 0.5       

if(use.binary) {data[[i]]$ho <- as.vector(tmp$index2) }
else { data[[i]]$ho <-  as.vector(tmp$index)  }         

keep.tally[[i]] <- tmp

#data[[i]]$HO_cnt <- c( length(which(data[[i]]$haulouts>=1)) , # full HO
#  length(which(data[[i]]$haulouts>0 & data[[i]]$haulouts<1)), # partial HO
#  length(data[[i]]$haulouts) )   # number timesteps
HO_cnt[i,] <-  c( length(which(data[[i]]$ho>=1)) , # full HO
  length(which(data[[i]]$ho>=0 & data[[i]]$ho<0.5)), # diving
  length(which(data[[i]]$ho==0.5)), # no info available
  length(data[[i]]$ho) )   # number timesteps
} # end tagID loop 
HO_cnt <-  as.data.frame(HO_cnt)
colnames(HO_cnt) <- c("inactive","active","no.info","total")
rownames(HO_cnt) <- names(data)

return(list(data=data, HO_cnt=HO_cnt,keep.tally)  )

} # end function



check.individ.dives <- function(tmp, dives) {
  is.na(tmp$index)

}



