plot.output <- function(fit, dat1,map_range=NULL,haulout4=FALSE) {

out <- fit[[1]]$summary

table(out$b.5)
table(dat1[[1]]$ho>0.5)

if(haulout4) {
  out$b.new <- 2
  out$b.new[out$b>2] <- 3
  out$b.new[out$b<1.5] <- 1
} else {
  out$b.new<- out$b.5
}

if(is.null(map_range)) {

 map_range <- c(range(out$lon)+c(-0.5,0.5) , range(out$lat)+c(-0.5,0.5)  )
}

#windows()
pdf(file=paste(dat1[[1]]$id[1],"haulout3.pdf",sep=""))

#windows() # MAP STATES
par(mar=c(1,1,1.5,1), oma=c(2,2,0,0))
plot(argosdat$lon[argosdat$id==as.character(dat1[[1]]$id[1])],argosdat$lat[argosdat$id==as.character(dat1[[1]]$id[1])],
     #xlim=c(78,82.5),ylim=c(-68.5,-66.5),
     xlim=map_range[1:2],ylim=map_range[3:4],
     type="b",col="gray",lty=3,
     xlab="",ylab="" ,yaxt="n",xaxt="n",main="", mgp=c(3,0.25,0))
coastplot()
points(out$lon,out$lat,col=c("black"),type="b")
points(out$lon,out$lat,col=c("blue","red","green")[out$b.new],pch=16)
points(out$lon[dat1[[1]]$ho==1],out$lat[dat1[[1]]$ho==1],col=c("black"))

degAxis(1, mgp=c(3,0.25,0),tcl=-0.25);
degAxis(2,at=c(seq(-68,-60,by=2)), mgp=c(3,0.25,0),tcl=-0.25) ;
mtext(dat1[[1]]$id[1],3,0.25,cex=1.5) #,adj=0.05
mtext("Latitude",2,2); mtext("Longitude",1,2)


# STATE TIME SERIES
plot(out$date,out$b,type="b",main=dat1[[1]]$id[1])
points(out$date,out$b,col=c("blue","red","green")[out$b.new],pch=16)

#windows()
# GAMMA and THETA hist
par(mfrow=c(3,2),mar=c(2.5,2.5,1,1));
layout(matrix(1:6,3,2,byrow=FALSE))#, widths=c(3,3,3,3), heights=rep(1,5))
hist(fit[[1]]$mcmc$gamma[1,,],40,xlim=c(0,1),main="Gamma1")
hist(fit[[1]]$mcmc$gamma[2,,],40,xlim=c(0,1),main="Gamma2")
hist(fit[[1]]$mcmc$gamma[3,,],40,xlim=c(0,1),main="Gamma3")

hist(fit[[1]]$mcmc$theta[1,,],40,xlim=c(-1,2*pi),main="Theta1")
hist(fit[[1]]$mcmc$theta[2,,],40,xlim=c(-1,2*pi),main="Theta2")
hist(fit[[1]]$mcmc$theta[3,,],40,xlim=c(-1,2*pi),main="Theta3")
range(fit[[1]]$mcmc$theta)


#windows() # CHAINS
plot(as.mcmc.list(fit[[1]]$mcmc$gamma),main="Gamma")
#windows()
plot(as.mcmc.list(fit[[1]]$mcmc$theta),main="Theta")

plot(as.mcmc.list(fit[[1]]$mcmc$Sigma),main="Sigma")

dev.off()
}