`ssm` =
function (loc.list, adapt, samples, thin, chains, LM=FALSE, HO=FALSE, res,
	buffer, ...)
{
  #browser()
    ssm1 = function(input) {

      y = input$y
        idx = input$idx
        RegN = input$RegN
        j = input$j
        x = rbind(y[idx[-length(idx)], ], y[max(idx) - 1, ])
        id = input$id

       	## create land mask if specified
		if(LM){
			print("Creating land mask...")
			data(worldLLhigh)
			x.min = floor(min(ifelse(y[,1]< 0, y[,1]+360, y[,1]), na.rm=TRUE)) - 1
			x.max = ceiling(max(ifelse(y[,1]< 0, y[,1]+360, y[,1]), na.rm=TRUE)) + 1
			y.min = floor(min(y[,2], na.rm=TRUE)) - 1
			y.max = ceiling(max(y[,2], na.rm=TRUE)) + 1

			lm = clipPolys(worldLLhigh, xlim=c(x.min,x.max), ylim=c(y.min,y.max))
			grid = expand.grid(x=seq(x.min,x.max,by=res), y=seq(y.min,y.max,by=res))
			ev = as.EventData(data.frame(EID=1:nrow(grid), X=grid[,1],Y=grid[,2]),
				projection="LL")
			f.ev = findPolys(ev,lm, maxRows=1e6)
			grid$water = 1
			grid$water[f.ev$EID] = 0
			water = matrix(grid$water, length(unique(grid$y)), length(unique(grid$x)),
				byrow=TRUE)
			print("Land mask created")
			}

        row.na = which(is.na(x[, 1]))
        x[row.na, 1] = approx(seq(nrow(x)), x[, 1], xout = row.na, rule = 2)$y
        row.na = which(is.na(x[, 2]))
        x[row.na, 2] = approx(seq(nrow(x)), x[, 2], xout = row.na, rule = 2)$y
        row.na = which(is.na(y[, 1]) | is.na(y[, 2]))
        if (length(row.na) == 0) {
            y.init = y
            y.init[!is.na(y.init)] = NA
        	}
        else {
            y.init = y
            y.init[!is.na(y.init)] = 9999
            yinit.xna = approx(seq(nrow(y)), y[, 1], xout = row.na, rule = 2)$y
            yinit.yna = approx(seq(nrow(y)), y[, 2], xout = row.na, rule = 2)$y
            y.init[row.na, ] = cbind(yinit.xna, yinit.yna)
            y.init[y.init == 9999] = NA
        	}

	    if(LM){
    		## make sure no x's are on land
			test.x = diag(water[trunc((x[,2]-y.min+res)/res),
				trunc((x[,1]-(x.min-360)+res)/res)])
			if(length(test.x[test.x==0])){
				foo = which(test.x==0)
				for(i in foo){
					if(i==1) x[i, ] = x[which(test.x==1)[1],]
					else if(i>1){
						x[i,] = x[i-1,]	# crude but should only affect convergence time
						}
					}
				}
			## make sure no y.inits are on land
			test.y = diag(water[trunc((y.init[,2]-y.min+res)/res),
				trunc((y.init[,1]-(x.min-360)+res)/res)])
			if(length(test.y[test.y==0])){
				foo = which(test.y==0)
				for(i in foo){
					if(i==1) y.init[i, ] = y.init[which(test.y==1)[1],]
					else if(i>1){
						y.init[i,] = y.init[i-1,] # crude but should only affect convergence time
						}
					}
				}
	    	}

        start.date = input$first.date
        tstep.sec = input$tstep * 86400
        steplims = seq(start.date, by = paste(tstep.sec, "sec"), length = length(idx))
     	jags.data = list(y = y, itau2 = input$itau2, nu = input$nu, j = j, idx = idx,
     		RegN = RegN)
        iSigma = matrix(c(1, 0, 0, 1), 2, 2)
    	if(LM){
			xy.rng = matrix(c(x.min-360,x.max-360,y.min,y.max),2,2)
			jags.data = list(y = y, itau2 = input$itau2, nu = input$nu, j = j, idx = idx,
				RegN = RegN, xy.rng = xy.rng, res = res, buffer = buffer, water = water)
			}
      if (HO) {jags.data$ho=input$ho}

			jags.inits <- list(list(iSigma=iSigma, gamma = c(0.8, NA, NA),
            	dev = c(0.6, 0.2), tmp = c(0.45, 0.55, 0.5),
            	x = x, logpsi = runif(1,-1,1), y = y.init),
            	list(iSigma=iSigma, gamma = c(0.85, NA, NA),
            	dev = c(0.65,0.45), tmp = c(0.5, 0.5, 0.5),    #c(0.6, 0.4, 0.35)
            	x = x, logpsi = runif(1,-1,1), y = y.init),
            	list(iSigma=iSigma, gamma = c(0.75, NA, NA),
            	dev = c(0.7,0.3), tmp = c(0.5, 0.6, 0.7),
            	x = x, logpsi = runif(1,-1,1), y = y.init),
            	list(iSigma=iSigma, gamma = c(0.6, NA, NA),
            	dev = c(0.3, 0.4), tmp = c(0.6, 0.5, 0.3),
            	x = x, logpsi = runif(1,-1,1), y = y.init))
            	jags.params <- c("Sigma", "x", "theta", "gamma", "phi", "b", "h",
            	"psi", "zhat")

              model.file = "jags/haulout4b.txt"


	if(chains==1) jags.inits = jags.inits[[1]]
	if(chains==2) jags.inits = list(jags.inits[[1]],jags.inits[[2]])
	if(chains==3) jags.inits = list(jags.inits[[1]],jags.inits[[2]],jags.inits[[3]])

#	model.file = paste(system.file('jags',package='bsam'), "/", model, ".txt", sep="")

	burn = jags.model(model.file, jags.data, jags.inits, n.chains=chains, n.adapt=adapt/2)
	update(burn, n.iter=adapt/2)
	psamples = jags.samples(burn, jags.params, n.iter=samples, thin=thin)

	lon = apply(psamples$x[,1,,],1, mean)
	lat = apply(psamples$x[,2,,],1, mean)
	lon.q = apply(psamples$x[,1,,],1, quantile, c(0.025, 0.5, 0.975))
	lat.q = apply(psamples$x[,2,,],1, quantile, c(0.025, 0.5, 0.975))

	         b = apply(psamples$b, 1, mean)
	         b.5 = apply(psamples$b, 1, median)
	         h = apply(psamples$h,1,mean)
	         summary = data.frame(id=as.character(id), date =
	                                as.POSIXct(as.numeric(steplims), origin="1970-01-01 00:00:00", tz="GMT"),
	                              lon, lat, lon.025=lon.q[1,], lon.5=lon.q[2,], lon.975=lon.q[3,],
	                              lat.025=lat.q[1,], lat.5=lat.q[2,], lat.975=lat.q[3,], b, b.5, h)

        step = with(summary, difftime(date[2], date[1],
                                      units="hours"))
        y = data.frame(lon=y[,1], lat=y[,2])

	out = list(summary=summary, timestep=step, data=y, N=RegN,
		mcmc=psamples)
	out
    }
    lapply(loc.list, ssm1)
}