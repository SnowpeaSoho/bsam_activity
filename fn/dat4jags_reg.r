dat4jags   <- function (indata, tstep = 1, ...)
{
    tstep.sec <- tstep * 86400
    tstep.hr <- tstep * 24
    datetime <- as.POSIXct(indata[, 2], format = "%Y-%m-%d %H:%M:%S",tz = "GMT")
    if (ncol(indata) == 5) {
        dat <- data.frame(id = indata[, 1], time = datetime,
            lc = indata[, 3], lon = indata[, 4], lat = indata[,
                5])
    }
    else if (ncol(indata) == 7) {
        dat <- data.frame(id = indata[, 1], time = datetime,
            lc = indata[, 3], lon = indata[, 4], lat = indata[,
                5], lonerr = indata[, 6], laterr = indata[, 7])
    }
    dat[order(dat[, 1], dat[, 2]), ]
    
    # add an (empty) first position to fix the start time at regular interval
    restack <- function(x) {
    for (i in 1:length(x)) {
    	 if (i==1) {x2 <-x[[1]]}
    	 if (i>1) {x2 <- rbind(x2,x[[i]])}
    	}
      x2}

    fix.start <- function(k) {
      k <- rbind(k[1, ], k) # insert a repeat of 1st pos
      k$lc[1] <- "Z"
      date1 <- trunc(k$time[1], units = "days") + c(0:c(24/tstep.hr))*tstep.hr*60*60
      k$time[1] <- date1[which(date1>k$time[1])[1]-1]
      k}

    dat2 <- by(dat, dat$id, fix.start)
    dat <- restack(dat2)

    dat$lc <- factor(dat$lc, levels = c("3", "2", "1", "0", "A",
        "B", "Z", "F", "G"), ordered = TRUE)
    sigma.lon <- c(0.289866, 0.3119293, 0.9020423, 2.1625936,
        0.507292, 4.2050261, 4.2050261, 0.01)
    sigma.lat <- c(0.1220553, 0.2605126, 0.4603374, 1.607056,
        0.5105468, 3.041276, 3.041276, 0.01)
    nu.lon <- c(3.070609, 1.220822, 2.298819, 0.9136517, 0.786954,
        1.079216, 1.079216, 1e+05)
    nu.lat <- c(2.075642, 6.314726, 3.896554, 1.010729, 1.057779,
        1.331283, 1.331283, 1e+05)
    nu.lon <- nu.lon[as.numeric(dat$lc)]
    nu.lat <- nu.lat[as.numeric(dat$lc)]
    sigma.lon <- (sigma.lon/6366.71 * 180)/pi
    sigma.lat <- (sigma.lat/6366.71 * 180)/pi
    itau2.lon <- sigma.lon[as.numeric(dat$lc)]^-2
    itau2.lat <- sigma.lat[as.numeric(dat$lc)]^-2
    if (ncol(dat) == 7) {
        itau2.lon[which(dat$lc == "G")] = dat$lonerr[which(dat$lc ==
            "G")]^-2
        itau2.lat[which(dat$lc == "G")] = dat$laterr[which(dat$lc ==
            "G")]^-2
        nu.lat[which(dat$lc == "G")] = nu.lon[which(dat$lc ==
            "G")] = 1e+05
    }
    dat$itau2.lon <- itau2.lon
    dat$itau2.lat <- itau2.lat
    dat$nu.lon <- nu.lon
    dat$nu.lat <- nu.lat
    dostep <- function(k) {
        tt <- k$time
        tst.nona <- cut(tt, paste(tstep.sec, "sec"), labels = FALSE,
            incl = TRUE)
        tst.seq <- seq(1, max(tst.nona))
        tst.rle <- rle(tst.nona)
        tst.val <- tst.rle$values
        tst.isna <- !tst.seq %in% tst.val
        k.idx <- tst.all <- sort(c(tst.nona, tst.seq[tst.isna]))
        tstall.rle <- rle(tst.all)
        idx <- cumsum(c(1, tstall.rle$lengths))
        steplims <- seq(tt[1], by = paste(tstep.sec, "sec"),
            length = length(idx))
        tstall.isna <- tst.all %in% tst.seq[tst.isna]
        k.idx[tstall.isna] <- NA
        k.idx[!tstall.isna] <- seq(nrow(k))
        k.new <- k[k.idx, ]
        step.frac <- as.numeric(difftime(k.new$time, steplims[tst.all],
            units = "sec"))/tstep.sec
        step.frac[is.na(step.frac)] <- 0.5
        itau2lon.isna <- is.na(k.new$itau2.lon)
        itau2lat.isna <- is.na(k.new$itau2.lat)
        k.new$itau2.lon[itau2lon.isna] <- min(k.new$itau2.lon,
            na.rm = TRUE)
        k.new$itau2.lat[itau2lat.isna] <- min(k.new$itau2.lat,
            na.rm = TRUE)
        k.new$nu.lon[k.new$nu.lon < 2 | is.na(k.new$nu.lon)] <- 2
        k.new$nu.lat[k.new$nu.lat < 2 | is.na(k.new$nu.lat)] <- 2
        list(id = k.new$id[1], y = cbind(k.new$lon, k.new$lat),
            itau2 = cbind(k.new$itau2.lon, k.new$itau2.lat),
            nu = cbind(k.new$nu.lon, k.new$nu.lat), idx = idx,
            j = step.frac, RegN = length(idx), first.date = k.new$time[1],
            dts = k.new$time, tstep = tstep)
    }
    by(dat, dat$id, dostep)
}
