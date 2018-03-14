setClass(
  "SDPMig",
  slots = c(Name    = "character",
            Init    = "list",
            Species = "list",
            Sites   = "list", 
            Results = "list")
)

CalcTR <- function(x, t) {
  TR <- approx(parms$xFTReward, parms$yFTReward, t, rule = 2)$y
  s  <- parms$w * (x - parms$xc) 
  
  if(x == 0) { SR <- 0 } else {
    SR <- (((exp(s) - exp(-s)) / (exp(s) + exp(-s))) + 1)/2
  }
  TR * SR + parms$B0
}


makeSDPmig <- function(parms, project = "") {
  
  sparms <- read.csv(parms$path, skip = 2, sep = ";", dec = ".")
  r_gain.x <- as.matrix(sparms[,which(substring(names(sparms), 1,1)=="x")])
  r_gain.y <- as.matrix(sparms[,which(substring(names(sparms), 1,1)=="y")])
  r_gain.p <- matrix(sparms$stoch, ncol = ncol(r_gain.x), nrow = nrow(r_gain.y))
  r_dist   <- geosphere::distm(cbind(sparms$Lon, sparms$Lat), cbind(sparms$Lon, sparms$Lat))/1000
  
  new(
    "SDPMig",
    Name = project,
    Init = list(
      MaxT   = parms$MaxT,
      NSites = parms$NSites,
      MaxX   = parms$MaxX
    ),
    Species = list(
      B0    = parms$B0,
      w     = parms$w,
      xc    = parms$xc,
      c     = parms$c,
      speed = parms$speed,
      WindAssist = parms$WindAssist,
      WindProb   = parms$WindProb,
      ZStdNorm   = parms$ZStdNorm,
      PStdNorm   = parms$PStdNorm,
      xFTReward  = parms$xFTReward,
      yFTReward  = parms$yFTReward,
      decError   = parms$decError
    ),
    Sites = list(
      crds = cbind(sparms$Lon, sparms$Lat),
      dist = r_dist,
      b0   = sparms$b0,
      b1   = sparms$b1,
      b2   = sparms$b2,
      expend  =  sparms$expend,
      pred_a1 =  parms$pred_a1,
      pred_a2 =  parms$pred_a2,
      gain    = list(gain.x = r_gain.x, 
                     gain.y = r_gain.y, 
                     gain.p = r_gain.p)
    ),
    Results = list(
      FitnessMatrix     = NA,
      DecisionMatrix    = NA,
      ProbMatrix        = NA
    )
  )

}


### Backward iteration

bwdIteration <- function(obj) {
  
  Init(obj@Init$MaxT, obj@Init$NSites, obj@Init$MaxX,
       obj@Species$w,obj@Species$xc,obj@Species$B0,obj@Sites$b0,obj@Sites$b1,obj@Sites$b2,obj@Sites$pred_a1,obj@Sites$pred_a2,
       obj@Species$c,obj@Species$speed,obj@Species$WindAssist,obj@Species$WindProb,
       obj@Species$ZStdNorm,obj@Species$PStdNorm,
       obj@Species$xFTReward,obj@Species$yFTReward,
       obj@Species$decError,
       obj@Sites$dist,
       obj@Sites$gain[[1]],obj@Sites$gain[[2]],obj@Sites$gain[[3]],obj@Sites$expend)
  
  out <- BackwardIteration()
  
  obj@Results$FitnessMatrix <- out[[1]]
  DM <- array(dim = c(dim(out[[2]]),2))
  DM[,,,1] <- out[[2]] 
  DM[,,,2] <- out[[3]]
  obj@Results$DecisionMatrix <- DM  
  PM <- array(dim = c(dim(out[[4]]),2))
  PM[,,,1] <- out[[4]] 
  PM[,,,2] <- out[[5]] 
  obj@Results$ProbMatrix <- PM
  
  obj 
}



############################
#### Forward Simulation#####
############################

MigSim <- function(obj, NrInd, start_t, start_site, start_x) {
  
  InitSim(obj@Init$MaxT, obj@Init$NSites, obj@Init$MaxX,
          obj@Species$w,obj@Species$xc,obj@Species$B0,obj@Sites$b0,obj@Sites$b1,obj@Sites$b2,obj@Sites$pred_a1,obj@Sites$pred_a2,
          obj@Species$c,obj@Species$speed,obj@Species$WindAssist,obj@Species$WindProb,
          obj@Species$ZStdNorm,obj@Species$PStdNorm,
          obj@Species$xFTReward,obj@Species$yFTReward,
          obj@Species$decError,
          obj@Sites$dist,
          obj@Sites$gain[[1]],obj@Sites$gain[[2]],obj@Sites$gain[[3]],obj@Sites$expend,
          obj@Results$FitnessMatrix,
          obj@Results$DecisionMatrix[,,,1], obj@Results$DecisionMatrix[,,,2],
          obj@Results$ProbMatrix[,,,1], obj@Results$ProbMatrix[,,,2])

  forwardSim(NrInd, start_t-1, start_site-1, start_x)
  
}


############################
### Plotting ###############
############################


fitnessPlot <- function(x, cond = c(20, 80)) {
  
  opar <- par(oma = c(0,1,1,4), mar = c(4,5,2,1))
  opar <- par(fig=c(0,0.5,0.5,1))
  plot(NA, xlim = c(0.5, x@Init$NSites+0.5), ylim = c(1, x@Init$MaxT), yaxt = "n", ylab = "time", xlab = "", xaxs = "i", yaxs = "i")
  image(x = 1:(x@Init$NSites+1), y = 0:x@Init$MaxT,
        t(x@Results$FitnessMatrix[,,cond[1]])[,nrow(x@Results$FitnessMatrix[,,cond[1]]):1],
        breaks = seq(x@Species$B0-0.1, max(x@Results$FitnessMatrix[,,x@Init$MaxX+1]+0.1, na.rm = T), length = 100),
        col = rev(heat.colors(99)), add = T)
  axis(2, at = seq(1, x@Init$MaxT, length = 10), labels = seq(x@Init$MaxT, 1, length = 10), las = 1)
  mtext(paste0("Fitness: x = ", cond[1]), 3, line = 1)
  
  opar <- par(fig=c(0.5,1,0.5,1), new = T)
  plot(NA, xlim = c(0.5, x@Init$NSites+0.5), ylim = c(1, x@Init$MaxT), yaxt = "n", ylab = "", xlab = "", xaxs = "i", yaxs = "i")
  image(x = 1:(x@Init$NSites+1), y = 0:x@Init$MaxT,
        t(x@Results$FitnessMatrix[,,cond[2]])[,nrow(x@Results$FitnessMatrix[,,cond[2]]):1],
        breaks = seq(x@Species$B0-0.1, max(x@Results$FitnessMatrix[,,x@Init$MaxX+1]+0.1, na.rm = T), length = 100),
        col = rev(heat.colors(99)), add = T)
  mtext(paste0("Fitness: x = ", cond[2]), 3, line = 1)
  axis(2, at = seq(1, x@Init$MaxT, length = 10), labels = seq(x@Init$MaxT, 1, length = 10), las = 1)
  image.plot(t(x@Results$FitnessMatrix[,,cond[2]])[,nrow(x@Results$FitnessMatrix[,,cond[2]]):1],
             breaks = seq(x@Species$B0-0.1, max(x@Results$FitnessMatrix[,,x@Init$MaxX+1]+0.1, na.rm = T), length = 100),
             col = rev(heat.colors(99)), legend.only = T, legend.mar = 0)
  
  opar <- par(fig=c(0,1,0,0.5), new = T, mar = c(5, 5, 4, 1))
  plot(NA, xlim = c(1, x@Init$MaxT), ylim = c(0, max(x@Sites$gain[[2]], na.rm = T)), xlab = "t", ylab = "gain")
  matplot(x = t(x@Sites$gain[[1]]), y = t(x@Sites$gain[[2]]), type = "l", lwd = 1.4, lty = 1, col = rainbow(nrow(x@Sites$gain[[2]])),
          add = T, xlab = "", ylab = "")
  legend(1, x@Init$NSites+1, paste0("site ", c(1, x@Init$NSites-1)), pch = 16, col = rainbow(nrow(x@Sites$gain[[2]]))[c(1,x@Init$NSites-1)],
         horiz = T, bty = "n", xpd = T, cex = 0.8)
  par(opar)
}


fitnessSitePlot <- function(x, mfrow = c(3,3)) {
  
  opar <- par(mfrow = mfrow, mar = c(0,0,0,0), oma = c(4,4,1,4), new = F)
  for(i in 1:(x@Init$NSites+1)) {
    matplot(1:dim(x@Results$FitnessMatrix)[3], x@Results$FitnessMatrix[,i,-1], pch = 16, type = "p", cex = 1.25, 
            col = viridis::inferno(dim(x@Results$FitnessMatrix)[3], alpha = 0.8), ylim = c(0, max(x@Results$FitnessMatrix, na.rm = T)),
            xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    
    if(i==1) axis(2, las = 1)
    if(i==(mfrow^2)[1]) axis(1)
    
    par(new = T)
    if(i!=(x@Init$NSites+1)) {
      plot(x@Sites$gain[[1]][i,], x@Sites$gain[[2]][i,], type = "l", lty = 2, lwd = 2, col = "darkgreen",
           xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0, max(x@Sites$gain[[2]], na.rm = T)))
    }
    
    if(i==(mfrow)[1]) axis(4, las = 1)
  }
  
  mtext("Time",   1, outer= T, line = 2)
  mtext("Reward", 2, outer= T, line = 2)  
  mtext("Gain",   4, outer= T, line = 2)  
  par(opar)
  
}


decisionPlot <- function(x, time = c(5, 40, 80)) {
  
  opar <- par(fig=c(0,0.33,0.5,1), mar = c(.5,.5,.5,.5), oma = c(3,3,2,2), new = F)
  plot(NA, xlim = c(0, x@Init$MaxX), ylim = c(1, x@Init$NSites), xaxt = "n", yaxt = "n", xlab = "", ylab = "", xaxs = "i")
  
  tmp2 <- abs(t(x@Results$DecisionMatrix[,time[1],,1])+1)
  
  image(x = 0:x@Init$MaxX, y = 1:(x@Init$NSites+1), tmp2, 
        add = T, col = rev(heat.colors(99)), breaks = seq(0, 1, length = 100))
  
  ind  <- runif(length(as.vector(tmp2))) < t(x@Results$ProbMatrix[,time[1],,1])
  tmp1 <- t(x@Results$ProbMatrix[,time[1],,1])
  tmp1[ind] <- NA
  
  image(x = 0:x@Init$MaxX, y = 1:(x@Init$NSites+1), tmp1, 
       add = T, col = grey(0:8/8), breaks = seq(0, 1, length = 10))
  
  axis(1, labels = F)
  axis(2, at = 1:x@Init$NSites, labels = T, las = 1)
  mtext(paste("t =", time[1]), 3, line = 0.5)
  
  opar <- par(fig=c(0.33,0.66,0.5,1), mar = c(.5,.5,.5,.5), oma = c(3,3,2,2), new = T)
  plot(NA, xlim = c(0, x@Init$MaxX), ylim = c(1, x@Init$NSites), xaxt = "n", yaxt = "n", xlab = "", ylab = "", xaxs = "i")
  
  tmp2 <- abs(t(x@Results$DecisionMatrix[,time[2],,1])+1)
  
  image(x = 0:x@Init$MaxX, y = 1:(x@Init$NSites+1), tmp2, 
       add = T, col = rev(heat.colors(99)), breaks = seq(0, 1, length = 100))
  
  ind  <- runif(length(as.vector(tmp1))) < t(x@Results$ProbMatrix[,time[2],,1])
  tmp1 <- t(x@Results$ProbMatrix[,time[2],,1])
  tmp1[ind] <- NA
   
  image(x = 0:x@Init$MaxX, y = 1:(x@Init$NSites+1), tmp1, 
       add = T, col = grey(0:8/8), breaks = seq(0, 1, length = 10))
   
  
  axis(1, labels = F)
  axis(2, at = 1:x@Init$NSites, labels = F)
  mtext(paste("t =", time[2]), 3, line = 0.5)
  
  opar <- par(fig=c(0.66,1,0.5,1), mar = c(.5,.5,.5,.5), oma = c(3,3,2,2), new = T)
  plot(NA, xlim = c(0, x@Init$MaxX), ylim = c(1, x@Init$NSites), xaxt = "n", yaxt = "n", xlab = "", ylab = "", xaxs = "i")
  
  tmp2 <- abs(t(x@Results$DecisionMatrix[,time[3],,1])+1)
  
  image(x = 0:x@Init$MaxX, y = 1:(x@Init$NSites+1), tmp2, 
        add = T, col = rev(heat.colors(99)), breaks = seq(0, 1, length = 100))
  
  ind  <- runif(length(as.vector(tmp1))) < t(x@Results$ProbMatrix[,time[3],,1])
  tmp1 <- t(x@Results$ProbMatrix[,time[3],,1])
  tmp1[ind] <- NA
  
  image(x = 0:x@Init$MaxX, y = 1:(x@Init$NSites+1), tmp1, 
       add = T, col = grey(0:8/8), breaks = seq(0, 1, length = 10))
   
  axis(1, labels = F)
  axis(2, at = 1:x@Init$NSites, labels = F)
  mtext(paste("t =", time[3]), 3, line = 0.5)
  
  opar <- par(fig=c(0,0.33,0,0.5), new = T)
  plot(NA, xlim = c(0, x@Init$MaxX), ylim = c(1, x@Init$NSites+1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", xaxs = "i")
  image(x = 0:x@Init$MaxX, y = 1:(x@Init$NSites+1), t(x@Results$DecisionMatrix[,time[1],,2]+1), 
        add = T, col = rainbow(x@Init$NSites), breaks = c(0:x@Init$NSites+1))
  axis(1, labels = T)
  axis(2, at = 1:x@Init$NSites)
  
  opar <- par(fig=c(0.33,0.66,0,0.5), new = T)
  plot(NA, xlim = c(0, x@Init$MaxX), ylim = c(1, x@Init$NSites+1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", xaxs = "i")
  image(x = 0:x@Init$MaxX, y = 1:(x@Init$NSites+1), t(x@Results$DecisionMatrix[,time[2],,2])+1, 
        add = T, col = rainbow(x@Init$NSites+1), breaks = c(0:(x@Init$NSites+1)))
  axis(1, labels = F)
  axis(2, at = 1:x@Init$NSites, labels = F)
  
  opar <- par(fig=c(0.66,1,0,0.5), new = T)
  plot(NA, xlim = c(0, x@Init$MaxX), ylim = c(1, x@Init$NSites+1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", xaxs = "i")
  image(x = 0:x@Init$MaxX, y = 1:(x@Init$NSites+1), t(x@Results$DecisionMatrix[,time[3],,2])+1, 
        add = T, col = rainbow(x@Init$NSites+1), breaks = c(0:(x@Init$NSites+1)))
  axis(1, labels = T)
  axis(2, at = 1:x@Init$NSites, labels = F)
  
  mtext("site", 2, outer = T, line = 1.5, cex = 1.2)
  mtext("condition(x)", 1, outer = T, line = 1.5, cex = 1.2)
  
  legend("topright", as.character(1:x@Init$NSites), pch = 16, col = rainbow(x@Init$NSites))
  
  
  par(opar)
  
}


simuPlot <- function(simu, sdpM, fun = "median") {
  
  opar <- par(mfrow = c(3,1), mar = c(4,4,1,1))
  plot(NA, xlim = c(0, sdpM@Init$MaxT+1), ylim = c(1, sdpM@Init$NSites+1), bty = "n", xlab = "time", ylab = "site", 
       yaxt = "n", las = 1, cex.lab  = 1.2)
  axis(2, at = 1:sdpM@Init$NSites)
  for(i in 1:dim(simu)[1]) {
    tmp <- t(simu[i,,])
    if(sum(tmp[,5])>0) tmp <- tmp[1:(min(which(tmp[,5]==1))),]
    lines(tmp[,1], tmp[,2]+1, lwd = 2, col = rgb(0.4,0.6,0.9, alpha = 0.5))
    if(sum(tmp[,5])>0) {
      points(tmp[nrow(tmp),1], tmp[nrow(tmp),2]+1, pch = "x") 
    }  else {
      points(tmp[nrow(tmp),1], tmp[nrow(tmp),2]+1, col = rgb(0.4,0.6,0.9, alpha = 0.5), pch = 16)
    }
  }
  legend("topleft", "death", pch = "x", bty = "n")
  
  
  plot(NA, xlim = c(0, sdpM@Init$MaxT+1), ylim = c(1, sdpM@Init$MaxX), bty = "n", xlab = "time", ylab = "condition", 
       las = 1, cex.lab  = 1.2, bty = "o")
  
  for(i in 1:dim(simu)[1]) {
    tmp <- simu[i,,]
    tmp <- tmp[,tmp[5,]!=1]
    if(!is.null(nrow(tmp))) lines(tmp[1,], tmp[3,], col = adjustcolor("grey60", alpha.f = 0.6))
  }
  for(i in 1:dim(simu)[1]) {
    tmp <- simu[i,,]
    tmp <- tmp[,!is.na(tmp[2,]) & tmp[5,]!=1]
    if(!is.null(nrow(tmp))) points(tmp[1,], tmp[3,], col = c(rainbow(sdpM@Init$NSites+1, alpha = 0.6)[tmp[2,]+1]), pch = 16) else {
      points(tmp[1], tmp[3], col = c(rainbow(sdpM@Init$NSites, alpha = 0.6)[tmp[2]]), pch = 16) 
    }
  }
  
  legend("topleft", paste("site", 1:sdpM@Init$NSites), pch = 16, col = rainbow(sdpM@Init$NSites))
  
  res <- matrix(NA, ncol = sdpM@Init$NSites+1, nrow = dim(simu)[1])
    
  for(i in 1:nrow(res)) {
      tmp <- simu[i,,]
      if(sum(tmp[5,])>0) tmp <- tmp[,1:(min(which(tmp[5,]==1)))]
      res[i,] <- apply(cbind(0:sdpM@Init$NSites), 1, function(x) {
        ifelse(all(tmp[2,]<x) , NA, sum(tmp[2,]==x, na.rm = T))
      })
  }

  plot(NA, xlim = c(1, sdpM@Init$NSites), ylim = c(0, max(res, na.rm = T)), xlab = "time", ylab = "days")
  matplot(t(res[,-ncol(res)]), type = "b", pch = 16, cex = 1.5, lty = 1, lwd = 1, col = adjustcolor("grey60", alpha.f = 0.7), add = T)
  
  points(1:(sdpM@Init$NSites), apply(res[,-ncol(res)], 2, fun, na.rm = T), type = "b", cex = 3, pch = 16, col = adjustcolor("firebrick", alpha.f = 0.9))
  
  par(opar)
  
  res
}

