library(SDPmig) ## R package available on supplemental repository

###############################################
### Parameters ################################
###############################################

parms <- list(
  MaxT   = 100,  ## Maximum time
  MaxX   = 100,  ## Maximum body condition
  NSites = 8,    ## Excl. Breeding
  
  ### Species specific Parameters ##
  
  B0 = 3,       ## Future reproductive success
  w  = 0.028,   ## Parameters for sigmoidal TR function
  xc = 55,      ## 
  
  ## Flying capacities
  c     = 14776,
  speed = 1440,
  
  ## Wind
  WindAssist = 0,
  WindProb   = 1,
  
  ## Decision Error
  decError = 2000,
  
  ## Terminal Reward
  xFTReward = c(0, 86, 87, 97, 98, 100),
  yFTReward = c(0, 0,   2,  2,  0,   0),
  
  
  ### Site specific Parameters ###
  path = "Paramters/WFGparameter.csv",
  
  
  pred_a1  = 2,
  pred_a2  = 2,
  
  ### Accuracy
  ZStdNorm = c(-2.5, -2.0, -1.5, -1.0, -0.5,  0.0,  0.5,  1.0,  1.5,  2.0,  2.5),
  PStdNorm = c(0.0092, 0.0279, 0.0655, 0.1210, 0.1747, 0.2034, 0.1747, 0.1210, 0.0655, 0.0279, 0.0092)
  
) ## End parameter ####


#######################
### Simulation ########
#######################

dist <- as.matrix(read.table("Parameters/dist.txt"))


sdp0  <- makeSDPmig(parms, "GWF-Geese: Hunting")
sdp0@Sites$dist  <- as.matrix(read.table("Parameters/dist.txt"))



### Simulation with different hunting intensities at different sites

ind <- cbind(0:2, rep(c(1:8, 158), each = 3)) # 0:2 = hunting intensities (default, +, ++) , at specific sites column 2
res <- matrix(ncol = 11, nrow = nrow(ind))

for(s in 1:nrow(ind)) {
  
  sdp@Sites$b0 <- rep(1e-4, 9)
  sdp@Sites$b1 <- rep(1e-3, 9)
  sdp@Sites$b2 <- rep(1e-4, 9)
  
  
  sites <- as.numeric(strsplit(as.character(ind[s,2]), "")[[1]])
  sdp@Sites$b0[sites] <- ifelse(ind[s,1]==0, 1e-4, ifelse(ind[s,1]==1, 1e-3, 1e-2))
  sdp@Sites$b1[sites] <- ifelse(ind[s,1]==0, 1e-3, ifelse(ind[s,1]==1, 1e-2, 1e-1))
  sdp@Sites$b2[sites] <- ifelse(ind[s,1]==0, 1e-4, ifelse(ind[s,1]==1, 1e-3, 1e-2))
  
  sdpM <- bwdIteration(sdp)
  simu <- MigSim(sdpM, 100, 1, 1, c(33, 10))
  
  tmp0 <- matrix(NA, ncol = sdpM@Init$NSites+1, nrow = dim(simu)[1])
  
  for(i in 1:nrow(tmp0)) {
    tmp <- simu[i,,]
    if(sum(tmp[5,])>0) tmp <- tmp[,1:(min(which(tmp[5,]==1)))]
    tmp0[i,] <- apply(cbind(0:sdpM@Init$NSites), 1, function(x) {
      ifelse(all(tmp[2,]<x) , NA, sum(tmp[2,]==x, na.rm = T))
    })
  }
  
  tmp2 <- apply(simu, 1, function(x) {
    if(sum(x[5,])>0) {
      tmp <- x[3,min(which(x[5,]==1))]
      ifelse(tmp>0, 1, 2) } else NA
  })
  tmp3 <- apply(cbind(c(1,2)), 1, function(x) sum(!is.na(tmp2) & tmp2==x))
  
  
  out <- c(apply(tmp0, 2, mean, na.rm = T)[-ncol(tmp0)],
           sum(apply(tmp0, 1, function(x) sum(x, na.rm = T)<100))/100, tmp3)
  res[s,] <- out
}


### Plot
layout(matrix(1:(9*2), byrow = T, ncol = 2), width = c(1, 0.3))
opar <- par(mar = c(1,6,1,0), oma = c(3,0,0,0), las = 1)

for(i in 1:9) {
  
  tmp <- res[ind[,2]==unique(ind[,2])[i],]
  bp <- barplot(tmp[-1,1:8], beside = T, ylim = c(0, 60),
                col = c(adjustcolor("indianred4", alpha.f = 0.5), adjustcolor("indianred4", alpha.f = 0.9)), ylab = "", xlab = "")
  points(bp[1,], tmp[1,1:8], pch = 16, cex = 1.4)
  points(bp[2,], tmp[1,1:8], pch = 16, cex = 1.4)
  
  arrows(bp[1,], tmp[2,1:8], bp[1,], tmp[1,1:8], length = 0, lwd = 2.5)
  arrows(bp[2,], tmp[3,1:8], bp[2,], tmp[1,1:8], length = 0, lwd = 2.5)
  
  if(i == 9) axis(1, at = apply(bp, 2, mean), labels = c("NL", "D", "PL", "L/U", "E/T", "K/K", "Ark", "Nem"))
  
  bp2 <- barplot(tmp[,9], col = c("white", adjustcolor("indianred4", alpha.f = 0.5), adjustcolor("indianred4", alpha.f = 0.9)),
                 ylim = c(0,1), xlab = "", ylab = "")
  
  if(i == 9) axis(1, at = bp2, labels = c("0", "+", "++"))
  
}
mtext("Days at Site", 2, outer = T, line = -3, las = 3, cex = 1.2)
# mtext("Mortality", 2, outer = T, line = -39, las = 3, cex = 1.2)
# mtext("Reason for Mortality (shot/starvation)", 2, outer = T, line = -50, las = 3, cex = 1.2)

par(opar)





