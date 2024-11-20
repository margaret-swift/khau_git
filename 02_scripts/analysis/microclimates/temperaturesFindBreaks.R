findBestBreaks <- function(sp, ws) {
  aic.df <- data.frame(expand.grid(ws[[1]], ws[[2]]),
                       AIC=NA)
  for (i in 1:nrow(aic.df)) {
    breaks <- as.numeric(aic.df[i, c(1,2)])
    m <- modelMe(sp, breaks)
    aic.df$AIC[i] <- AIC(m)
  }
  inx <- which.min(aic.df$AIC)
  best <- aic.df[inx,]
  
  # handling edge cases
  if (which(ws[[1]] == best[,1]) %in% c(1, length(ws[[1]]))) {
    message('expand window one!')
  }
  if (which(ws[[2]] == best[,2]) %in% c(1, length(ws[[2]]))) {
    message('expand window two!')
  }
  return(best)
}



# windows <- list(seq(10, 20, by=0.1),
# seq(25, 35, by=0.1))
# findBestBreaks("Oryx", windows)

# windows <- list(seq(12, 16, by=0.5), 
#                 seq(30, 32, by=0.5))
# findBestBreaks("Roan", windows)