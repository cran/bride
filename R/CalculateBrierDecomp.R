#the main function
CalculateBrierDecomp <- function(p, y, n.bins) {

  #initialize main object to return
  bride.list <- list()
  bride.list$p <- p
  bride.list$y <- y

  #p bins
  p.breaks <- seq(0, 1, length.out = n.bins + 1) + 
               c(-0.01, rep(0, n.bins - 1), 0.01)
  ##get bin no. of each p-value
  p.binning <- cut(p, breaks=p.breaks, include.lowest=TRUE, ordered_result=TRUE)
  p.binning <- as.numeric(p.binning)
  # calculate the in-bin averages
  p.means <- lapply( seq(n.bins),
               function(z){ 
                 mean(p[which(p.binning == z)]) 
               }
             ) 
  p.means <- unlist(p.means)
  # replace p-values by in-bin averages
  p <- p.means[as.vector(unlist(p.binning))]

  #store
  bride.list$n.bins <- n.bins
  bride.list$p.binavgs <- p.means
  bride.list$p.binned <- p


  ##calculate cross table
  # p-bins
  c.names <- paste(seq(n.bins))
  # y-classes
  r.names <- paste(sort(unique(y)))
  # expand grid for each entry of cross table
  ctab.grid <- as.matrix(expand.grid(r.names, c.names))
  # calculate individual elements of cross table
  ctab.grid <- apply(ctab.grid, 1, 
                 function(z) { 
                   length( intersect( which( p.binning == z[2] ), 
                                      which( y == z[1] ) )) 
                 } 
               )
  c.table <- matrix(ctab.grid, ncol=n.bins)
  colnames(c.table) <- c.names
  rownames(c.table) <- r.names

  #store
  bride.list$c.table <- c.table

  n.od <- colSums(c.table)
  n <- sum(c.table)


  #biased decomposition
  inds <- which(n.od > 0)
  rel <- sum((n.od / n * (p.means - c.table["1", ] / n.od)^2)[inds])
  res <- sum((n.od / n * 
              (c.table["1", ] / n.od - sum(c.table["1", ]) / n)^2)[inds] )
  unc <- sum(c.table["0", inds]) * sum(c.table["1", inds]) / n^2

  #store
  bride.list$rel <- rel
  bride.list$res <- res
  bride.list$unc <- unc
  bride.list$br <- rel - res + unc

  #bias-corrected decomposition
  inds <- which(n.od > 1) #avoid division by zero
  r1 <- -1 * sum(((c.table["0", ] * c.table["1", ]) / 
        (n.od * (n.od - 1)))[inds] ) / n
  r2 <- sum(c.table["0", inds]) * sum(c.table["1", inds]) / 
        (n^2 * (n - 1))
  rel2a <- rel + r1
  res2a <- res + r1 + r2
  unc2 <- unc + r2

  rel2 <- max(c(rel2a, rel2a - res2a, 0))
  res2 <- max(c(res2a, res2a - rel2a, 0))

  #store
  bride.list$rel2 <- rel2
  bride.list$res2 <- res2
  bride.list$unc2 <- unc2

  #variances of reliability, resolution, uncertainty
  vrel <- 0.0
  vres <- 0.0
  vrelunb <- 0.0
  vresunb <- 0.0
  I4 <- y
  v4 <- var(I4)
  Fd <- sum(I4)
  dFresunb <- (n - 2 * Fd) / n^3 / (n - 1)

  for (d in seq(n.bins)) {
    I1 <- 1 * (p.binning == d)
    if (sum(I1) > 1) {
      I2 <- I1 * y
      I3 <- I1 * p
      v1 <- var(I1)
      v2 <- var(I2)
      v3 <- var(I3)
      c12 <- cov(I1, I2)
      c23 <- cov(I2, I3)
      c13 <- cov(I1, I3)
      c14 <- cov(I1, I4)
      c24 <- cov(I2, I4)
  
      Ad <- sum(I1)
      Bd <- sum(I2)
      Cd <- sum(I3)
  
      f <- Bd - Cd
      vrel <- vrel + (f / Ad)^4 * v1 + (2.0 * f / Ad)^2 * (v2 + v3) - 
                       4.0 * (f / Ad)^3 * (c12 - c13) - 
                       8.0 * (f / Ad)^2 * c23
  
      f1 <- Bd / Ad - Fd / n
      f2 <- Bd / Ad + Fd / n
      vres <- vres + f1 * f1 * f2 * f2 * v1 + 4.0 * f1 * f1 * v2 - 
                       4.0 * f1 * f1 * f2 * c12 
  
      dArelunb <- -1.0 / (n * Ad^2) * 
                  ((Bd - Cd)^2 - Ad * Bd / (Ad - 1) - 
                   Bd * (Bd - Ad) / ((Ad - 1)^2))
      dBrelunb <- (2 * Bd - 1) / (Ad - 1) / n - 2 * Cd / n / Ad
      dCrelunb <- -1. / (n * Ad) * 2 * (Bd - Cd)
      vrelunb <- vrelunb + n * (dArelunb^2 * v1 + 
                                dBrelunb^2 * v2 +
                                dCrelunb^2 * v3 +
                                2 * dArelunb * dBrelunb * c12 +
                                2 * dArelunb * dCrelunb * c13 +
                                2 * dBrelunb * dCrelunb * c23 )
      dAresunb <- -f1 * f2 / n +
                   Bd / n / Ad^2 / (Ad-1)^2 * ((Ad-Bd)^2 - Bd * (Bd-1))
      dBresunb <- 2 * f1 / n - (Ad - 2 * Bd) / Ad / n / (Ad - 1)
      vresunb <- vresunb + n * (dAresunb^2 * v1 +
                                  dBresunb^2 * v2 +
                                  2 * dAresunb * dBresunb * c12 +
                                  2 *dAresunb * dFresunb * c14 +
                                  2 * dBresunb * dFresunb * c24)
    }

  }
  vresunb <- vresunb + n * (dFresunb^2 * v4)
  vunc <- (1 - 2 * Fd / n)^2 / n * v4
  vuncunb <- (n - 2 * Fd)^2 / n / (n - 1)^2 * v4
  vrel <- vrel / n
  vres <- vres / n

  #store
  bride.list$relv <- vrel
  bride.list$resv <- vres
  bride.list$uncv <- vunc
  bride.list$rel2v <- vrelunb
  bride.list$res2v <- vresunb
  bride.list$unc2v <- vuncunb

  #return
  bride.list 
}

