### Credible Intervals

### Sets Up Libraries
devtools::load_all(".")
library(Matrix)
library(ggplot2)

iter <- 30
### Desired Set-up to Evaluate
p <- c(50, 50, 50)          # Matrix size for each dimension
d <- length(p)              # Number of tensor dimensions
r <- c(0.10, 0.20, 0.30)    # Sparity level for Random Matrices Only
b <- 103                    # Degrees of Freedom for Precision Matrices
n <- 10                     # Number of Observations
burnin <- 20
nmcmc <-1000                 # Number of MCMC samples

### Obtains the Coverage
preCov <- matrix(NA,
                 nrow = iter,
                 ncol = d + 1)
preLen <- matrix(NA,
                 nrow = iter,
                 ncol = d + 1)
for(ii in 1:iter){
  ### Reads the Matrices
  ### General Suffix
  suffix <- paste0(p[1], "-", p[2], "-", p[3], "-",
                   r[1], "-", r[2], "-", r[3], "-",
                   n, "-", nmcmc, "-", burnin, "-", ii)

  ### Loads R Objects
  filNam <- paste0("./out/BFit-", suffix, ".rds")
  BFit   <- readRDS(file = filNam)
  filNam <- paste0("./out/O-", suffix, ".rds")
  O      <- readRDS(file = filNam)
  filNam <- paste0("./out/E-", suffix, ".rds")
  E      <- readRDS(file = filNam)

  ### Changes the Precision matrices to the Correlation
  Ohat <- list()
  for(i in 1:d){
    Ohat[[i]] <- array(data = NA,
                       dim  = dim(BFit$samC[[i]]))
    for(j in 1:nmcmc){
      Ohat[[i]][,,j] <- cov2cor(BFit$samC[[i]][,,j])
    }
  }

  totSel <- 0
  totCov <- 0
  totLen <- 0
  for(i in 1:d){
    ### Computes the Quantiles
    uppQua <- apply(X      = Ohat[[i]],
                    MARGIN = c(1,2),
                    FUN    = quantile,
                    probs  = 0.975)
    lowQua <- apply(X      = Ohat[[i]],
                    MARGIN = c(1,2),
                    FUN    = quantile,
                    probs  = 0.025)
    ### Selects the Non Zero Off diagonal entries
    sel    <- matrix(((O[[i]] != 1) * (O[[i]] != 0)) == 1)
    numSel <- sum(sel)
    ### Computes the Coverage
    insCre       <- (uppQua[sel] >= O[[i]][sel]) * (lowQua[sel] <= O[[i]][sel])
    numCov       <- sum(insCre)
    cov          <- numCov / numSel
    sumLen       <- sum(uppQua[sel] - lowQua[sel])
    totSel       <- totSel + numSel
    totCov       <- totCov + numCov
    preCov[ii, i] <- cov
    preLen[ii, i] <- sumLen / numSel
    totLen        <- totLen + sumLen
  }
  preCov[ii, d + 1] <- totCov / totSel
  preLen[ii, d + 1] <- totLen / totSel
}

