### Predictive edge

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
preSel <- matrix(NA,
                 nrow = iter,
                 ncol = d + 1)

Ohat   <- list()
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

  totMat <- 0
  totSel <- 0
  for(i in 1:d){
    Ohat[[i]] <- BFit$samC[[i]]
    matched   <- 0
    ### Selects the Zero entries
    sel    <- matrix(O[[i]] == 0)
    numSel <- sum(sel)
    for(s in 1:nmcmc){
      matched <- matched + sum(Ohat[[i]][,,s][matrix(O[[i]] == 0)] == 0)
    }
    preSel[ii, i] <- matched / nmcmc / numSel
    totMat        <- totMat + matched
    totSel        <- totSel + numSel
  }
  preSel[ii, d + 1] <- totMat / totSel / nmcmc
}
