###############################################################################
###
### Simulated Data Generation for Tensor Graphical Models
###
###############################################################################

### Set-up
set.seed(240220)
devtools::load_all(".")
library(Matrix)
library(Tlasso)

### Parameter Settings
ii <- 1                     # Replication Number
p <- c(10, 10, 10)          # Matrix size for each dimension
d <- length(p)              # Number of tensor dimensions
r <- c(0.30, 0.20, 0.10)    # Sparity level for Random Matrices Only
b <- 103                    # Degrees of Freedom for Precision Matrices
n <- 10                     # Number of Observations
nmcmc <-100                 # Number of MCMC samples
burnin <- 10                # Burn-In Period
### General Suffix
suffix <- paste0(p[1], "-", p[2], "-", p[3], "-",
                 r[1], "-", r[2], "-", r[3], "-",
                 n, "-", nmcmc, "-", burnin, "-", ii)

### Generates the Precision, Covariance and Adjacency matrix
PSE <- spaPreMatGen(p        = p,
                    type     = "R",
                    sparsity = r)
### Asigns each Matrix List
O <- PSE$P
S <- PSE$S
E <- PSE$E

### Generates the Samples
print("Creating Samples")
### Estimation Sample
### As A List
lT <- TNormalSampler(n         = n,
                     SigmaList = S)
### As an Array
TT <- array(data = unlist(lT),
            dim  =c(p, n))
### Test Sample
### As A List
testlT <- TNormalSampler(n         = n,
                         SigmaList = S)
### As an Array
testTT <- array(data = unlist(testlT),
                dim  =c(p, n))

### Runs the Algorithm
print("Performing Bayesian Tensor Gaussian Graphical Model")
timBay <- Sys.time()
BFit <- BTGM(lT,
             burnin = 0,
             nmcmc  = nmcmc)
timBay <- difftime(time1 = Sys.time(),
                   time2 = timBay,
                   units = "secs")

### Runs TLasso
### Computes the Penalizing Parameter
CC     <- 20
lambda <- CC * sqrt(log(p) / (n * prod(p) / p))
print("Performing Tensor Lasso Model")
timFre <- Sys.time()
CFit <- Tlasso.fit(data       = TT,
                   T          = 1,
                   lambda.vec = lambda)
timFre <- difftime(time1 = Sys.time(),
                   time2 = timFre,
                   units = "secs")

### Saves R Objects
filNam <- paste0("./out/BFit-", suffix, ".rds")
saveRDS(BFit, file = filNam)
filNam <- paste0("./out/CFit-", suffix, ".rds")
saveRDS(CFit, file = filNam)
filNam <- paste0("./out/O-", suffix, ".rds")
saveRDS(O, file = filNam)
filNam <- paste0("./out/E-", suffix, ".rds")
saveRDS(E, file = filNam)
filNam <- paste0("./out/lT-", suffix, ".rds")
saveRDS(lT, file = filNam)
filNam <- paste0("./out/testlT-", suffix, ".rds")
saveRDS(testlT, file = filNam)

BMFit <- list()
for(i in 1:d){
  BMFit[[i]] <- apply(BFit$samC[[i]], c(1,2), median)
  BMFit[[i]] <- cov2cor(BMFit[[i]])
  BMFit[[i]] <- Matrix(data   = BMFit[[i]],
                       sparse = TRUE)
}

for(i in 1:d){
  CFit[[i]] <- cov2cor(CFit[[i]])
  CFit[[i]] <- Matrix(data   = CFit[[i]],
                      sparse = TRUE)
}

### Model Lists
models      <- list()
models[[1]] <- CFit
models[[2]] <- BMFit
models[[3]] <- O

print("Computing Precision Statistics")
preSta <- modComPre(modelList = models,
                    trueModel = models[[3]])

print("Computing Likelihood")
likSta <- modComLik(modelList = models,
                    tensors   = testlT)

BMFit <- list()
for(i in 1:d){
  BMFit[[i]] <- apply(BFit$samE[[i]], c(1,2), median)
  BMFit[[i]] <- Matrix(data   = BMFit[[i]],
                       sparse = TRUE)
}

### Computes the Adjacency Matrix of
for(j in 1:d){
  models[[1]][[j]] <- models[[1]][[j]] != 0
}

models[[2]] <- BMFit
models[[3]] <- BMFit

print("Computing Adjacency Statistics")
adjSta <- modComAdj(modelList = models,
                    trueModel = models[[3]])

gc()

timSta <- c(timFre, timBay, 0)

sta <- cbind(preSta$kroFro, preSta$kroInf, preSta$sinFro, preSta$sinInf,
             adjSta$acc, adjSta$tpr, adjSta$tnr, likSta, timSta)

sta <- data.frame(sta)
colnames(sta) <- c("Norm K-F", "Norm K-I", "Avg. Norm F", "Avg. Norm I",
                   "Accuracy", "TPR", "TNR", "Log_Lik", "Time")
row.names(sta) <- c("TLasso", "TBGGM", "Truth")


### Saves the Information
### Saves Performance Table
filNam <- paste0("./out/res-", suffix, ".txt")
print(sta)
write.table(x    = sta,
            file = filNam)
filNam <- paste0("./out/sta-", suffix, ".rds")
saveRDS(sta, file = filNam)

### Removes Variables
# remove(b)
# remove(d)
# remove(i)
# remove(n)
# remove(p)
# remove(r)
# remove(CC)
# remove(lambda)
# remove(j)
