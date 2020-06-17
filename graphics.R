### Graphic and Table Evaluation

### Sets Up Libraries
devtools::load_all(".")
library(Matrix)
library(ggplot2)

### Desired Set-up to Evaluate
p <- c(10, 10, 10)          # Matrix size for each dimension
d <- length(p)              # Number of tensor dimensions
r <- c(0.30, 0.20, 0.10)    # Sparity level for Random Matrices Only
b <- 103                    # Degrees of Freedom for Precision Matrices
n <- 10                     # Number of Observations
nmcmc <-100                 # Number of MCMC samples
burnin <- 10                # Burn-In Period
### General Suffix
ii <- 1
suffix <- paste0(p[1], "-", p[2], "-", p[3], "-",
                 r[1], "-", r[2], "-", r[3], "-",
                 n, "-", nmcmc, "-", burnin, "-", ii)

### Loads R Objects
filNam <- paste0("./out/BFit-", suffix, ".rds")
BFit   <- readRDS(file = filNam)
filNam <- paste0("./out/CFit-", suffix, ".rds")
CFit   <- readRDS(file = filNam)
filNam <- paste0("./out/O-", suffix, ".rds")
O      <- readRDS(file = filNam)
filNam <- paste0("./out/E-", suffix, ".rds")
E      <- readRDS(file = filNam)
filNam <- paste0("./out/sta-", suffix, ".rds")
sta    <- readRDS(file = filNam)

### Plots selected Pisteriors
sel <- list()
for (i in 1:d) {
  ### Gets the size of the Matrix
  p <- nrow(O[[i]])
  non <- sample(x    = seq(1:(p^2))[(E[[i]] * upper.tri(x    = E[[i]],
                                                        diag = FALSE)) == 1],
                size = 2)
  zer <- sample(x    = seq(1:p^2)[(!E[[i]]) * upper.tri(x    = E[[i]],
                                                        diag = FALSE) == 1],
                size = 2)
  sel[[i]] <- c(non, zer)
}

posTab <- datPos(posPre = BFit$samC,
                 reaPre = O,
                 sel    = sel)

posPlo <- ggplot(data    = posTab,
                 mapping = aes(x    = iter,
                               y    = sample)) +
  geom_line(data    = posTab,
            mapping = aes(x = iter, y = truth),
            color   = "red",
            size    = 0.5) +
  geom_linerange(data    = posTab,
                 mapping = aes(x    = iter,
                               ymin = pmin(sample, truth),
                               ymax = pmax(sample, truth))) +
  facet_wrap(~entry,
             scales = "free")

ggsave(filename = paste0("pos-", suffix, ".pdf"),
       plot     = posPlo,
       device   = "pdf",
       path     = "./out/")

# print(posPlo)

posPlo <- ggplot(data    = posTab,
                 mapping = aes(x    = iter,
                               y    = sample,
                               ymin = -1,
                               ymax =  1)) +
  geom_line(data    = posTab,
            mapping = aes(x = iter, y = truth),
            color   = "red",
            size    = 0.5) +
  geom_linerange(data    = posTab,
                 mapping = aes(x    = iter,
                               ymin = pmin(sample, truth),
                               ymax = pmax(sample, truth))) +
  facet_wrap(~entry)

ggsave(filename = paste0("pos-zoom-", suffix, ".pdf"),
       plot     = posPlo,
       device   = "pdf",
       path     = "./out/")

# print(posPlo)

### Adapts Results
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

prePlo <- ggplot(data = heatMat(ll         = models,
                                modelNames = c("TLasso", "TBGGM", "Truth")),
                 aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low   = "blue",
                       high   = "red",
                       mid    = "white",
                       limits = c(-1,1)) +
  facet_grid(rows = vars(matrix),
             cols = vars(model))

ggsave(filename = paste0("pre-", suffix, ".pdf"),
       plot     = prePlo,
       device   = "pdf",
       path     = "./out/")

# print(prePlo)

### Computes the Adjacency Matrix of
BMFit <- list()
for(j in 1:d){
  BMFit[[j]] <- apply(BFit$samE[[j]], c(1,2), median)
  BMFit[[j]] <- Matrix(data   = BMFit[[j]],
                       sparse = TRUE)
  models[[1]][[j]] <- models[[1]][[j]] != 0
  models[[3]][[j]] <- E[[j]] + t(E[[j]]) - diag(diag(E[[j]]))
}

models[[2]] <- BMFit


adjPlo <- ggplot(data = heatMat(ll         = models,
                                modelNames = c("TLasso", "TBGGM", "Truth")),
                 aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low   = "white",
                       high   = "blue",
                       limits = c(0,1)) +
  facet_grid(rows = vars(matrix),
             cols = vars(model))

ggsave(filename = paste0("adj-", suffix, ".pdf"),
       plot     = adjPlo,
       device   = "pdf",
       path     = "./out/")

# print(adjPlo)

print(sta)
