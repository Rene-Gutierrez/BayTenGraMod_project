### Graphical Lasso Iterator

library(glasso)

### Tlasso Replication
DFit <- list()
t <- 1
for (i in 1:d) {
  DFit[[i]] <- diag(p[i])
}
for (j in 1:t) {
  for (i in 1:d) {
    sK <- Sk(tensors    = lT,
             precisions = DFit,
             k          = i)

    sK <- sK / n / prod(p[-i])

    out <- huge::huge(x = sK,
                      lambda = lambda[i],
                      method = "glasso",
                      verbose = FALSE)
    DFit[[i]] <- as.matrix(out$icov[[1]]) / norm(as.matrix(out$icov[[1]]),
                                                 type = "F")
  }
}

### Tlasso Replication Constant Norm
t <- 1
DFit <- list()
for (i in 1:d) {
  DFit[[i]] <- diag(p[i])
}
dNor <- norm(DFit[[1]], type = 'F')
for (j in 1:t) {
  for (i in 1:d) {
    sK <- Sk(tensors    = lT,
             precisions = DFit,
             k          = i)

    sK <- sK / n / prod(p[-i])

    out <- huge::huge(x = sK,
                      lambda = lambda[i],
                      method = "glasso",
                      verbose = FALSE)
    DFit[[i]] <- as.matrix(out$icov[[1]]) / norm(as.matrix(out$icov[[1]]),
                                                 type = "F") * dNor
  }
}
CFit <- DFit

### Oracle
t <- 1
EFit <- list()
for (i in 1:d) {
  EFit[[i]] <- diag(p[i])
}
dNor <- norm(EFit[[1]], type = 'F')
for (j in 1:t) {
  for (i in 1:d) {
    sK <- Sk(tensors    = lT,
             precisions = O,
             k          = i)

    sK <- matrix(sK / n / prod(p[-i]), p[i], p[i])

    out <- huge::huge(x = sK,
                      lambda = lambda[i],
                      method = "glasso",
                      verbose = FALSE)
    EFit[[i]] <- as.matrix(out$icov[[1]]) / norm(as.matrix(out$icov[[1]]),
                                                 type = "F")
  }
}

CFit <- EFit

### Adapted Oracle
t <- 1
EFit <- list()
for (i in 1:d) {
  EFit[[i]] <- diag(p[i])
}
dNor <- norm(EFit[[1]], type = 'F')
for (j in 1:t) {
  for (i in 1:d) {
    sK <- Sk(tensors    = lT,
             precisions = O,
             k          = i)

    sK <- matrix(sK / n / prod(p[-i]), p[i], p[i])

    out <- huge::huge(x = sK,
                      lambda = lambda[i],
                      method = "glasso",
                      verbose = FALSE)
    EFit[[i]] <- as.matrix(out$icov[[1]]) / norm(as.matrix(out$icov[[1]]),
                                                 type = "F") * dNor
  }
}

CFit <- EFit

### Oracle with adaptative lambda
t <- 1
FFit <- list()
for (i in 1:d) {
  FFit[[i]] <- diag(p[i])
}
dNor <- norm(FFit[[1]], type = 'F')
for (j in 1:t) {
  for (i in 1:d) {
    sK <- Sk(tensors    = lT,
             precisions = O,
             k          = i)

    sK <- matrix(sK / n / prod(p[-i]), p[i], p[i])

    out <- huge::huge(x = sK,
                      method = "glasso",
                      verbose = FALSE)
    FFit[[i]] <- as.matrix(out$icov[[1]]) / norm(as.matrix(out$icov[[1]]),
                                                 type = "F")
  }
}

CFit <- FFit
