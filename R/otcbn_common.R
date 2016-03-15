
from_to_values <- function(from, to) {
  seq(from, to , length.out = max(0, to-from+1))  
}

my.topological.sort <- function(poset) {
  p = ncol(poset)
  sorted.list = rep(0, p)
  for(i in 1:p) {
    sorted.list[i] = min.degree = which.min(apply(poset, 2, sum))  
    poset[min.degree, ] = 0
    poset[, min.degree] = Inf
  }
  
  sorted.list
  
}


# poset = matrix(0, p, p)
# 
# poset[2, 1] = 1
# poset[1, 3] = 1
# 
# my.topological.sort(poset)

mixedRadixGeneration <- function (d, m)
{
  N <- prod(m)
  T <- matrix(0, N, d)
  a <- matrix(0, 1, d + 1)
  m <- cbind(2, m)
  for (i in 1:N) {
    T[i, ] = a[2:(d + 1)]
    j <- d
    while (a[j + 1] == m[j + 1] - 1) {
      a[j + 1] <- 0
      j <- j - 1
    }
    a[j + 1] = a[j + 1] + 1
  }
  return(T)
}


hyperCube <- function (n)
{
  m <- matrix(2, nrow = 1, ncol = n)
  H <- mixedRadixGeneration(n, m)
  return(H)
}


# The source code is copied partly from icbn package

#' orderIdeals
#' @export
orderIdeals <- function (poset)
{
  H <- hyperCube(ncol(poset))
  N_H <- nrow(H)
  p <- ncol(H)
  G <- matrix(1, nrow = N_H, ncol = 1)
  for (j1 in 1:p) {
    for (j2 in 1:p) {
      if (poset[j1, j2]) {
        G <- as.numeric(G & (H[, j1] >= H[, j2]))
      }
    }
  }
  list(G=G, H=H, lattice_size=sum(G==1))
}


#' compatible_genotypes
#' @export
compatible_genotypes <- function (X, poset)
{
  # H <- hyperCube(p)
  G <- orderIdeals(poset)$G
  
  compatible_indexes = c()
  N <- nrow(X)
  p <- ncol(X)
  E <- p - (1:p)
  N_H <- 2^p
  D <- rep(0, N_H)
  for (i in 1:N) {
    idx = 1 + sum(X[i, ] * 2^E)
    D[idx] = D[idx] + 1
    
    if(G[idx] == 1)
    compatible_indexes <- c(compatible_indexes, i)
  }

  list( fraction=sum(D[as.logical(G)])/sum(D), compatible_indexes=compatible_indexes)
}