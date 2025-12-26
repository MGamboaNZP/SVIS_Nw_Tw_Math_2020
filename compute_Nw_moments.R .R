# =================================================================================
# compute_Nw_moments.R
#
# Compute the k-th order moments of Nw in a stochastic SIV model with imperfect
# vaccine and warning vaccination level.
#
# Based on:
# Gamboa, M., & Lopez-Herrero, M. J. (2020)
# "The Effect of Setting a Warning Vaccination Level on a Stochastic SIVS Model
# with Imperfect Vaccine." Mathematics, 8(7), 1136.
# https://doi.org/10.3390/math8071136
# =================================================================================

fqvi <- function(N, i, v, beta, si, gama) {
  valor <- (N - i - v) * ((beta * i / N) + si) +
    (beta * q * i * v / N) +
    (q * v * si) +
    (gama * i)
  return(valor)
}

tridiag <- function(upper, lower, main){
  n <- length(main)
  out <- matrix(0, n, n)
  diag(out) <- main
  indx <- seq.int(length(upper))
  out[cbind(indx + 1, indx)] <- lower
  out[cbind(indx, indx + 1)] <- upper
  return(out)
}

# --------------------------
# Model parameters
# --------------------------

N <- 100   # total population
h <- 45     # warning vaccination level
v0 <- 90    # initial vaccinated
i0 <- 1    # initial infected
beta <- 1.15
si <- 0.01
gama <- 1
q <- 0.1 # vaccine failure probability
kmax <- 4     # maximum order of moments to compute


# --------------------------
# Function to build tridiagonal matrix for a given vaccinated state
# --------------------------
matrizQvv <- function(v) {
  if (si == 0) {
    show("No external reinfection: model not valid")
  } else {
    c <- c()
    for (k in 1:(N - v)) {
      c[k] <- -k * gama
    }
    b <- c()
    for (k in 0:(N - v)) {
      b[k + 1] <- fqvi(N, k, v, beta, si, gama)
    }
    a <- c()
    for (k in 1:(N - v)) {
      a[k] <- -(N - v - k + 1) * ((beta * (k - 1) / N) + si)
    }
    QVV <- tridiag(a, c, b)
  }
  return(QVV)
}

# --------------------------
# Initialize matrices
# --------------------------
matrizKmenos1 <- matrix(1, nrow = N + 1, ncol = N + 1)

# Containers for first two moments
matrizesperanzas <- NULL
matrizmomentofact2 <- NULL

# --------------------------
# Compute k-th order moments
# --------------------------
for (k in 1:kmax) {
  matrizK <- matrix(nrow = N + 1, ncol = N + 1)
  
  for (i in 1:(N + 1)) {
    for (j in 1:(h + 1)) {
      matrizK[i, j] <- 0
    }
  }
  
  bvk <- c()
  
  for (j in (h + 2):N) {
    for (i in 1:(N - j + 1)) {
      lambda <- ((((i-1) * beta) / N) + si) * (N - j - i + 2)
      if (j == (h + 2)) {
        bvk[i] <- k * lambda * matrizKmenos1[i + 1, j]
      } else {
        parentesis <- (((i-1) * beta) / N) + si
        bvk[i] <- (j - 1) * q * parentesis * matrizK[i + 1, j - 1] +
          k * lambda * matrizKmenos1[i, j]
      }
      i <- N - j + 2
      parentesis <- (((i - 1) * beta) / N) + si
      bvk[i] <- (j - 1) * q * parentesis * matrizK[i + 1, j - 1]
      if (j == (h + 2)) {
        bvk[i] <- 0
      }
    }
    QVV <- matrizQvv(j - 1)
    m <- solve(QVV, bvk)
    for (i in 1:(N - j + 3)) {
      matrizK[i, j] <- m[i]
    }
    bvk <- c()
  }
  
  j <- N + 1
  if (j == (h + 2)) {
    bvk <- 0
  } else {
    bvk <- (N * q * si * matrizK[2, j - 1])
  }
  matrizK[1, N + 1] <- bvk / fqvi(N, 0, N, beta, si, gama)
  
  if (k == 1) {
    matrizesperanzas <- matrizK
  }
  if (k == 2) {
    matrizmomentofact2 <- matrizK
  }
  
  matrizKmenos1 <- matrizK
}

# --------------------------
# Extract results
# --------------------------
esperanza <- matrizesperanzas[i0 + 1, v0 + 1]
momentofactorial2 <- matrizmomentofact2[i0 + 1, v0 + 1]
varianza <- momentofactorial2 + esperanza - esperanza^2
momentofactorialk <- matrizK[i0 + 1, v0 + 1]

# --------------------------
# Example output
# --------------------------
cat("Expected value E[Nw] =", esperanza, "\n")
cat("Second factorial moment =", momentofactorial2, "\n")
cat("Variance =", varianza, "\n")
cat("k-th factorial moment (k =", kmax, ") =", momentofactorialk, "\n")

