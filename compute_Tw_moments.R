# =================================================================================
# compute_Tw_moments.R
#
# Compute the k-th order moments of Tw (time until warning vaccination level is reached)
# in a stochastic SIV model with imperfect vaccine.
#
# Based on:
# Gamboa, M., & Lopez-Herrero, M. J. (2020)
# "The Effect of Setting a Warning Vaccination Level on a Stochastic SIVS Model
# with Imperfect Vaccine." Mathematics, 8(7), 1136.
# https://doi.org/10.3390/math8071136
# =================================================================================

#---------------------------------------
# Auxiliary function: transition values
#---------------------------------------
fqvi <- function(N, i, v, beta, si, gama) {
  valor <- (N - i - v) * ((beta * i / N) + si) + (beta * q * i * v / N) + (q * v * si) + (gama * i)
  return(valor)
}

#---------------------------------------
# Auxiliary function: tridiagonal matrix
#---------------------------------------
tridiag <- function(upper, lower, main) {
  out <- matrix(0, length(main), length(main))
  diag(out) <- main
  indx <- seq.int(length(upper))
  out[cbind(indx + 1, indx)] <- lower
  out[cbind(indx, indx + 1)] <- upper
  return(out)
}

#=================== PARAMETERS ===================
h <- 70        # Warning vaccination level
N <- 100       # Total population
v0 <- 90       # Number of vaccinated individuals where we calculate Tw
beta <- 4     # Transmission rate
si <- 0.01     # External infection
gama <- 1      # Recovery rate
q <- 0.1       # Vaccine failure probability

#=================== TRANSITION MATRIX FUNCTION ===================
matrizQvv <- function(v) {
  if (si == 0) {
    show("Not an SIV model with external reinfection")
  } else {
    c <- c()
    for (k in 1:(N - v)) { c[k] <- -k * gama }
    b <- c()
    for (k in 0:(N - v)) { b[k + 1] <- fqvi(N, k, v, beta, si, gama) }
    a <- c()
    for (k in 1:(N - v)) { a[k] <- -(N - v - k + 1) * ((beta * (k - 1) / N) + si) }
    QVV <- tridiag(a, c, b)
  }
  return(QVV)
}

#=================== INITIALIZATION ===================
matrizesperanza <- matrix(1, nrow = N + 1, ncol = N + 1)
matrizvarianza <- matrix(1, nrow = N + 1, ncol = N + 1)
matrizKmenos1 <- matrix(1, nrow = N + 1, ncol = N + 1)

kmax <- 2  # Calculate first and second moments
for (k in 1:kmax) {
  matrizK <- matrix(nrow = N + 1, ncol = N + 1)
  for (i in 1:(N + 1)) {
    for (j in 1:(h + 1)) { matrizK[i, j] <- 0 }
  }
  bvk <- c()
  
  for (j in (h + 2):N) {
    for (i in 1:(N - j + 2)) {
      if (j == (h + 2)) {
        bvk[i] <- k * matrizKmenos1[i, j]
      } else {
        bvk[i] <- (j - 1) * q * ((((i-1) * beta) / N) + si) * matrizK[i + 1, j - 1] + k * matrizKmenos1[i, j]
      }
    }
    QVV <- matrizQvv(j - 1)
    m <- solve(QVV, bvk)
    for (i in 1:(N - j + 3)) { matrizK[i, j] <- m[i] }
    bvk <- c()
  }
  
  j <- N + 1
  if (j == (h + 2)) { bvk <- k * matrizKmenos1[1, j] } 
  else { bvk <- (N * q * si * matrizK[2, j - 1] + k * matrizKmenos1[1, j]) }
  matrizK[1, N + 1] <- bvk / fqvi(N, 0, N, beta, si, gama)
  
  if (k == 1) matrizesperanza <- matrizK
  if (k == 2) matrizvarianza <- matrizK
  
  matrizKmenos1 <- matrizK
}

#=================== RESULTS ===================
mvarianza <- matrizvarianza - matrizesperanza^2
i0 <- 1
esperanza <- matrizesperanza[i0 + 1, v0 + 1]
varianza <- mvarianza[i0 + 1, v0 + 1]

cat("Expected Tw (i0 =", i0, ", v0 =", v0, "):", esperanza, "\n")
cat("Variance of Tw (i0 =", i0, ", v0 =", v0, "):", varianza, "\n")

