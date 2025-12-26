# =================================================================================
# compute_Nw_mass_function.R
#
# Compute the mass probability function of Nw in a stochastic SIV model with imperfect
# vaccine and warning vaccination level.
#
# Based on:
# Gamboa, M., & Lopez-Herrero, M. J. (2020)
# "The Effect of Setting a Warning Vaccination Level on a Stochastic SIVS Model
# with Imperfect Vaccine." Mathematics, 8(7), 1136.
# https://doi.org/10.3390/math8071136
# =================================================================================

# Function to compute the rate term in the SIV model
fqvi <- function(N, i, v, beta, si, gama){
  valor <- (N - i - v) * ((beta * i / N) + si) + (beta * q * i * v / N) + (q * v * si) + (gama * i)
  return(valor)
}

# Function to create a tridiagonal matrix
tridiag <- function(upper, lower, main){
  out <- matrix(0, length(main), length(main))
  diag(out) <- main
  indx <- seq.int(length(upper))
  out[cbind(indx + 1, indx)] <- lower
  out[cbind(indx, indx + 1)] <- upper
  return(out)
}

#----------------------------------------------
# Parameters
#----------------------------------------------
h <- 65          # Warning level
N <- 100         # Total population
beta <- 1.15     # Transmission rate
si <- 0.01       # External infection rate
gama <- 1        # Recovery rate
q <- 0.1         # Vaccine failure probability
v0 <- 90         # Initial vaccinated individuals for conditional probability
i0 <- 1          # Initial number of infected individuals

funcionprobabilidad <- c()  # Vector to store the probability function

#----------------------------------------------
# Function to construct the tridiagonal matrix Q for given v
#----------------------------------------------
matrizQvv <- function(v){
  if (si == 0){
    show("No es un Modelo SIV con reinfección externa")
  } else {
    # Lower diagonal
    c <- -gama * seq(1, N - v)
    # Main diagonal
    b <- sapply(0:(N - v), function(k) fqvi(N, k, v, beta, si, gama))
    # Upper diagonal
    a <- rep(0, N - v)
    # Construct tridiagonal matrix
    QVV <- tridiag(a, c, b)
  }
  return(QVV)
}

#----------------------------------------------
# Initialization of matrices
#----------------------------------------------
bvk <- c()
matrizKmenos1 <- matrix(nrow = N + 1, ncol = N + 1)
matrizKmenos1[, h + 1] <- 1   # Base case for k = 0

#----------------------------------------------
# Construct the probability matrix for k = 0
#----------------------------------------------
for (j in (h + 2):N){
  for (i in 1:(N - j + 2)){
    bvk[i] <- q * (j - 1) * (beta * (i - 1) / N + si) * matrizKmenos1[i + 1, j - 1]
  }
  QVV <- matrizQvv(j - 1)
  a <- solve(QVV, bvk)
  matrizKmenos1[1:length(a), j] <- a
  bvk <- c()
}
matrizKmenos1[1, N + 1] <- matrizKmenos1[2, N] * q * N * si / fqvi(N, 0, N, beta, si, gama)

#----------------------------------------------
# Iterative computation of the probability function
#----------------------------------------------
k <- 1
suma <- matrizKmenos1[i0 + 1, v0 + 1]
matrizK <- matrix(nrow = N + 1, ncol = N + 1)
matrizK[i0 + 1, v0 + 1] <- 0

while (suma < 0.999 & matrizK[i0 + 1, v0 + 1] >= 0) {
  
  # Reset the probability matrix
  matrizK[,] <- 0
  
  bvk <- c()
  
  # Compute contributions for each vaccination level
  for (j in (h + 2):N){
    for (i in 1:(N - j + 2)){
      if (i == (N - j + 2)){
        bvk[i] <- q * (j - 1) * (beta * (i - 1) / N + si) * matrizK[i + 1, j - 1]
      } else {
        bvk[i] <- q * (j - 1) * (beta * (i - 1) / N + si) * matrizK[i + 1, j - 1] +
          (beta * (i - 1) / N + si) * (N - (j - 1) - (i - 1)) * matrizKmenos1[i + 1, j]
      }
    }
    QVV <- matrizQvv(j - 1)
    mm <- solve(QVV, bvk)
    matrizK[1:length(mm), j] <- mm
    bvk <- c()
  }
  
  # Special case for the last column
  matrizK[1, N + 1] <- q * N * si * matrizK[2, N] / fqvi(N, 0, N, beta, si, gama)
  
  # Update matrices
  matrizKmenos1 <- matrizK
  funcionprobabilidad <- c(funcionprobabilidad, matrizK[i0 + 1, v0 + 1])
  suma <- sum(funcionprobabilidad)
  k <- k + 1
}

#----------------------------------------------
# Output the probability function
#----------------------------------------------
cat("Probability mass function for Nw (v0 =", v0, "):\n")
print(funcionprobabilidad)
