# SVIS Stochastic Model R Code

This repository contains R scripts associated with the publication:

> **Gamboa, M., & Lopez-Herrero, M. J. (2020).**  
> The Effect of Setting a Warning Vaccination Level on a Stochastic SIVS Model with Imperfect Vaccine.  
> *Mathematics*, 8(7), 1136.  
> https://doi.org/10.3390/math8071136

## Repository Structure

- `R/compute_Nw_moments.R` : Computes the moments of the variable Nw.
- `R/compute_Nw_distribution.R` : Computes the probability mass function of Nw.
- `R/compute_Tw_moments.R` : Computes the moments of the time to reach the warning level (Tw).


## Example Usage

```r
# Load scripts
source("R/compute_Nw_moments.R")
source("R/compute_Nw_distribution.R")
source("R/compute_Tw_moments.R")


# Compute Nw moments and distribution
Nw_moments <- compute_Nw_moments(N=100, v0=90, h=80, beta=1.15, si=0.01, gama=1, q=0.1)
Nw_distribution <- compute_Nw_distribution(N=100, v0=90, h=80, beta=1.15, si=0.01, gama=1, q=0.1)

# Compute Tw moments 
Tw_moments <- compute_Tw_moments(N=100, v0=90, h=80, beta=1.15, si=0.01, gama=1, q=0.1)

