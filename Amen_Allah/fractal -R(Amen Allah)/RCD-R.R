# Load necessary libraries
library(pracma)  # for percentile calculations

# Load model parameters
pfile <- readMat("C:/Users/amena/OneDrive - ESPRIT/Bureau/fractal -R/B11.mat")
c1 <- pfile$c1
c2 <- pfile$c2
c3 <- pfile$c3
c4 <- pfile$c4

# Define alim and slom (make sure to replace these with your actual values)
alim <- pfile$alim  # Load from .mat if available
# or set it manually
# alim <- some_value

slom <- pfile$slom  # Load from .mat if available
# or set it manually
# slom <- matrix(c(...), nrow = ..., ncol = ...)
ac <- pfile$ac 
# Initialize matrices
am <- array(0, dim = c(nstep, 4))  # a-values
pxxm <- array(0, dim = c(nstep, 4, nvc))  # model probabilities
p01m <- array(0, dim = c(nstep, 4, nvc))
p10m <- array(0, dim = c(nstep, 4, nvc))

# Define function to compute probabilities
compute_probabilities <- function(cs, parn, c1, c2, c3, c4, slom, alim) {
  for (vc in 1:nvc) {
    scale <- cstlog2[cs]
    int <- c1[parn, 3] + c2[parn, 3] * scale
    am[cs, parn] <- ifelse(scale > alim, 2^ac[parn], 2^(c3[parn] + c4[parn] * scale))
    
    pxxm[cs, parn, vc] <- int + slom[parn, 3] * vc
    pxxm[cs, parn, vc] <- ifelse(pxxm[cs, parn, vc] > 1, 1, pxxm[cs, parn, vc])
    p01m[cs, parn, vc] <- (1 - pxxm[cs, parn, vc]) / 2
    p10m[cs, parn, vc] <- (1 - pxxm[cs, parn, vc]) / 2
  }
}

print(dim(c1))
# Compute probabilities for each position class
for (parn in 1:3) {
  compute_probabilities(cs, parn, c1, c2, c3, c4, slom, alim)
}


# Load necessary library
library(R.matlab)

# Read the .mat file
data <- readMat("C:/Users/amena/OneDrive - ESPRIT/Bureau/fractal -R/C11.mat")

# Print the structure of the data to understand its format
str(data)

# Assuming the relevant data is in a list under a specific name
matrix_data <- data$A  

# Check the number of rows and columns
nr <- nrow(matrix_data)  # Get number of rows
nc <- ncol(matrix_data)  # Get number of columns
print(paste("Number of rows (nr):", nr))  
print(paste("Number of columns (nc):", nc))  

# Validate the number of rows
if (is.null(nr) || nr <= 0) {
  stop("Error: Invalid number of rows in the extracted data.")
}

nstep <- 10  # Set the number of steps (replace with actual value)

# Ensure nstep is valid
if (nstep > nc) {
  stop("Error: nstep exceeds the number of columns in matrix_data.")
}

max_par_limit <- 10000  # Set a maximum limit for npar

# Preallocate parval matrix
parval <- matrix(0, max_par_limit, 4)

# Initialize data array
for (nt in 1:10) {
  for (cs in 1:nstep) {
    npar <- nr * (2^(cs - 1))  # Number of "parents"
    nkid <- nr * (2^cs)        # Number of "kids"
    
    # Debug: Print npar and nkid values
    print(paste("Iteration:", nt, "Step:", cs, "npar:", npar, "nkid:", nkid))
    
    # Enforce maximum limit for npar
    if (npar > max_par_limit) {
      npar <- max_par_limit
    }
    
    parn <- numeric(4)          # Types of parents
    parinfo <- matrix(0, npar, 2)  # Parent information
    
    # Loop through each parent to classify and assign values
    for (ap in 1:npar) {
      matrix_data[ap, cs] <- max(0, matrix_data[ap, cs])  # Ensure non-negative
      lef <- if (ap == 1) 0 else matrix_data[ap - 1, cs]
      rig <- if (ap == npar) 0 else matrix_data[ap + 1, cs]
      
      if (matrix_data[ap, cs] > 0) {
        if ((lef == 0) && (rig == 0)) {
          parn[1] <- parn[1] + 1
          parval[parn[1], 1] <- matrix_data[ap, cs]
          parinfo[ap, 1] <- 1  # Isolated
        } else if ((lef == 0) && (rig > 0)) {
          parn[2] <- parn[2] + 1
          parval[parn[2], 2] <- matrix_data[ap, cs]
          parinfo[ap, 1] <- 2  # Starting
        } else if ((lef > 0) && (rig > 0)) {
          parn[3] <- parn[3] + 1
          parval[parn[3], 3] <- matrix_data[ap, cs]
          parinfo[ap, 1] <- 3  # Enclosed
        } else {
          parn[4] <- parn[4] + 1
          parval[parn[4], 4] <- matrix_data[ap, cs]
          parinfo[ap, 1] <- 4  # Stopping
        }
      }
    }
    
    # Add volume class in parinfo
    for (pc in 1:4) {
      parpos <- which(parinfo[, 1] == pc)
      for (pp in 1:parn[pc]) {
        if (matrix_data[parpos[pp], cs] <= vlim[pc, 1]) {
          parinfo[parpos[pp], 2] <- 1
        } else if (matrix_data[parpos[pp], cs] <= vlim[pc, 2]) {
          parinfo[parpos[pp], 2] <- 2
        } else {
          parinfo[parpos[pp], 2] <- 3
        }
      }
    }
    
    # Disaggregation model
    for (ap in 1:npar) {
      pos1 <- (2 * ap) - 1
      pos2 <- 2 * ap
      
      if (ap <= nr) {  # Ensure we are within bounds
        if (matrix_data[ap, cs] == 0) {
          if (pos1 <= nr && cs + 1 <= nc) {
            matrix_data[pos1, cs + 1] <- 0
          }
          if (pos2 <= nr && cs + 1 <= nc) {
            matrix_data[pos2, cs + 1] <- 0
          }
        } else {
          pc <- parinfo[ap, 1]
          vc <- parinfo[ap, 2]
          ran <- runif(1)
          
          # Debug: Print random number and parameters
          print(paste("Random number (ran):", ran))
          print(paste("pc:", pc, "vc:", vc))
          
          if (matrix_data[ap, cs] == 1) {
            p01ms <- p01m[cs, pc, vc] / (p01m[cs, pc, vc] + p10m[cs, pc, vc])
            p10ms <- p10m[cs, pc, vc] / (p01m[cs, pc, vc] + p10m[cs, pc, vc])
            
            # Debug: Print values of p01ms and p10ms
            print(paste("p01ms:", p01ms, "p10ms:", p10ms))
            
            if (is.nan(p01ms) || is.nan(p10ms) || (p01m[cs, pc, vc] == 0 && p10m[cs, pc, vc] == 0)) {
              stop("Error: p01m or p10m contains NA or both are zero.")
            }
            
            if (ran < p01ms) {
              if (pos1 <= nr && cs + 1 <= nc) {
                matrix_data[pos1, cs + 1] <- 0
              }
              if (pos2 <= nr && cs + 1 <= nc) {
                matrix_data[pos2, cs + 1] <- 1
              }
            } else {
              if (pos1 <= nr && cs + 1 <= nc) {
                matrix_data[pos1, cs + 1] <- 1
              }
              if (pos2 <= nr && cs + 1 <= nc) {
                matrix_data[pos2, cs + 1] <- 0
              }
            }
          } else {
            if (ran < p01m[cs, pc, vc]) {
              if (pos1 <= nr && cs + 1 <= nc) {
                matrix_data[pos1, cs + 1] <- 0
              }
              if (pos2 <= nr && cs + 1 <= nc) {
                matrix_data[pos2, cs + 1] <- matrix_data[ap, cs]
              }
            } else if (ran < (p01m[cs, pc, vc] + p10m[cs, pc, vc])) {
              if (pos1 <= nr && cs + 1 <= nc) {
                matrix_data[pos1, cs + 1] <- matrix_data[ap, cs]
              }
              if (pos2 <= nr && cs + 1 <= nc) {
                matrix_data[pos2, cs + 1] <- 0
              }
            } else {
              ran2 <- rgamma(2, am[cs, pc])  # Ensure am is defined
              ran3 <- ran2 / sum(ran2)
              
              if (pos1 <= nr && cs + 1 <= nc) {
                matrix_data[pos1, cs + 1] <- round(ran3[1] * matrix_data[ap, cs])
              }
              if (pos2 <= nr && cs + 1 <= nc) {
                matrix_data[pos2, cs + 1] <- round(ran3[2] * matrix_data[ap, cs])
              }
              
              # Ensure we have non-negative values after rounding
              if (pos1 <= nr && cs + 1 <= nc && is.na(matrix_data[pos1, cs + 1])) {
                matrix_data[pos1, cs + 1] <- 1
                if (pos2 <= nr && cs + 1 <= nc) {
                  matrix_data[pos2, cs + 1] <- matrix_data[pos2, cs + 1] - 1
                }
              } else if (pos2 <= nr && cs + 1 <= nc && is.na(matrix_data[pos2, cs + 1])) {
                matrix_data[pos2, cs + 1] <- 1
                if (pos1 <= nr && cs + 1 <= nc) {
                  matrix_data[pos1, cs + 1] <- matrix_data[pos1, cs + 1] - 1
                }
              }
            }
          }
        }
      }
    }
  }
}

# Optional: Save or visualize the results
write.csv(matrix_data, "C:/Users/amena/OneDrive - ESPRIT/Bureau/fractal -R/D11.csv", row.names = FALSE)

    
    