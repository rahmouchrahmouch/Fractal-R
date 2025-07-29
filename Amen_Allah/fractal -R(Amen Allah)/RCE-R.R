
# Load data
A <- readMat("C:/Users/amena/OneDrive - ESPRIT/Bureau/fractal -R/A.mat")
F11 <- A$A[, 11]  # Extract the 11th column from the matrix A

# Save F11 data for later use if needed
save(F11, file='F11.RData')

# Set cascade steps and parameters
nstep <- 7
tres <- 900  # Time resolution of the data in seconds
nvc <- 3     # Volume classes
boxs <- c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)
boxt <- boxs * tres
cst <- sapply(1:nstep, function(cs) mean(boxt[cs:(cs+1)]))

# Prepare data structure
data <- matrix(NA, nrow = length(F11), ncol = nstep + 1)
data[, 1] <- F11
nr <- nrow(data)

# Cascade evaluation
for (cs in 1:nstep) {
  npar <- nr / (2^cs)
  for (ap in 1:npar) {
    ao <- ap * 2
    data[ap, cs + 1] <- data[ao - 1, cs] + data[ao, cs]
  }
}

# Initialize 'parval' as a list with 4 empty lists (for the 4 parent categories)
parval <- vector("list", 4)
parinfo <- vector("list", length(F11))  # Initialize parinfo

for (cs in 1:nstep) {
  npar <- nr / (2^cs)
  parn <- rep(0, 4)
  
  for (ap in 2:(npar - 1)) {
    if (data[ap, cs + 1] > 0) {
      if (data[ap - 1, cs + 1] == 0 && data[ap + 1, cs + 1] == 0) {
        parn[1] <- parn[1] + 1
        parval[[1]] <- c(parval[[1]], data[ap, cs + 1])  # Append to list for parent category 1
        parinfo[[ap]] <- c(1, NA)
      } else if (data[ap - 1, cs + 1] == 0 && data[ap + 1, cs + 1] > 0) {
        parn[2] <- parn[2] + 1
        parval[[2]] <- c(parval[[2]], data[ap, cs + 1])  # Append to list for parent category 2
        parinfo[[ap]] <- c(2, NA)
      } else if (data[ap - 1, cs + 1] > 0 && data[ap + 1, cs + 1] > 0) {
        parn[3] <- parn[3] + 1
        parval[[3]] <- c(parval[[3]], data[ap, cs + 1])  # Append to list for parent category 3
        parinfo[[ap]] <- c(3, NA)
      } else {
        parn[4] <- parn[4] + 1
        parval[[4]] <- c(parval[[4]], data[ap, cs + 1])  # Append to list for parent category 4
        parinfo[[ap]] <- c(4, NA)
      }
      
      ao <- ap * 2
      if (data[ao - 1, cs] == 0) {
        parinfo[[ap]][3] <- 1
      } else if (data[ao, cs] == 0) {
        parinfo[[ap]][3] <- 2
      } else {
        parinfo[[ap]][3] <- 3
        parinfo[[ap]][4] <- data[ao - 1, cs] / data[ap, cs + 1]
      }
    }
  }
}

# Based on three volume classes, make volume limits
lim <- matrix(NA, 4, 2)
for (i in 1:4) {
  lim[i, ] <- quantile(parval[[i]], probs = c(0.33, 0.67))
}

# Add volume class in parinfo
for (pc in 1:4) {
  parpos <- which(sapply(parinfo, function(x) !is.null(x) && x[1] == pc))
  for (pp in parpos) {
    if (!is.na(data[pp, cs + 1]) && !is.na(lim[pc, 1]) && !is.na(lim[pc, 2])) {
      if (data[pp, cs + 1] <= lim[pc, 1]) {
        parinfo[[pp]][2] <- 1
      } else if (data[pp, cs + 1] <= lim[pc, 2]) {
        parinfo[[pp]][2] <- 2
      } else {
        parinfo[[pp]][2] <- 3
      }
    } else {
      parinfo[[pp]][2] <- NA  # Handle missing values
    }
  }
}

# Probability calculation (p-statistics)
prob <- array(NA, dim = c(nstep, 4, nvc, 3))
for (cs in 1:nstep) {
  for (pc in 1:4) {
    ppos <- which(sapply(parinfo, function(x) !is.null(x) && x[1] == pc))
    pinfo <- do.call(rbind, lapply(parinfo[ppos], function(x) x))
    
    for (vc in 1:nvc) {
      vpos <- which(pinfo[, 2] == vc)
      pvinfo <- pinfo[vpos, ]
      pvn <- nrow(pvinfo)
      
      for (dt in 1:3) {
        n <- length(which(pvinfo[, 3] == dt))
        prob[cs, pc, vc, dt] <- n / pvn
      }
    }
  }
}

# Saving the results
save(prob, file = "probabilities.RData")

# Plotting the probabilities
plot_pxx <- function(prob, nvc) {
  for (pc in 1:4) {
    for (dt in 1:3) {
      pplot <- matrix(NA, nvc, nstep)
      for (vc in 1:nvc) {
        pplot[vc, ] <- prob[, pc, vc, dt]
      }
      matplot(t(pplot), type = "l", main = paste("PC:", pc, "DT:", dt),
              xlab = "Cascade Step", ylab = "Probability")
      legend("topright", legend = paste("VC", 1:nvc), col = 1:nvc, lty = 1:nvc)
    }
  }
}

plot_pxx(prob, nvc)

# Histogram fitting and Beta distribution fitting
for (pc in 1:4) {
  if (length(parval[[pc]]) > 0) {
    # Step 1: Histogram Fitting
    hist(parval[[pc]], breaks = 20, main = paste("Histogram for Parent Category", pc),
         xlab = "Values", col = "lightblue", border = "black")
    
    # Step 2: Scale data for Beta Distribution Fitting
    scaled_data <- (parval[[pc]] - min(parval[[pc]])) / (max(parval[[pc]]) - min(parval[[pc]]))
    
    # Fit Beta distribution
    beta_fit <- fitdist(scaled_data, "beta", start = list(shape1 = 1, shape2 = 1))
    
    # Print summary of the fit
    print(summary(beta_fit))
    
    # Step 3: Plot histogram with fitted Beta distribution
    plot(beta_fit)
    title(main = paste("Beta Fit for Parent Category", pc))
    
    # Overlay the fitted Beta distribution
    x_fit <- seq(0, 1, length.out = 100)
    y_fit <- dbeta(x_fit, shape1 = beta_fit$estimate['shape1'], shape2 = beta_fit$estimate['shape2'])
    
    # Rescale to the original data range
    y_fit <- y_fit * (max(parval[[pc]]) - min(parval[[pc]])) * length(parval[[pc]]) / 100
    
    # Add the Beta distribution line to the histogram
    lines(y_fit, col = "red", lwd = 2)
  }
}
# Additional Histograms for Parent Categories 2 
for (pc in c(2, 3)) {  # Specify the parent categories for which you want additional histograms
  if (length(parval[[pc]]) > 0) {
    # Create a histogram for the specified parent category
    hist(parval[[pc]], breaks = 20, main = paste("Histogram for Parent Category", pc),
         xlab = "Values", col = "lightgreen", border = "black")
    
    # Step 1: Scale data for Beta Distribution Fitting
    scaled_data <- (parval[[pc]] - min(parval[[pc]])) / (max(parval[[pc]]) - min(parval[[pc]]))
    
    # Fit Beta distribution
    beta_fit <- fitdist(scaled_data, "beta", start = list(shape1 = 1, shape2 = 1))
    
    # Print summary of the fit
    print(summary(beta_fit))
    
    # Step 2: Plot histogram with fitted Beta distribution
    plot(beta_fit, main = paste("Beta Fit for Parent Category", pc), xlab = "Values")
    
    # Overlay the fitted Beta distribution
    x_fit <- seq(0, 1, length.out = 100)
    y_fit <- dbeta(x_fit, shape1 = beta_fit$estimate['shape1'], shape2 = beta_fit$estimate['shape2'])
    
    # Rescale to the original data range
    y_fit <- y_fit * (max(parval[[pc]]) - min(parval[[pc]])) * length(parval[[pc]]) / 100
    
    # Add the Beta distribution line to the histogram
    lines(y_fit, col = "red", lwd = 2)
  }
}

# Histogram and Beta Fit for Parent Category 3
pc <- 3  # Specify the parent category

if (length(parval[[pc]]) > 0) {
  # Create a histogram for Parent Category 3
  hist(parval[[pc]], breaks = 20, main = paste("Histogram for Parent Category", pc),
       xlab = "Values", col = "lightgreen", border = "black")
  
  # Step 1: Scale data for Beta Distribution Fitting
  scaled_data <- (parval[[pc]] - min(parval[[pc]])) / (max(parval[[pc]]) - min(parval[[pc]]))
  
  # Fit Beta distribution
  beta_fit <- fitdist(scaled_data, "beta", start = list(shape1 = 1, shape2 = 1))
  
  # Print summary of the fit
  print(summary(beta_fit))
  
  # Step 2: Plot histogram with fitted Beta distribution
  plot(beta_fit, main = paste("Beta Fit for Parent Category", pc), xlab = "Values")
  
  # Overlay the fitted Beta distribution
  x_fit <- seq(0, 1, length.out = 100)
  y_fit <- dbeta(x_fit, shape1 = beta_fit$estimate['shape1'], shape2 = beta_fit$estimate['shape2'])
  
  # Rescale to the original data range
  y_fit <- y_fit * (max(parval[[pc]]) - min(parval[[pc]])) * length(parval[[pc]]) / 100
  
  # Add the Beta distribution line to the histogram
  lines(y_fit, col = "red", lwd = 2)
}

