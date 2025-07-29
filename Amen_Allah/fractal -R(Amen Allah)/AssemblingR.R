
# Load the .mat file
data <- readMat("C:/Users/amena/OneDrive - ESPRIT/Bureau/fractal -R/B11.mat")

# Extract the relevant matrix (assuming it's named 'data' in the .mat file)
data_matrix <- data$data  # Adjust this line if the variable name is different

# Get the number of rows and columns
rows <- nrow(data_matrix)
cols <- ncol(data_matrix)

# Ensure the number of rows is an integer
num_days <- floor(rows / 96)
new <- numeric(num_days)  # Initialize as a numeric vector

# Sum data for each day (every 96 entries)
for (k in 1:num_days) {
  start_index <- (k - 1) * 96 + 1
  end_index <- start_index + 95
  new[k] <- sum(data_matrix[start_index:end_index, 1])  # Sum for the first column
}

# Normalize by dividing by 0.0625
C11 <- new / 0.0625

# Save the result to a .RData file
save(C11, file = "C:/Users/amena/OneDrive - ESPRIT/Bureau/fractal -R/C11.RData")
print(C11)
# Save C11 to a CSV file
write.csv(C11, file = "C:/Users/amena/OneDrive - ESPRIT/Bureau/fractal -R/C11.csv", row.names = FALSE)

