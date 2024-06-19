# Example data and function overlapsAny (assuming it checks overlap with some condition)
library(tidyverse)
r <- 10
cons <- data.frame(A = c(1:5), B = c(7:11), C = c(8:12))  # Example data frame with columns A, B, C

# Function to check overlap with any column
overlapsAny <- function(r, column) {
  # Example logic to check overlap (replace with your actual logic)
  r %in% column
}

# Number input (example: 3)
num <- 3

# Initialize empty variable to store result
clus <- "NA"

# Loop through letters from A to Z
for (i in 1:num) {
  letter <- LETTERS[i]  # Get the ith letter from A to Z
  if (overlapsAny(r, cons[[letter]])) {
    clus <- letter  # Set clus to the letter if overlap condition is met
    break          # Exit loop as soon as condition is met
  }
}

# Print the result
print(clus)
