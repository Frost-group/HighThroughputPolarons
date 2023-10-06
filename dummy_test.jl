using Statistics
using LinearAlgebra
# 1x3n vector
row_vector = [1, 2, 3, 4, 5, 6, 7, 8, 9]

# Reshape into a 3xn matrix
matrix = transpose(reshape(row_vector, :, 3))

# Print the result
println(matrix)


# Reshaped 3xn matrix
matrix = [1 2 3; 4 5 6; 7 8 9; 10 11 12]

# Calculate the mean along each column
means = mean(matrix, dims=2)

# Print the result
println("Mean of each element:", means)