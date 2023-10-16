using Statistics
using LinearAlgebra
using PolaronMobility
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

##

MAPI = material(4.9, 24.1, 0.12, 2.25)

p = polaron(MAPI, verbose = true)

p_2 = polaron(MAPI, [1, 100, 200], [1,2,3,4], verbose = true)