using Statistics
using LinearAlgebra
using PolaronMobility
using Gnuplot
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

##
x = 1:0.1:10
@gp    "set grid" "set key right" "set logscale y"
@gp :- "set title 'Plot title'" "set label 'X label'" "set xrange [0:*]"
@gp :- x x.^0.5 "w l tit 'Pow 0.5' dt 2 lw 2 lc rgb 'red'"
@gp :- x x      "w l tit 'Pow 1'   dt 1 lw 3 lc rgb 'blue'"
@gp :- x x.^2   "w l tit 'Pow 2'   dt 8 lw 2 lc rgb 'purple'"
##

using DataFrames
# Sample data arrays
arr_a = [1, 2, 3, 4, 5]
arr_b = [6, 7, 8, 9, 10]
arr_c = [11, 12, 13, 14, 15]
arr_d = [16, 17, 18, 19, 20]
arr_e = [21, 22, 23, 24, 25]

# Create DataFrame without column names
df = DataFrame([arr_a, arr_b, arr_c, arr_d, arr_e], :auto)
##
