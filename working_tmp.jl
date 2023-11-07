using PolaronMobility
using Statistics
using LinearAlgebra
using JSON
using CSV
using Unitful
using Gnuplot
using Plots
using DataFrames
include("supplimentary_functions.jl")
include("multimode_workflow.jl")
#%%
phonon_files = readdir("LiegeDataset/Repository/phonon")
mass_files = readdir("LiegeDataset/Repository/eff_masses")
mass_abinit_files = readdir("LiegeDataset/Repository/abinit_effmasses")

#%%

phonon_data_1 = JSON.parse(read("LiegeDataset/Repository/phonon/" * phonon_files[1], String))
#%%
structure_string = phonon_data_1["metadata"]["structure"]
split_structure = split(structure_string)
 # Split and remove leading/trailing whitespaces

# Find the index of "_cell_volume" in the split elements
where = findfirst(isequal("_cell_volume"), split_structure)

# Check if "_cell_volume" was found
if where !== nothing
    # Convert the element at where + 1 to a floating-point number
    result = parse(Float64, split_structure[where + 1])

    # Print the result
    println(result)
else
    println("'_cell_volume' not found in the structure.")
end
# Print the result
#%%




#%%
mpid = "mp-22922"
#mpid = "mp-66"
what = 2
print(get_mstar(mpid, what))
#print(read_abinit_effmasses("LiegeDataset/Repository/abinit_effmasses/$mpid" * "_output")[what])
#%%
data = CSV.File("LiegeDataset/Results/StandardFeynmanFrohlich/conduction/standard_feynman_data_conduction.csv") |> DataFrame
println(names(data))
position = findfirst(x -> x == rpad(mpid, 21), data[:,1])
println(position)
#%%
mass = data[position,6]
#%%
print(data[4,1])
#mass = data[position]["MSTAR"]
#%%

files = readdir("LiegeDataset/Results/GeneralizedFrohlich/conduction")
println(files[2084])
#%%
data = CSV.File("LiegeDataset/Results/GeneralizedFrohlich/conduction/" * files[1397], delim = "\t") |> DataFrame
println(data[:,2])
#%%
(0.022822 * 4 * pi)^(-2) 
(0.008419 * 4 * pi)^(-2)
(0.007460 * 4 * pi * 1.240512)^(-2)
(0.009067 * 4 * pi)^(-2)

#%%
0.188836 / 0.025282
0.926796 / 0.016920
543.14425 / 7.653690
36.606386 / 0.351985
#%%
(0.175198 * 4 * pi)^(-2) 

print(0 / 0)
#%%
data = CSV.File("multi_mode.tsv") |> DataFrame
Error = (data[:,3] - data[:,5])./data[:,5] * 100
column_name = "Error"  # Replace with the actual column name
data.new_column = Error
CSV.write("multi_mode.tsv", DataFrame(data), delim='\t', append=false)
