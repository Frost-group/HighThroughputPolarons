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
mpid = "mp-5175"
#mpid = "mp-66"
what = 2
print(get_mstar(mpid, what))
#print(read_abinit_effmasses("LiegeDataset/Repository/abinit_effmasses/$mpid" * "_output")[what])