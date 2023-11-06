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
#%%



function get_dielectric(file)
    eps_ionic = dielectric_sort(file["dielectric"]["eps_electronic"])
    eps_static = dielectric_sort(file["dielectric"]["eps_total"])
    return eps_ionic, eps_static
end

function loop_material()
    phonon_files = readdir("LiegeDataset/Repository/phonon")
    dummy_name = []
    dummy_eps_ionic = []
    dummy_eps_static = []
    dummy_volume = []
    dummy_mass = []
    for i in 1:length(phonon_files)
        
        phonon_data = JSON.parse(read("LiegeDataset/Repository/phonon/" * phonon_files[i], String))
        name = phonon_files[i][1:end-5]
        mass = get_mstar(name)
        if mass !== nothing
            eps_ionic, eps_static = get_dielectric(phonon_data)

            structure_string = phonon_data["metadata"]["structure"]
            split_structure = split(structure_string)
            # Split and remove leading/trailing whitespaces

            # Find the index of "_cell_volume" in the split elements
            where_ = findfirst(isequal("_cell_volume"), split_structure)
            # Check if "_cell_volume" was found
            if where_ !== nothing
                # Convert the element at where + 1 to a floating-point number
                volume = parse(Float64, split_structure[where_ + 1])

                # Print the result
            else
                println("'_cell_volume' not found in the structure.")
            end


            push!(dummy_name, name)
            append!(dummy_eps_ionic, eps_ionic)
            append!(dummy_eps_static, eps_static)
            append!(dummy_volume, float(volume))
            append!(dummy_mass, mass)
        end

    end
    column_names = ["Name", "mass", "eps_ionic", "eps_static", "volume"]
    df_material = DataFrame([dummy_name, dummy_mass, dummy_eps_ionic, dummy_eps_static, dummy_volume], column_names)                
    CSV.write("Data/material_data.tsv", df_material, delim='\t', quotechar='"', header=true)
    return df_material
end
#%%