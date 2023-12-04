using PolaronMobility
using Statistics
using LinearAlgebra
using CSV
using Unitful
using Gnuplot
using DataFrames
using Plots
#%%

mobility_data = data = CSV.File("Data/single_mode_mobility_1.tsv") |> DataFrame

material_data = CSV.File("Data/material_data.tsv") |> DataFrame
mass_data = [material_data[:, "Name"], material_data[:, "mass"]]

merged_df = innerjoin(mobility_data, material_data, on=:Name)
#%%
plot(merged_df[:,"Alpha"], merged_df[:,"Mobility"], line=:scatter, xscale=:log10, yscale=:log10, zcolor=merged_df[:,"mass"])
#%%
CSV.write("Data/merged_files/mobility_1_material.tsv", merged_df, delim='\t', quotechar='"', header=true)
#%%
