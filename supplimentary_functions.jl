using PolaronMobility
using Statistics
using LinearAlgebra
using JSON
using CSV
using Unitful
using Gnuplot
using Plots
using DataFrames

function load_file(file)
    """
    load_file(file)
    Reads the contents of a file and returns the data.

    # Arguments
    - `file` (string): The path to the file to be read.

    # Returns
    - `data` (string): The contents of the file.
    """
    open(file) do f
        data = read(f)
        return data
    end
    
end


function dielectric_sort(dielectric)
    determinant_value = det(hcat(dielectric...))
    final_dielectric = cbrt(determinant_value)
    return final_dielectric
end

function get_mstar(mpid)
    mass = 0
    try
        data = CSV.File("LiegeDataset/Results/StandardFeynmanFrohlich/conduction/standard_feynman_data_conduction.csv") |> DataFrame
        position = findfirst(x -> x == rpad(mpid, 21), data[:,1])
        mass = data[position,6]
    catch 
    end
    
    if mass == 0
        try
            data = CSV.File("LiegeDataset/Results/StandardFrohlich/conduction/standard_froelich_data_conduction.csv") |> DataFrame
            position = findfirst(x -> x == rpad(mpid, 21), data[:,1])
            mass = data[position,6]
        catch 
        end
    end

    if mass == 0
        println("ohno  ", mpid)
        return
        
    else
        return mass
    end
end