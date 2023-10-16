using PolaronMobility
using Statistics
using LinearAlgebra
using JSON

##
# Functions to load files
function load_file(dirs, files)
    dummy = []
    for i in files
    data = read(dirs * i, String)
    push!(dummy, JSON.parse(data))
    end
    return dummy
end

# Functions to load data
function JSON_read_mass(file, n)
    return file[n]["cond_eff_mass"][1]["n"]["300"]["1e+18"]
end

function JSON_read_die(file, n)
    return file[n]["dielectric"]["eps_total"]
end

function JSON_read_fre(file, n)
    return file[n]["phonon"]["ph_bandstructure"][1]
end
##

# Load data
atom = JSON.parsefile("LiegeDataset/Repository/atoms/atomic_data.json")
phonon_files = readdir("LiegeDataset/Repository/phonon")
mass_files = readdir("LiegeDataset/Repository/eff_masses")

phonon = load_file("LiegeDataset/Repository/phonon/", phonon_files)
mass = load_file("LiegeDataset/Repository/eff_masses/", mass_files)
flush(stdout)

##

print(mass_files[46])

##
# Trials of data
mass1 = JSON_read_mass(mass, 1)

die1 = JSON_read_die(phonon, 1)

fre1 = JSON_read_fre(phonon, 151)

##

# Functions to sort out data
function mass_sort(mass)
    #final_mass = mean(mass)
    final_mass = prod(mass)^(1/length(mass))

    return final_mass
end

function die_sort(die)
    determinant_value = det(hcat(die...))
    final_die = cbrt(determinant_value)

    return final_die
end

function fre_sort(fre)
    
    fre_new = transpose(reshape(fre, 3, :))
    println(mean(fre_new, dims=2))
    final_fre = maximum(mean(fre_new, dims=2))

    return final_fre
end

##
die_final = die_sort(die1)
##

# Functions to loop the data
function loop_mass(data)
    list_of_data = []
    for i in 1:length(data)

        if !isempty(mass[i]["cond_eff_mass"][1])
            mass_i = JSON_read_mass(mass, i)
            mean_mass_i = mass_sort(mass_i)
            append!(list_of_data, mean_mass_i)
            """
            Add saving data to csv functions
            """
        else
            println(i)
        end
    end
    return list_of_data
end

final_mass = loop_data(mass)
##

function loop_mass(data)
    list_of_data = []
    for i in 1:length(data)

        if !isempty(mass[i]["cond_eff_mass"][1])
            mass_i = JSON_read_mass(mass, i)
            mean_mass_i = mass_sort(mass_i)
            append!(list_of_data, mean_mass_i)
            """
            Add saving data to csv functions
            """
        else
            println(i)
        end
    end
    return list_of_data
end

final_mass = loop_data(mass)

##
function loop_mass(data)
    list_of_data = []
    for i in 1:length(data)

        if !isempty(mass[i]["cond_eff_mass"][1])
            mass_i = JSON_read_mass(mass, i)
            mean_mass_i = mass_sort(mass_i)
            append!(list_of_data, mean_mass_i)
            """
            Add saving data to csv functions
            """
        else
            println(i)
        end
    end
    return list_of_data
end

final_mass = loop_mass(mass)

##

function loop_die(data)
    list_of_data = []
    for i in 1:length(data)

        if !isempty(mass[i]["cond_eff_mass"][1])
            mass_i = JSON_read_mass(mass, i)
            mean_mass_i = mass_sort(mass_i)
            append!(list_of_data, mean_mass_i)
            """
            Add saving data to csv functions
            """
        else
            println(i)
        end
    end
    return list_of_data
end

final_mass = loop_mass(phonon)

##
println(phonon[2]["dielectric"])

##

println(phonon[80]["phonon"]["ph_bandstructure"][1])