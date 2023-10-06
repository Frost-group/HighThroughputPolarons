using PolaronMobility
using Statistics
using LinearAlgebra
##
MAPI = material(4.9, 24.1, 0.12, 2.25)

p = polaron(MAPI, verbose = true)

p_2 = polaron(MAPI, [1, 100, 200], [1,2,3,4], verbose = true)



##
using JSON

function load_file(dirs, files)
    dummy = []
    for i in files
    data = read(dirs * i, String)
    push!(dummy, JSON.parse(data))
    end
    return dummy
end

function JSON_read_mass(file, n)
    return file[n]["cond_eff_mass"][1]["p"]["300"]["1e+18"]
end

function JSON_read_die(file, n)
    return file[n]["dielectric"]["eps_total"]
end

function JSON_read_fre(file, n)
    return file[n]["phonon"]["ph_bandstructure"][1]
end
##

atom = JSON.parsefile("LiegeDataset/Repository/atoms/atomic_data.json")
phonon_files = readdir("LiegeDataset/Repository/phonon")
mass_files = readdir("LiegeDataset/Repository/eff_masses")

phonon = load_file("LiegeDataset/Repository/phonon/", phonon_files)
mass = load_file("LiegeDataset/Repository/eff_masses/", mass_files)
flush(stdout)

##
mass1 = JSON_read_mass(mass, 1)

die1 = JSON_read_die(phonon, 1)

fre1 = JSON_read_fre(phonon, 151)
##

function mass_sort(mass)
    final_mass = mean(mass)

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
#mass_final = mass_sort.(mass)