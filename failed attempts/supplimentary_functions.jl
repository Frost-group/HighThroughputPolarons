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

function geometric_average(arr)
    n = length(arr)
    prod_val = prod(arr)
    geometric_avg = prod_val^(1/n)
    return geometric_avg
end

function dielectric_sort(dielectric)
    determinant_value = det(hcat(dielectric...))
    final_dielectric = cbrt(determinant_value)
    return final_dielectric
end
function read_abinit_effmasses(filename)
    effmass_flag = false
    val_idx = -1
    cond_idx = -1
    flag_occ = true
    band = 0
    deg_bands = []
    val_bands = []
    val_tensor = []
    val_eigen = []
    cond_bands = []
    cond_tensor = []
    cond_eigen = []
    open(filename, "r") do f
        for lines in eachline(f)
            println("cool", deg_bands, cond_idx)
            if occursin("  occ  ", lines) && flag_occ
                is_number = true
                finished = false
                l_split = split(lines)[2:end]
                val_idx = 0
                println(l_split)
                while is_number && !finished
                    for ls in l_split
                        if parse(Float64, ls) > 0.0
                            val_idx += 1
                        else
                            finished = true
                            break
                        end
                    end
                    lines = readline(f)
                    l_split = split(lines)
                    try
                        is_number = parse(Bool, l_split[1])
                    catch
                        is_number = false
                    end
                end
                cond_idx = val_idx + 1
            end
            if occursin("CALCULATION OF EFFECTIVE MASSES", lines)
                effmass_flag = true
            end
            if occursin("END OF EFFECTIVE MASSES SECTION", lines)
                effmass_flag = false
            end
            if effmass_flag
                if occursin("At k-point", lines)
                    
                    l_split = split(lines)
                    deg_bands = parse(Int, l_split[end-2]):parse(Int, l_split[end])
                    #println("cool", deg_bands, cond_idx)
                    if val_idx in deg_bands
                        val_tensor = zeros(length(deg_bands), 3, 3)
                        val_eigen = zeros(length(deg_bands), 3)
                        val_bands = deg_bands
                        cond_bands = []
                    elseif cond_idx in deg_bands
                        cond_tensor = zeros(length(deg_bands), 3, 3)
                        cond_eigen = zeros(length(deg_bands), 3)
                        cond_bands = deg_bands
                        println(band, cond_bands)
                    end
                end
                if occursin("K-point", lines)
                    l_split = split(lines)
                    band = parse(Int, l_split[end])
                    if val_idx == band && !(:val_tensor in names(Main))
                        val_tensor = zeros(1, 3, 3)
                        val_eigen = zeros(1, 3)
                        val_bands = [band]
                        cond_bands = []
                    end
                    if cond_idx == band && !(:cond_tensor in names(Main))
                        cond_tensor = zeros(1, 3, 3)
                        cond_eigen = zeros(1, 3)
                        cond_bands = [band]
                        println(band, cond_bands)
                    end
                end
                #println(band)
                if occursin("mass tensor", lines) && occursin(":", lines)
                    if occursin("eigenvalues", lines)
                        mdir = 1
                    else
                        mdir = 3
                    end
                    for idir in 1:mdir
                        lines = readline(f)
                        l_split = split(lines)

                        if band in val_bands
                            idx = findfirst(x -> x == band, val_bands)
                            if occursin("SADDLE", lines)
                                val_tensor[idx, :, :] .= [NaN, NaN, NaN]
                                val_eigen[idx, :] .= [NaN, NaN, NaN]
                                break
                            end
                            if mdir == 3
                                val_tensor[idx, idir, :] .= [parse(Float64, l_split[1]), parse(Float64, l_split[2]), parse(Float64, l_split[3])]
                            elseif mdir == 1
                                val_eigen[idx, :] .= [parse(Float64, l_split[1]), parse(Float64, l_split[2]), parse(Float64, l_split[3])]
                            end
                        elseif band in cond_bands
                            idx = findfirst(x -> x == band, cond_bands)
                            if occursin("SADDLE", lines)
                                cond_tensor[idx, :, :] .= [NaN, NaN, NaN]
                                cond_eigen[idx, :] .= [NaN, NaN, NaN]
                                break
                            end
                            if mdir == 3
                                cond_tensor[idx, idir, :] .= [parse(Float64, l_split[1]), parse(Float64, l_split[2]), parse(Float64, l_split[3])]
                            elseif mdir == 1
                                cond_eigen[idx, :] .= [parse(Float64, l_split[1]), parse(Float64, l_split[2]), parse(Float64, l_split[3])]
                            end
                        end
                    end
                end
            end
        end
    end

    return [[-val_tensor, abs.(val_eigen), val_bands], [cond_tensor, cond_eigen, cond_bands]]
end

function get_mstar(mpid, what)
    m_tensor = zeros(3, 3)
    m_eigen = [-1.0, -1.0, -1.0]
    m_deg = [-1.0]
    m_star_s = false
    not_working = true
    opt = "all"
    if opt == "mp" || opt == "all"
        try
            data = JSON.parsefile("LiegeDataset/Repository/eff_masses/$mpid.json")
            if what == 1
                m_star_x, m_star_y, m_star_z = data["cond_eff_mass"][1]["n"]["300"]["1e+18"]
            else
                m_star_x, m_star_y, m_star_z = data["cond_eff_mass"][1]["p"]["300"]["1e+18"]
            end
            m_eigen = abs([m_star_x, m_star_y, m_star_z])
            m_tensor = [m_tensor]
            m_deg = [1]
            m_star_s = "mp"
            not_working = false
        catch
        end
    end

    if opt == "abinit" || (opt == "all" && not_working)
        try
            tensor, eigen, deg = read_abinit_effmasses("LiegeDataset/Repository/abinit_effmasses/$mpid" * "_output")[what]
            m_tensor, m_eigen, m_deg = tensor, eigen, deg
            m_star_s = "abinit"
            not_working = false
        catch
        end
    end

    if not_working
        println("$mpid not in any effective mass database")
        return -1
    else
        println("$mpid effective mass calculated from $m_star_s database")
        return m_tensor, m_eigen, m_deg, m_star_s
    end
end