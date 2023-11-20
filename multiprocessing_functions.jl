using PolaronMobility
using Statistics
using LinearAlgebra
using CSV
using Unitful
using Gnuplot
using Plots
using DataFrames


files = readdir("LiegeDataset/Results/GeneralizedFrohlich/conduction")
function general_val(data)
    return data[:, 1], data[:, 3], data[:, 4], data[:, 5], data[:, 6]
end

function trial(file)
    ħ = 1.054571817e-34
    kB = 1.380649e-23
    data = CSV.File("LiegeDataset/Results/GeneralizedFrohlich/conduction/" * file, delim = "\t") |> DataFrame
    mode, freq, res_alpha, res_paper, area = general_val(data)
    freq = res_paper ./ res_alpha
    non_zero_index = .!isnan.(freq) .& .!isinf.(freq) .& (freq .!= 0)
    #non_zero_index = Bool[0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1]
    #non_zero_index = Bool[0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1]
    for i in 1:length(res_alpha)
        if res_alpha[i] < 1e-8
            non_zero_index[i] = 0
        end
    end
    freq = freq[non_zero_index]
    res_alpha = res_alpha[non_zero_index]
    #res_paper = res_paper[non_zero_index]
    freq_actual = -freq * 0.2417990504024
    println(file)

    if length(freq) == 1
        #p = polaron(res_alpha[1], ω = freq_actual[1], β0 = ħ/kB*1e12*2π, v_guesses = 3.0001, w_guesses = 2.999)
        p = polaron(res_alpha[1], 0, ω = freq_actual[1], β0 = ħ/kB*1e12*2π, v_guesses = 3.0001, w_guesses = 2.999)
    elseif length(freq) == 0
        return
    else
        #p = polaron(res_alpha', ω = freq_actual, β0 = ħ/kB*1e12*2π, v_guesses = 3.0001, w_guesses = 2.999)
        p = polaron(res_alpha', 0, ω = freq_actual, β0 = ħ/kB*1e12*2π, v_guesses = 3.0001, w_guesses = 2.999)
    end
    addunits!(p)
    ZPR = p.F0 |> u"meV"
    ZPR = ustrip(ZPR)
    parts = split(file, '-')
    name = join(parts[1:2], '-')
        
    column_names = ["Name", "Frequency [meV]", "ZPR [meV]", "Reference_Alpha", "Reference_ZPR [meV]"]
    df_General = DataFrame([[name], [round(sum(freq_actual), digits=6)], [round(sum(ZPR), digits=6)], [round(sum(res_alpha), digits=6)], [round(sum(res_paper), digits=6)]], column_names)                
    #CSV.write("Data/multi_mode.tsv", df_General, delim='\t', quotechar='"', header=true)
    #CSV.write("Data/multi_mode_2/$name.tsv", df_General, delim='\t', quotechar='"',  header=true)
    #CSV.write("Data/multi_mode_zero_K/$name.tsv", df_General, delim='\t', quotechar='"',  header=true)
end

function trial_with_vw(file)
    ħ = 1.054571817e-34
    kB = 1.380649e-23
    data = CSV.File("LiegeDataset/Results/GeneralizedFrohlich/conduction/" * file, delim = "\t") |> DataFrame
    mode, freq, res_alpha, res_paper, area = general_val(data)
    freq = res_paper ./ res_alpha
    non_zero_index = .!isnan.(freq) .& .!isinf.(freq) .& (freq .!= 0)
    #non_zero_index = Bool[0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1]
    #non_zero_index = Bool[0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1]
    for i in 1:length(res_alpha)
        if res_alpha[i] < 1e-8
            non_zero_index[i] = 0
        end
    end
    freq = freq[non_zero_index]
    res_alpha = res_alpha[non_zero_index]
    #res_paper = res_paper[non_zero_index]
    freq_actual = -freq * 0.2417990504024
    println(file)

    if length(freq) == 1
        p = polaron(res_alpha[1], ω = freq_actual[1], β0 = ħ/kB*1e12*2π, v_guesses = 3.0001, w_guesses = 2.999)
        #p = polaron(res_alpha[1], 0, ω = freq_actual[1], β0 = ħ/kB*1e12*2π, v_guesses = 3.0001, w_guesses = 2.999)
    elseif length(freq) == 0
        return
    else
        p = polaron(res_alpha', ω = freq_actual, β0 = ħ/kB*1e12*2π, v_guesses = 3.0001, w_guesses = 2.999)
        #p = polaron(res_alpha', 0, ω = freq_actual, β0 = ħ/kB*1e12*2π, v_guesses = 3.0001, w_guesses = 2.999)
    end
    addunits!(p)
    ZPR = p.F0 |> u"meV"
    ZPR = ustrip(ZPR)
    parts = split(file, '-')
    name = join(parts[1:2], '-')
    v = p.v
    w = p.w
    column_names = ["Name", "Frequency [meV]", "ZPR [meV]", "Reference_Alpha", "Reference_ZPR [meV]", "v", "w"]
    df_General = DataFrame([[name], [round(sum(freq_actual), digits=6)], [round(sum(ZPR), digits=6)], [round(sum(res_alpha), digits=6)], [round(sum(res_paper), digits=6)], [round(v, digits=6)], [round(w, digits=6)]], column_names)                
    CSV.write("Data/multi_mode_vw_2/$name.tsv", df_General, delim='\t', quotechar='"',  header=true)
    #CSV.write("Data/multi_mode_zero_K/$name.tsv", df_General, delim='\t', quotechar='"',  header=true)
end