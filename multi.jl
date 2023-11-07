using Distributed
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
dummy_name = []
dummy_freq = []
dummy_alpha = []

dummy_ZPR = []
dummy_result = []

#%%
# Number of processes to spawn

num_processes = 24
start_time = time()
# Start local worker processes
addprocs(num_processes)


#%%
# Display elapsed time
# Function to be executed in parallel
#= @everywhere function my_parallel_function(id)
    println("Task $id is running on process ", myid())
    return id^2
end =#
function trial(file)
    ħ = 1.054571817e-34
    kB = 1.380649e-23
    data = CSV.File("LiegeDataset/Results/GeneralizedFrohlich/conduction/" * file, delim = "\t") |> DataFrame
    mode, freq, res_alpha, res_paper, area = general_val(data)
    freq = res_paper ./ res_alpha
    non_zero_index = .!isnan.(freq) .& .!isinf.(freq) .& (freq .!= 0)
    println(non_zero_index)
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
    println(non_zero_index)
    println(file, freq_actual, res_alpha, res_paper)

    if length(freq) == 1
        p = polaron(res_alpha[1], ω = freq_actual[1], β0 = ħ/kB*1e12*2π, v_guesses = 3.0001, w_guesses = 2.999)

    elseif length(freq) == 0
        return
    else
        p = polaron(res_alpha', ω = freq_actual, β0 = ħ/kB*1e12*2π, v_guesses = 3.0001, w_guesses = 2.999)
    end
    addunits!(p)
    ZPR = p.F0 |> u"meV"
    ZPR = ustrip(ZPR)
    parts = split(file, '-')
    name = join(parts[1:2], '-')
    push!(dummy_name, name)
    append!(dummy_alpha, sum(res_alpha))
    append!(dummy_freq, sum(freq))
    append!(dummy_result, ZPR)
    append!(dummy_ZPR, sum(res_paper))
        
    column_names = ["Name", "Frequency [meV]", "ZPR", "Reference_Alpha", "Reference_ZPR"]
    df_General = DataFrame([[name], [sum(freq)], [ZPR], [sum(res_alpha)], [sum(res_paper)]], column_names)                
    #CSV.write("Data/multi_mode.tsv", df_General, delim='\t', quotechar='"', header=true)
    CSV.write("Data/multi_mode.tsv", df_General, delim='\t', quotechar='"', append=true, header=false)
end

# Generate tasks

# Fetch results from each task

    # Code that may raise exceptions
tasks = [trial(file) for file in files[2084: end]]
results = fetch(tasks)
rmprocs(workers())


# Display results
elapsed_time = time() - start_time
println("Elapsed time: $elapsed_time seconds")
#%%
