using PolaronMobility
using Statistics
using LinearAlgebra
using CSV
using Unitful
using Gnuplot
using Plots
using DataFrames
function load_file(file)
    open(file) do f
        data = read(f)
        return data
    end
    
end

function variable_cal(data)
    return data[1], data[2], data[3], data[4], data[5] * 0.2417990504024, data[6], data[7], data[8]
end


function looping(data)
    dummy_name = []
    dummy_chem = []
    dummy_alpha = []
    dummy_ZPR = []
    dummy_mass = []
    dummy_freq = []
    dummy_res = []
    dummy_alpha_paper = []
    for i in 1:length(data)
        name, chem, eps_s, eps_o, freq, mass, res_alpha, res_paper = variable_cal(data[i])
        mass = mass * mass
        mat = material(eps_o, eps_s, mass, freq)
        try
            p = polaron(mat, v_guesses = 3.0001, w_guesses = 2.999 ,verbose = false)
            addunits!(p)
            ZPR = p.F0 |> u"meV"
            alpha = p.α
            freq = freq / 0.2417990504024
            push!(dummy_name, name)
            push!(dummy_chem, chem)
            append!(dummy_alpha, alpha)
            append!(dummy_ZPR, ZPR)
            append!(dummy_mass, mass)
            append!(dummy_freq, freq)
            append!(dummy_res, res_paper)
            append!(dummy_alpha_paper, res_alpha)
        catch e
            # Catch and handle the error
            println("An error occurred: ", e, i, name)
        end
        end
    return dummy_name, dummy_chem, dummy_alpha, dummy_ZPR, dummy_mass, dummy_freq, dummy_alpha_paper, dummy_res
end

function Feynman_Data()
    data = CSV.File("LiegeDataset/Results/StandardFeynmanFrohlich/conduction/standard_feynman_data_conduction.csv")
    name_arr, chem_arr, alpha_arr, ZPR_arr, mass_arr, freq_arr, alpha_paper_arr, res_arr = looping(data)
    error = (ustrip.(ZPR_arr) - res_arr)./res_arr * 100
    column_names = ["Name", "Formula", "Mass", "Frequency [meV]", "Alpha", "Reference_Alpha", "ZPR", "Reference_ZPR", "Error"]
    df_Feynman = DataFrame([name_arr, chem_arr,  mass_arr, freq_arr, alpha_arr,
                            alpha_paper_arr, ustrip.(ZPR_arr), res_arr, error], column_names)                
    CSV.write("Data/Feynman.tsv", df_Feynman, delim='\t', quotechar='"', header=true)
    return df_Feynman
end

function Standard_Data()
    data_standard = CSV.File("LiegeDataset/Results/StandardFrohlich/conduction/standard_froelich_data_conduction")
    name_arr, chem_arr, alpha_arr, ZPR_arr, mass_arr, freq_arr, alpha_standard_arr, res_standard_arr = looping(data_standard)
    error = (ustrip.(ZPR_arr) - res_standard_arr)./res_standard_arr * 100
    column_names = ["Name", "Formula", "Mass", "Frequency [meV]", "Alpha", "Reference_Alpha", "ZPR", "Reference_ZPR", "Error"]
    df_Standard = DataFrame([name_arr, chem_arr,  mass_arr, freq_arr, alpha_arr,
                            alpha_standard_arr, ustrip.(ZPR_arr), res_standard_arr, error], column_names)                
    CSV.write("Data/Standard.tsv", df_Standard, delim='\t', quotechar='"', header=true)
    return df_Standard
end

function load_general_data()
    files = readdir("LiegeDataset/Results/GeneralizedFrohlich/conduction")
    function general_val(data)
        return data[1][1], data[1][2], data[1][4], data[1][5], data[1][6], data[1][7]
    end
    dummy_name = []
    dummy_chem = []
    dummy_mass = []
    dummy_freq = []
    dummy_alpha = []
    dummy_ZPR = []
    for i in 2:1396
        data = CSV.File("LiegeDataset/Results/GeneralizedFrohlich/conduction/" * files[i])
        name, chem, freq, mass, res_alpha, res_paper = general_val(data)
        mass = mass * mass
        push!(dummy_name, name)
        push!(dummy_chem, chem)
        append!(dummy_alpha, res_alpha)
        append!(dummy_mass, mass)
        append!(dummy_freq, freq)
        append!(dummy_ZPR, res_paper)
        end
    return dummy_name, dummy_chem, dummy_mass, dummy_freq, dummy_alpha, dummy_ZPR
end
function General_Data()
    name_arr, chem_arr, mass_arr, freq_arr, alpha_arr, ZPR_arr = load_general_data()
    column_names = ["Name", "Formula", "Mass", "Frequency [eV]", "Reference_Alpha", "Reference_ZPR"]
    df_General = DataFrame([name_arr, chem_arr,  mass_arr, freq_arr, alpha_arr,
                            ZPR_arr], column_names)                
    CSV.write("Data/General.tsv", df_General, delim='\t', quotechar='"', header=true)
    return df_General
end


function General_Comparison_Data(df_General)

    dummy_name = []
    dummy_chem = []
    dummy_result = []
    dummy_reference = []
    dummy_alpha_GFr = []
    dummy_alpha_variational = []
    for i in df_General[:, "Name"]
        for j in df_Feynman[:, "Name"]
            if i == j
                push!(dummy_name, i)
                push!(dummy_chem, df_General[df_General.Name.==i, "Formula"][1])
                append!(dummy_alpha_variational, df_Feynman[df_Feynman.Name.==i, "Alpha"])
                append!(dummy_alpha_GFr, df_General[df_General.Name.==i, "Reference_Alpha"])
                append!(dummy_result, df_Feynman[df_Feynman.Name.==i, "ZPR"])
                append!(dummy_reference, df_General[df_General.Name.==i, "Reference_ZPR"])
            end
        end
    end

    Error = (dummy_result - dummy_reference)./dummy_reference * 100
    column_names = ["Name", "Formula", "Alpha", "Reference_Alpha", "ZPR", "Reference_ZPR", "Error"]
    df_General_final = DataFrame([dummy_name, dummy_chem, dummy_alpha_variational, dummy_alpha_GFr, dummy_result, dummy_reference, Error], column_names)
    CSV.write("Data/General_comparison.tsv", df_General_final, delim='\t', quotechar='"', header=true)
    return df_General_final
end


##
df_Feynman = Feynman_Data()
df_Standard = Standard_Data()
df_General = General_Data()
df_General_final = General_Comparison_Data(df_General)



#= function plot_Feynman()
    scatter(df_Feynman[:, "Alpha"], df_Feynman[:, "Error"], legend=false, xticks = 0:1:30, yguidefontsize=10, xguidefontsize=10, titlefontsize = 10, xlims=(0, 30), ylims=(-5, 5))
    xlabel!("α")
    ylabel!("ZPR percentage difference")
    display(Plots.plot!())
    savefig("ZPR, α comparison Feynman.png")
end
plot_Feynman()

df_Standard = CSV.read("Standard.tsv", DataFrame)
scatter(df_Standard[:, "Alpha"], df_Standard[:, "Error"], legend = false, xticks = 0:1:30, yguidefontsize=10, xguidefontsize=10, titlefontsize = 10, xlims=(0, 30), ylims=(-5, 100))
xlabel!("α")
ylabel!("ZPR percentage difference")
savefig("ZPR, α comparison standard.png")
display(Plots.plot!())

#scatter(df_General_final[:, "Alpha"], df_General_final[:, "Error"], legend = false,yguidefontsize=10, xguidefontsize=10, titlefontsize = 10, xlims = [0,30], ylims = [-100, 100])
scatter(df_General_final[:, "Reference_Alpha"], df_General_final[:, "Error"], legend = false,yguidefontsize=10, xguidefontsize=10, titlefontsize = 10, xlims = [0,30], ylims = [-100, 200])
xlabel!("α (Feynman)")
ylabel!("ZPR percentage difference")
savefig("ZPR, α comparison general Feynman α.png")
display(Plots.plot!())

scatter(df_General_final[:, "Alpha"], df_General_final[:, "Error"], legend = false,yguidefontsize=10, xguidefontsize=10, titlefontsize = 10, xlims = [0,30], ylims = [-100, 200])
xlabel!("α (General)")
ylabel!("ZPR percentage difference")
savefig("ZPR, α comparison general general α.png")
display(Plots.plot!())

scatter(df_General_final[:, "Alpha"], df_General_final[:, "Reference_Alpha"], legend = false,yguidefontsize=10, xguidefontsize=10, titlefontsize = 10, xlims = [0,30], ylims = [0,30])
xlabel!("α (Feynman)")
ylabel!("α (General)")
savefig("α comparison.png")
display(Plots.plot!()) =#

