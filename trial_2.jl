using PolaronMobility
using Statistics
using LinearAlgebra
using CSV
using Unitful
using Gnuplot
using Plots
##


function load_file(file)
    open(file) do f
        data = read(f)
        return data
    end
    
end

#data_file = load_file("LiegeDataset/Results/StandardFeynmanFrohlich/conduction/standard_feynman_data_conduction")
data = CSV.File("LiegeDataset/Results/StandardFeynmanFrohlich/conduction/standard_feynman_data_conduction.csv")
##
function variable_cal(data)
    return data[1], data[2], data[3], data[4], data[5] * 0.2417990504024, data[6], data[7], data[8]
end
##
name, chem, eps_s, eps_o, freq, mass, res_alpha, res_paper = variable_cal(data[864]) #56, 864
println(variable_cal(data[864]))
#mass = mass * mass
MAPI = material(eps_o, eps_s, mass, freq)
p = polaron(MAPI, verbose = false)
println(data[864][7])

##
addunits!(p)
p.α

##
function looping(data)
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
    return dummy_alpha, dummy_ZPR, dummy_mass, dummy_freq, dummy_alpha_paper, dummy_res
end
##
alpha_arr, ZPR_arr, mass_arr, freq_arr, alpha_paper_arr, res_arr = looping(data)

##
scatter(alpha_arr, ustrip.(ZPR_arr), yguidefontsize=10, xguidefontsize=10, titlefontsize = 10, xlims=(0, 30), ylims=(-3000,3000))
xlabel!("α")
ylabel!("Ground State Energy")
title!("α Versis Ground State Energy", fontsize = 8)
display(Plots.plot!())
#savefig("alpha_ZPR_2.png")

##
scatter(freq_arr, ustrip.(ZPR_arr), yguidefontsize=10, xguidefontsize=10, titlefontsize = 10, xlims=(0, 200), ylims=(-3000,3000))
xlabel!("ωeff")
ylabel!("Ground State Energy")
title!("ωeff Versis Ground State Energy", fontsize = 8)
display(Plots.plot!())
savefig("freq_ZPR.png")

##
scatter(mass_arr, ustrip.(ZPR_arr), yguidefontsize=10, xguidefontsize=10, titlefontsize = 10, xlims=(0, 30), ylims=(-3000,3000))
xlabel!("m⋆")
ylabel!("Ground State Energy")
title!("m⋆ Versis Ground State Energy", fontsize = 8)
display(Plots.plot!())
savefig("MASS_ZPR_2.png")
#@gp alpha_arr ZPR_arr
##

scatter(alpha_arr, alpha_paper_arr, yguidefontsize=10, xguidefontsize=10, titlefontsize = 10, xlims=(0, 30), ylims=(0,30))
xlabel!("α_ours")
ylabel!("α_paper")
display(Plots.plot!())
savefig("alpha_comparison.png")
##
scatter(ustrip.(ZPR_arr), res_arr, yguidefontsize=10, xguidefontsize=10, titlefontsize = 10, xlims=(-3000, 0), ylims=(-3000,0))
xlabel!("ZPR_ours")
ylabel!("ZPR_paper")
display(Plots.plot!())
savefig("ZPR_comparison.png")

##
error = (ustrip.(ZPR_arr) - res_arr)./res_arr
##
scatter(alpha_arr, error, xticks = 0:1:30, yguidefontsize=10, xguidefontsize=10, titlefontsize = 10, xlims=(0, 30), ylims=(-0.05, 0.05))
xlabel!("α")
ylabel!("ZPR percentage difference")
display(Plots.plot!())
savefig("ZPR, α comparison.png")
##
println(length(ZPR_arr))

##
data_standard = CSV.File("LiegeDataset/Results/StandardFrohlich/conduction/standard_froelich_data_conduction")
alpha_arr, ZPR_arr, mass_arr, freq_arr, alpha_standard_arr, res_standard_arr = looping(data_standard)
##
error_standard = (ustrip.(ZPR_arr) - res_standard_arr)./res_standard_arr
scatter(alpha_arr, error_standard, xticks = 0:1:30, yguidefontsize=10, xguidefontsize=10, titlefontsize = 10, xlims=(0, 30), ylims=(-0.05, 1))
xlabel!("α")
ylabel!("ZPR percentage difference")
savefig("ZPR, α comparison standard.png")
display(Plots.plot!())
##
alpha_arr = convert(Array{Float64}, alpha_arr)
##
@gp alpha_arr error_standard :-
@gp :- title = "Fit function to values stored in an array"
save(term="pngcairo size 1100,700 fontscale 1.6", output="output.png")
