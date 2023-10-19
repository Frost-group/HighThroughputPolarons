using PolaronMobility
using Statistics
using LinearAlgebra
using CSV
using Unitful
using Gnuplot
using Plots
##
kB = 1.380649e-23
ħ = 1.054571817e-34
##
files = readdir("LiegeDataset/Results/GeneralizedFrohlich/conduction")
##
function variable_cal(data)
    return data[1][1], data[1][2], data[1][4] * 0.2417990504024 * 1000, data[1][5], data[1][6], data[1][7]
end
##
data = CSV.File("LiegeDataset/Results/GeneralizedFrohlich/conduction/" * files[2])
name, chem, freq, mass, res_alpha, res_paper = variable_cal(data)
println(name)
println(chem)
println(freq)
println(mass)
println(res_alpha)
println(res_paper)

##
p = polaron(res_alpha, ω=freq, ωeff=freq, mb = mass, β0 = ħ/kB*1e12)
addunits!(p)
ZPR = p.F0 |> u"meV"
##
function looping(files)
    dummy_alpha = []
    dummy_ZPR = []
    dummy_mass = []
    dummy_freq = []
    dummy_res = []
    dummy_alpha_paper = []
    for i in 2:1396
        data = CSV.File("LiegeDataset/Results/GeneralizedFrohlich/conduction/" * files[i])
        name, chem, freq, mass, res_alpha, res_paper = variable_cal(data)
        mass = mass * mass
        println(i)


        try
            p = polaron(res_alpha, ω=freq, ωeff=freq, mb = mass, v_guesses = 3.0001, w_guesses = 2.999, β0 = ħ/kB*1e12*2π, verbose = false)
            #p = polaron(res_alpha ,verbose = false)
            addunits!(p)
            ZPR = p.F0 |> u"meV"
            alpha = p.α
            freq = freq / 0.2417990504024 / 1000
            append!(dummy_ZPR, ZPR)
            append!(dummy_mass, mass)
            append!(dummy_freq, freq)
            append!(dummy_res, res_paper)
            append!(dummy_alpha_paper, res_alpha)
        catch e
            # Catch and handle the error
            println("An error occurred: ", e, i)
        end
        end
    return dummy_ZPR, dummy_mass, dummy_freq, dummy_alpha_paper, dummy_res
end

ZPR_arr, mass_arr, freq_arr, alpha_general_arr, res_general_arr = looping(files)
##

error_general = (res_general_arr - ustrip.(ZPR_arr))./res_general_arr * 100
##
#scatter(alpha_general_arr, error_general, legend = false, xticks = 0:1:30, yguidefontsize=10, xguidefontsize=10, titlefontsize = 10, xlims=(0, 30), ylims=(-5, 100))
scatter(alpha_general_arr, error_general, legend = false, yguidefontsize=10, xguidefontsize=10, titlefontsize = 10, xlims=(0, 30), ylims = (-10000, 200))
xlabel!("α")
ylabel!("ZPR percentage difference")
savefig("ZPR, α comparison general.png")
display(Plots.plot!())
##
dummy = []
name = "abcde"
push!(dummy, name)
print(dummy)