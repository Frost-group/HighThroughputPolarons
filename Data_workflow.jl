using PolaronMobility
using Statistics
using LinearAlgebra
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

function variable_cal(data)
    """
    variable_cal(data)
    The `variable_cal` function takes in a data array as input and returns a modified version of the data array.

    # Arguments
    - `data` (array): An array containing data elements.

    # Returns
    A tuple containing the modified elements of the input data array.
    """

    return data[1], data[2], data[3], data[4], data[5] * 0.2417990504024, data[6], data[7], data[8] # constant for unit conversion form meV to THz
end

function looping(data)
    """
    looping(data)

    The `looping` function takes in an array of data and iterates over each element.
    For each element, it calls the `variable_cal` function to calculate some values.
    It then performs some calculations and creates a material object.
    If the calculations are successful, it adds the calculated values to separate arrays.
    If an error occurs during the calculations, it prints an error message.
    Finally, it returns the arrays containing the calculated values.

    # Arguments
    - `data` (array): An array containing data elements.

    # Returns
    - `dummy_name` (array): An array containing the names of the elements in the `data` array.
    - `dummy_chem` (array): An array containing the chemical properties of the elements in the `data` array.
    - `dummy_alpha` (array): An array containing the calculated alpha values.
    - `dummy_ZPR` (array): An array containing the calculated ZPR values.
    - `dummy_mass` (array): An array containing the calculated mass values.
    - `dummy_freq` (array): An array containing the calculated frequency values.
    - `dummy_alpha_paper` (array): An array containing the calculated alpha paper values.
    - `dummy_res` (array): An array containing the calculated res_paper values.
    """

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
            name = replace(name, " " => "")
            chem = replace(chem, " " => "")
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
    """
    Feynman_Data() -> DataFrame

    The `Feynman_Data` function reads data from Feynman Frohlich Model data file and passes it to the `looping` function.
    It then performs calculations on the data and creates a DataFrame.
    Finally, it writes the DataFrame to a TSV file and returns it.
    """

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
    """
    Standard_Data() -> DataFrame

    The `Standard_Data` function reads data from Standard Frohlich Model data file and passes it to the `looping` function.
    It then performs calculations on the data and creates a DataFrame.
    Finally, it writes the DataFrame to a TSV file and returns it.
    """

    data_standard = CSV.File("LiegeDataset/Results/StandardFrohlich/conduction/standard_froelich_data_conduction.csv")
    name_arr, chem_arr, alpha_arr, ZPR_arr, mass_arr, freq_arr, alpha_standard_arr, res_standard_arr = looping(data_standard)
    error = (ustrip.(ZPR_arr) - res_standard_arr)./res_standard_arr * 100
    column_names = ["Name", "Formula", "Mass", "Frequency [meV]", "Alpha", "Reference_Alpha", "ZPR", "Reference_ZPR", "Error"]
    df_Standard = DataFrame([name_arr, chem_arr,  mass_arr, freq_arr, alpha_arr,
                            alpha_standard_arr, ustrip.(ZPR_arr), res_standard_arr, error], column_names)                
    CSV.write("Data/Standard.tsv", df_Standard, delim='\t', quotechar='"', header=true)
    return df_Standard
end

function General_Data()
    """
    General_Data() -> DataFrame

    The `General_Data` function reads data from General Frohlich Model data file and passes it to the `looping` function.
    It then sorts the data and creates a DataFrame.
    Finally, it writes the DataFrame to a TSV file and returns it.
    """

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
        name = replace(name, " " => "")
        chem = replace(chem, " " => "")
        push!(dummy_name, name)
        push!(dummy_chem, chem)
        append!(dummy_alpha, res_alpha)
        append!(dummy_mass, mass)
        append!(dummy_freq, freq)
        append!(dummy_ZPR, res_paper)
        end
    column_names = ["Name", "Formula", "Mass", "Frequency [meV]", "Reference_Alpha", "Reference_ZPR"]
    df_General = DataFrame([dummy_name, dummy_chem, dummy_mass, dummy_freq, dummy_alpha, dummy_ZPR], column_names)                
    CSV.write("Data/General.tsv", df_General, delim='\t', quotechar='"', header=true)
    return df_General
end

function General_Comparison_Data(df_General)
    """
    General_Comparison_Data() -> DataFrame

    This function compares the data of GFM and Variational Approach.
    Then it writes the DataFrame to a TSV file and returns it.
    """

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

function multi_mode()
    """
    multi_mode() -> DataFrame

    The `multi_mode` function reads data from General Frohlich Model data file and passes it to the multi-mode polaron calculation (similar to "looping" function).
    It writes the DataFrame to a TSV file and returns it.
    """
    ħ = 1.054571817e-34
    kB = 1.380649e-23
    files = readdir("LiegeDataset/Results/GeneralizedFrohlich/conduction")
    function general_val(data)
        return data[:, 1], data[:, 3], data[:, 4], data[:, 5], data[:, 6]
    end
    dummy_name = []
    dummy_freq = []
    dummy_alpha = []
    dummy_ZPR = []
    dummy_result = []
    for i in 1396:length(files)
    #mp-10086
    #for i in 2084
        data = CSV.File("LiegeDataset/Results/GeneralizedFrohlich/conduction/" * files[i], delim = "\t") |> DataFrame
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
        println(i, files[i], freq_actual, res_alpha, res_paper)
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
        parts = split(files[i], '-')
        name = join(parts[1:2], '-')
        push!(dummy_name, name)
        append!(dummy_alpha, sum(res_alpha))
        append!(dummy_freq, (sum(freq)*(4 * pi * area[1]))^(-2))
        append!(dummy_result, ZPR)
        append!(dummy_ZPR, sum(res_paper))

    end
    
    column_names = ["Name", "Frequency [meV]", "ZPR [meV]", "Reference_Alpha", "Reference_ZPR [meV]"]
    df_General = DataFrame([[name], [round(sum(freq_actual), digits=6)], [round(sum(ZPR), digits=6)], [round(sum(res_alpha), digits=6)], [round(sum(res_paper), digits=6)]], column_names)         
    CSV.write("Data/multi_mode.tsv", df_General, delim='\t', quotechar='"', header=true)
    return df_General
end

function general_val_multi(data)
    return data[1][1], data[1][2], data[1][3], data[1][4], data[1][5]
end

function multi_mode_sorting(mode = 0)
    if mode == 0

        files = readdir("Data/multi_mode")

    elseif mode == 1

        files = readdir("Data/multi_mode_zero_K")
    end

    dummy_name = []
    dummy_freq = []
    dummy_ZPR = []
    dummy_alpha = []
    dummy_ZPR_paper = []
    for i in 1:length(files)
        
        if mode == 0

            data = CSV.File("Data/multi_mode/" * files[i])

        elseif mode == 1

            data = CSV.File("Data/multi_mode_zero_K/" * files[i])
        end
    name, freq, ZPR, alpha, ZPR_paper = general_val_multi(data)
    push!(dummy_name, name)
    append!(dummy_ZPR, ZPR)
    append!(dummy_alpha, alpha)
    append!(dummy_freq, freq)
    append!(dummy_ZPR_paper, ZPR_paper)
    end

    column_names = ["Name", "Frequency [meV]", "ZPR [meV]", "Reference_Alpha", "Reference_ZPR [meV]"]
    df_multi_mode = DataFrame([dummy_name, dummy_freq, dummy_ZPR, dummy_alpha, dummy_ZPR_paper], column_names)      
    if mode == 0

        CSV.write("Data/multi_mode.tsv", df_multi_mode, delim='\t', quotechar='"', header=true)

    elseif mode == 1

        CSV.write("Data/multi_mode_zero_K.tsv", df_multi_mode, delim='\t', quotechar='"', header=true)
    end   

    return df_multi_mode
end

function multi_single_comparison()
    data_multi = CSV.File("Data/multi_mode.tsv")
    data_single = CSV.File("Data/Standard.tsv")

    dummy_name = []
    dummy_ZPR_single = []
    dummy_alpha_single = []
    dummy_ZPR_multi = []
    dummy_alpha_multi = []
    dummy_error = []
    for i in data_multi
        for j in data_single
            if i[1] == j[1]

                name, ZPR_multi, ZPR_single, α_multi, α_single = i[1], i[3], j[7], i[4], j[5]
                error = (ZPR_multi - ZPR_single) / ZPR_single * 100
                push!(dummy_name, name)
                append!(dummy_ZPR_single, ZPR_single)
                append!(dummy_ZPR_multi, ZPR_multi)
                append!(dummy_alpha_multi, α_multi)
                append!(dummy_alpha_single, α_single)
                append!(dummy_error, error)
            end
        end
    end
    column_names = ["Name", "ZPR_multi", "ZPR_single", "α_multi", "α_single", "Error"]
    df = DataFrame([dummy_name, dummy_ZPR_single, dummy_ZPR_multi, dummy_alpha_single, dummy_alpha_multi, dummy_error], column_names)
    CSV.write("Data/single_multi_comparison.tsv", df, delim='\t', quotechar='"', header=true)
    return df
end
# Commands for saving data

#df_multimode = multi_mode() # Too long to calculate! Use multiprocessing unit to do so.