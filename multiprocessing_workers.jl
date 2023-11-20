using Distributed
#@everywhere cd("F:/OneDrive - Imperial College London/My Stuff/Imperial College/Year 5/Polaron/HighThroughputPolarons")
@everywhere include("multiprocessing_functions.jl")
num_processes = 10

# Start local worker processes
addprocs(num_processes)


@time @sync for file in files[1397: end]
    trial_with_vw(file)
end
rmprocs(workers())

