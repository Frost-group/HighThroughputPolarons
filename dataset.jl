using PolaronMobility

MAPI = material(4.9, 24.1, 0.12, 2.25)

p = polaron(MAPI, verbose = true)

p_2 = polaron(MAPI, [1, 100, 200], [1,2,3,4], verbose = true)

##
using JSON

atom = JSON.parsefile("LiegeDataset/Repository/atoms/atomic_data.json")

phonon_files = readdir("LiegeDataset/Repository/phonon")

function load_file(dirs, files)
dummy = []
for i in files
data = read(dirs * i, String)
push!(dummy, JSON.parse(data))
end
return dummy
end


phonon = load_file("LiegeDataset/Repository/phonon/", phonon_files)

for (key, value) in pairs(phonon[1])
println("$key: $value\n")
end

mass_files = readdir("LiegeDataset/Repository/eff_masses")
mass = load_file("LiegeDataset/Repository/eff_masses/", mass_files)