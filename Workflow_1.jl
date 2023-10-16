using PolaronMobility
using Statistics
using LinearAlgebra
using CSV
using Unitful
using Gnuplot
##
x = range(-2pi, stop=2pi, length=1000);
y = sin.(x)
name = "\$MyDataSet1"
@gp name=>(x, y) "plot $name w l lc rgb 'black'" "pl $name u 1:(1.5*\$2) w l lc rgb 'red'"      