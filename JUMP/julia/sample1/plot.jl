using Plots
using CSV
using DataFrames

data = CSV.read("./file.txt", DataFrame)

n = data.r
t = data.t
r = data.r
θ = data.θ
u = data.u
v = data.v
γ = data.γ

plot(t,r)
plot!(t, θ)
plot!(t, u)
plot!(t, v)
plot!(t, γ)