
function f(t::Float64, N::Int64, a_0::Float64, a::Vector{Float64})::Float64
    tmp::Float64 = a_0
    for i in 1:N
        tmp += a[i]*t^i
    end
    return tmp
end

function g(t::Float64, N::Int64, a_0::Float64, a::Vector{Float64})::Float64
    tmp::Float64 = 0.0
    for i in N:1
        tmp = tmp*t + a[i]
    end
    return tmp*t + a_0
end

N_array = [10^i for i in 5:8]
y_ticks = [10^i for i::Float64 in -8:1]
result_f = Float64[]
runtime_f = Float64[]
result_g = Float64[]
runtime_g = Float64[]
t = 0.9::Float64

for i_N in N_array
    local a_0::Float64 = rand()
    local a::Vector{Float64} = rand(Float64, i_N)

    result, runtime = @timed f(t, i_N, a_0, a)
    # print(typeof(result))
    push!(result_f, result)
    push!(runtime_f, runtime)

    result, runtime = @timed g(t, i_N, a_0, a)
    push!(result_g, result)
    push!(runtime_g, runtime)
end

print(runtime_f)
print(runtime_g)

using Plots
plot(N_array, runtime_f, scale=:log10, grid=true, xticks=N_array, yticks=y_ticks, label="N(N+1)/2")
plot!(N_array, runtime_g, scale=:log10, grid=true, xticks=N_array, yticks=y_ticks, label="N")

savefig("./plot.png")

N_len = length(N_array)

open("output.txt","w") do out #ファイルに書く
    for i in 1:N_len
        println(out, N_array[i], " ", result_f[i], " ", runtime_f[i])
    end
    println(out, "")
    for i in 1:N_len
        println(out, N_array[i], " ", result_g[i], " ", runtime_g[i])
    end
end