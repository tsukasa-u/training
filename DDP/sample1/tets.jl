
# using ForwardDiff

# l(x, u) = begin
#     r, v_r, v_t = x
#     u_rp, u_rm, u_tp, u_tm = u
#     return u_rp + u_rm + u_tp + u_tm
# end
# lx(x, u) = ForwardDiff.gradient(x -> l(x, u), x)

# function main()
#     x = [Float64[0.0, 0.0, 0.0] for _ in 1:5]
#     u = [Real[0.0, 0.0, 0.0, 0.0] for _ in 1:5]
#     for i in 1:5
#         println(x[i])
#         println(typeof(lx(x[i], u[i])))
#     end
# end
using LinearAlgebra
# main()

# cat([[1, 3], [1, 4]]...; dims=1)
# x = [1.0, 2.0, 3.0]
# ([1.0, 1.0, 1.0] .- [4.0, 0.0, 0.5]./x).^2.0

# C = [2.0323617306118305 -1.0323617306118305 0.2518891687657431 -0.2518891687657431; -1.0323617306118305 2.0323617306118305 -0.2518891687657431 0.2518891687657431; 0.2518891687657431 -0.2518891687657431 0.13024277325133896 0.8697572267486611; -0.2518891687657431 0.2518891687657431 0.8697572267486611 0.13024277325133896]
# det(C)

# a = [[2,2, 2]]
# b = [[3,5, 6]]
# abs.((a-b)...)
# findmin(abs.(a[i] - b[i]) for i in 1:1)

# minimum(minimum(abs.(a[i] - b[i]) for i in 1:1))

# a = [5, 2, -1]
# b = [[4, 4, 4], [0, 0, 0]]
# min.(max.(a, b[2]), b[1])
# a = [4, 0, -0.5]
# sum((a .- [4.0, 0.0, 0.5]).^2.0)

f(x, u) = begin
    r, v_r, v_t = x
    u_rp, u_rm, u_tp, u_tm = u
    return [
        (v_r),
        (v_t^2.0/r - 1.0/r^2.0 + u_rp - u_rm),
        # (-v_r*v_t/r + u_tp - u_tm)
    ]
end
fx(x::Vector{Real}, u::Vector{Real}) = ForwardDiff.jacobian(x -> f(x, u), x)
fu(x::Vector{Real}, u::Vector{Real}) = ForwardDiff.jacobian(u -> f(x, u), u)
fxx(x, u) = ForwardDiff.jacobian(x -> ForwardDiff.jacobian(x -> f(x, u), x), x)
fux(x, u) = ForwardDiff.jacobian(x -> ForwardDiff.jacobian(u -> f(x, u), u), x)
fuu(x, u) = ForwardDiff.jacobian(u -> ForwardDiff.jacobian(u -> f(x, u), u), u)
fxx([1.0, 2.0, 3.0], [4.0, 0.0, 0.5, 0.0])
fux([1.0, 2.0, 3.0], [4.0, 0.0, 0.5, 0.0])
fuu([1.0, 2.0, 3.0], [4.0, 0.0, 0.5, 0.0])

A  = reshape(fxx([1.0, 2.0, 3.0], [4.0, 0.0, 0.5, 0.0]), (2, 3, 3))
A[1, :, :]
A[2, :, :]