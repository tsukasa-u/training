using LinearAlgebra
using ForwardDiff
using Plots
using CPUTime


c(x, u) = begin
    rx, ry, ϕ = x
    u_, = u
    return [
        u_ - 1.5,
        -u_ - 1.5,
        ry - 1.0,
        -ry - 1.0,
        1.0 - ((rx + 5.0)^2 + (ry + 1.0)^2),
        0.5 - ((rx + 8.0)^2 + (ry - 0.2)^2),
        1.5 - ((rx + 2.5)^2 + (ry - 1.0)^2),
    ]
end
cx(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(x -> c(x, u), x)
cu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(u -> c(x, u), u)
cxx(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(x -> ForwardDiff.jacobian(x -> c(x, u), x), x)
cux(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(x -> ForwardDiff.jacobian(u -> c(x, u), u), x)
# fxu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(u -> ForwardDiff.jacobian(x -> c(x, u), x), u)
cuu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(u -> ForwardDiff.jacobian(u -> c(x, u), u), u)

f(x, u) = begin
    rx, ry, ϕ = x
    u_, = u
    h = 0.1
    v = 1.5
    return [
        h*v*cos(ϕ),
        h*v*sin(ϕ),
        h*u_
    ] + x
end
fx(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(x -> f(x, u), x)
fu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(u -> f(x, u), u)
fxx(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(x -> ForwardDiff.jacobian(x -> f(x, u), x), x)
fux(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(x -> ForwardDiff.jacobian(u -> f(x, u), u), x)
# fxu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(u -> ForwardDiff.jacobian(x -> f(x, u), x), u)
fuu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(u -> ForwardDiff.jacobian(u -> f(x, u), u), u)

l(x, u) = begin
    rx, ry, ϕ = x
    u_, = u
    return 0.1*(rx^2+ry^2+ϕ^2+0.1*u_^2)
end
lx(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.gradient(x -> l(x, u), x)
lu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.gradient(u -> l(x, u), u)
lxx(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.hessian(x -> l(x, u), x)
lux(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(x -> ForwardDiff.gradient(u -> l(x, u), u), x)
# lxu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(u -> ForwardDiff.gradient(x -> l(x, u), x), u)
luu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.hessian(u -> l(x, u), u)

lf(x) = begin
    rx, ry, ϕ = x
    return 0.1*(rx^2+ry^2+ϕ^2)
end

vN(x) = begin
    return lf(x)
end
vNx(x::Vector{T}) where T<:Real = ForwardDiff.gradient(x -> vN(x), x)
vNxx(x::Vector{T}) where T<:Real = ForwardDiff.hessian(x -> vN(x), x)

mutable struct struct_V{T<:Real, S<:Integer}
    nx::S
    
    ΔV::T

    Vx::Vector{T}

    Vxx::Matrix{T}

    function struct_V(_nx::S, T::DataType) where S
        self = new{T, S}()

        self.nx = _nx

        self.Vx = zeros(_nx)
        self.Vxx = zeros(_nx, _nx)

        return self
    end
end

mutable struct struct_Q{T<:Real, S<:Integer}
    ns::S
    nx::S
    nu::S

    Qx::Vector{T}
    Qu::Vector{T}
    Qs::Vector{T}

    Qsx::Matrix{T}
    Qsu::Matrix{T}
    Qss::Matrix{T}
    Qxx::Matrix{T}
    Qxu::Matrix{T}
    Quu::Matrix{T}

    function struct_Q(_nx::S, _nu::S, _ns::S, T::DataType) where S
        self = new{T, S}()

        self.nx = _nx
        self.nu = _nu
        self.ns = _ns

        self.Qx = zeros(_nx)
        self.Qu = zeros(_nu)
        self.Qs = zeros(_ns)

        self.Qsx = zeros(_ns, _nx)
        self.Qsu = zeros(_ns, _nu)
        self.Qss = zeros(_ns, _ns)
        self.Qxx = zeros(_nx, _nx)  
        self.Qxu = zeros(_nx, _nu)
        self.Quu = zeros(_nu, _nu)

        return self
    end
end

mutable struct tuple_w{T<:Real, S<:Integer}
    ns::S
    nx::S
    nu::S

    x::Vector{T}
    u::Vector{T}
    s::Vector{T}
    y::Vector{T}

    function tuple_w(_nx::S, _nu::S, _ns::S, T::DataType) where S
        self = new{T, S}()

        self.nx = _nx
        self.nu = _nu
        self.ns = _ns

        self.x = ones(_nx)
        self.u = zeros(_nu)
        self.s = ones(_ns)
        self.y = ones(_ns)

        return self
    end
end

mutable struct tuple_r{T<:Real, S<:Integer}
    ns::S

    rd::Vector{T}
    rp::Vector{T}
    rhat::Vector{T}

    function tuple_r(_ns::S, T::DataType) where S
        self = new{T, S}()

        self.ns = _ns

        self.rd = zeros(_ns)
        self.rp = zeros(_ns)
        self.rhat = zeros(_ns)

        return self
    end
end

mutable struct struct_coefficients{T<:Real, S<:Integer}
    ns::S
    nx::S
    nu::S

    α::Vector{T}
    η::Vector{T}
    χ::Vector{T}
    β::Matrix{T}
    θ::Matrix{T}
    ζ::Matrix{T}

    function struct_coefficients(_nx::S, _nu::S, _ns::S, T::DataType) where S
        self = new{T, S}()

        self.nx = _nx
        self.nu = _nu
        self.ns = _ns

        self.α = zeros(_nu)
        self.η = zeros(_ns)
        self.χ = zeros(_ns)

        self.β = zeros(_nu, _nx)
        self.θ = zeros(_ns, _nx)
        self.ζ = zeros(_ns, _nx)

        return self
    end
end

function init_V!(_V::struct_V{T, S}, w::tuple_w{T, S}) where {T<:Real, S<:Integer}

    # _V.ΔV = 0.0
    _V.Vx[:] .= vNx(w.x)
    _V.Vxx[:, :] .= vNxx(w.x)

end

function compute_Q!(_Q::struct_Q{T, S}, _V::struct_V{T, S}, w::tuple_w{T, S}) where {T<:Real, S<:Integer}

    _Q.Qs[:] .= c(w.x, w.u)
    _Q.Qsx[:, :] .= cx(w.x, w.u)
    _Q.Qsu[:, :] .= cu(w.x, w.u)
    _Q.Qss[:, :] .= zeros(_Q.ns, _Q.ns)

    _fx::Matrix{T} = fx(w.x, w.u)
    _fu::Matrix{T} = fu(w.x, w.u)
    _fx_T::Matrix{T} = transpose(_fx)
    _fu_T::Matrix{T} = transpose(_fu)
    _fxx::Array{T, 3} = reshape(fxx(w.x, w.u), (_Q.nx, _Q.nx, _Q.nx))
    _fux::Array{T, 3} = reshape(fux(w.x, w.u), (_Q.nx, _Q.nu, _Q.nx))
    _fuu::Array{T, 3} = reshape(fuu(w.x, w.u), (_Q.nx, _Q.nu, _Q.nu))
    _cxx::Array{T, 3} = reshape(cxx(w.x, w.u), (_Q.ns, _Q.nx, _Q.nx))
    _cux::Array{T, 3} = reshape(cux(w.x, w.u), (_Q.ns, _Q.nu, _Q.nx))
    _cuu::Array{T, 3} = reshape(cuu(w.x, w.u), (_Q.ns, _Q.nu, _Q.nu))

    _Q.Qx[:] .= lx(w.x, w.u) + _fx_T*_V.Vx
    _Q.Qu[:] .= lu(w.x, w.u) + _fu_T*_V.Vx

    _Q.Qxx[:, :] .=           lxx(w.x, w.u)  + _fx_T*_V.Vxx*_fx + sum(_V.Vx[i] *           _fxx[i, :, :]  for i in 1:_Q.nx)
    # _Q.Qxu[:, :] .= lxu(w.x, w.u) + _fx_T*_V.Vxx*_fu + sum(_V.Vx[i] * _fxu[i, :, :] for i in 1:_Q.nx)
    _Q.Qxu[:, :] .= transpose(lux(w.x, w.u)) + _fx_T*_V.Vxx*_fu + sum(_V.Vx[i] * transpose(_fux[i, :, :]) for i in 1:_Q.nx)
    _Q.Quu[:, :] .=           luu(w.x, w.u)  + _fu_T*_V.Vxx*_fu + sum(_V.Vx[i] *           _fuu[i, :, :]  for i in 1:_Q.nx)


    _Q.Qx[:] .+= transpose(_Q.Qsx)*w.s
    _Q.Qu[:] .+= transpose(_Q.Qsu)*w.s
    _Q.Qxx[:, :] .+= sum(w.s[i] *           _cxx[i, :, :]  for i in 1:_Q.ns)
    _Q.Qxu[:, :] .+= sum(w.s[i] * transpose(_cux[i, :, :]) for i in 1:_Q.ns)
    _Q.Quu[:, :] .+= sum(w.s[i] *           _cuu[i, :, :]  for i in 1:_Q.ns)

end

φ(x, w::tuple_w, coeff::struct_coefficients) = w.u + coeff.α + coeff.β*(x - w.x)
ψ(x, w::tuple_w, coeff::struct_coefficients) = w.s + coeff.η + coeff.θ*(x - w.x)
ξ(x, w::tuple_w, coeff::struct_coefficients) = w.y + coeff.χ + coeff.ζ*(x - w.x)

function compute_coeff!(coeff::struct_coefficients, w::tuple_w, r::tuple_r, Q::struct_Q, μ::Vector{T}) where T<:Real
    _S::Matrix{T} = diagm(w.s)
    Y::Matrix{T} = diagm(w.y)
    # I = ones(w.ns, w.ns)
    # zero_us = zeros(w.nu, w.ns)
    # zero_su = zeros(w.ns, w.nu)
    # zero_ss = zeros(w.ns, w.ns)
    # zero_sx = zeros(w.ns, w.nx)

    r.rp[:] .= c(w.x, w.u) + w.y
    r.rd[:] .= _S*w.y - μ
    r.rhat[:] .= _S*r.rp - r.rd

    # P = [Q.Quu transpose(Q.Qsu) zero_us; Q.Qsu zero_ss I; zero_su Y S]

    # ret = -inv(P)*[Q.Qu Q.Qux; rp Q.Qsx; rd zero_sx]

    # coeff.α .= vec(ret[1:w.nx, 1])
    # coeff.η .= vec(ret[(w.nx + 1):(w.nx + w.nu), 1])
    # coeff.χ .= vec(ret[(w.nx + w.nu + 1):end, 1])

    # coeff.β .= ret[1:w.nx, 2:end]
    # coeff.θ .= ret[(w.nx + 1):(w.nx + w.nu), 2:end]
    # coeff.ζ .= ret[(w.nx + w.nu + 1):end, 2:end]

    Y_1::Matrix{T} = inv(Y)
    SY_1::Matrix{T} = _S*Y_1

    R::Matrix{T} = (Q.Quu + transpose(Q.Qsu)*SY_1*Q.Qsu)
    ret::Matrix{T} = -inv(R)*[Q.Qu + transpose(Q.Qsu)*Y_1*r.rhat transpose(Q.Qxu) + transpose(Q.Qsu)*SY_1*Q.Qsx]
    
    coeff.α[:] .= vec(ret[:, 1])

    coeff.β[:, :] .= ret[:, 2:end]

    coeff.η[:] .= Y_1*(r.rhat + _S*Q.Qsu*coeff.α)
    coeff.θ[:, :] .= SY_1*(Q.Qsx + Q.Qsu*coeff.β)

    coeff.χ[:] .= -r.rp - Q.Qsu*coeff.α
    coeff.ζ[:, :] .= -Q.Qsx - Q.Qsu*coeff.β

end

function update_V!(Q::struct_Q{T, S}, V::struct_V{T, S},  coeff::struct_coefficients{T, S}) where {T<:Real, S<:Integer}

    # _V.ΔV += 0.0
    _β_T = transpose(coeff.β)
    V.Vx[:] .= Q.Qx + _β_T*Q.Qu + _β_T*Q.Quu*coeff.α + Q.Qxu*coeff.α
    V.Vxx[:, :] .= Q.Qxx + Q.Qxu*coeff.β + _β_T*transpose(Q.Qxu) + _β_T*Q.Quu*coeff.β
end

function update_Q!(Q::struct_Q{T, S}, w::tuple_w{T, S}, r::tuple_r{T, S}) where {T<:Real, S<:Integer}
    Y_1::Matrix{T} = inv(diagm(w.y))
    _S::Matrix{T} = diagm(w.s)
    SY_1::Matrix{T} = _S*Y_1

    Qsx_T::Matrix{T} = transpose(Q.Qsx)
    Qsu_T::Matrix{T} = transpose(Q.Qsu)

    Q.Qx[:] .+= Qsx_T*Y_1*r.rhat
    Q.Qu[:] .+= Qsu_T*Y_1*r.rhat
    Q.Qxx[:, :] .+= Qsx_T*SY_1*Q.Qsx
    Q.Qxu[:, :] .+= Qsx_T*SY_1*Q.Qsu
    Q.Quu[:, :] .+= Qsu_T*SY_1*Q.Qsu

    # Q.Qs[:] .= c(w.x, w.u)
    # Q.Qsx[:, :] .= cx(w.x, w.u)
    # Q.Qsu[:, :] .= cu(w.x, w.u)
    # Q.Qss[:, :] .= zeros(Q.ns, Q.ns)
end

function FFP!(nw::S, list_w::Array{tuple_w{T, S}, 1}, list_coeff::Array{struct_coefficients{T, S}, 1}) where {S<:Integer, T<:Real}
    x::Vector{T} = copy(list_w[1].x)
    u::Vector{T} = zeros(list_w[1].nu)
    s::Vector{T} = zeros(list_w[1].ns)
    y::Vector{T} = zeros(list_w[1].ns)
    for i in 1:nw
        u[:] .= φ(x, list_w[i], list_coeff[i])
        s[:] .= ψ(x, list_w[i], list_coeff[i])
        y[:] .= ξ(x, list_w[i], list_coeff[i])

        list_w[i].u[:] .= u
        list_w[i].s[:] .= s
        list_w[i].y[:] .= y
        list_w[i].x[:] .= x

        # check_constraints(list_w[i])

        x[:] .= f(x, u)
    end
    
end

function BFP!(nw::S, list_w::Array{tuple_w{T, S}, 1}, list_r::Array{tuple_r{T, S}, 1}, list_coeff::Array{struct_coefficients{T, S}, 1}, μ::Vector{T}) where {S<:Integer, T<:Real}
    
    Q::struct_Q = struct_Q(list_w[1].nx, list_w[1].nu, list_w[1].ns, T)
    V::struct_V = struct_V(list_w[1].nx, T)
    init_V!(V, list_w[nw])
    for i in reverse(1:nw)
        compute_Q!(Q, V, list_w[i])
        update_Q!(Q, list_w[i], list_r[i])
        compute_coeff!(list_coeff[i], list_w[i], list_r[i], Q, μ)
        update_V!(Q, V, list_coeff[i])
    end
end

function check_constraints(w::tuple_w{T, S}) where {S<:Integer, T<:Real}
    w.s[:] .= max.(1E-3::T, w.s)
    w.y[:] .= max.(1E-3::T, w.y)
    w.u[:] .= max.(0.0::T, min.(0.01, w.u))
    w.x[:] .= max.(1E-3::T, w.x)
end

function loop!(n::S, nw::S, list_w::Array{tuple_w{T, S}, 1}, list_r::Array{tuple_r{T, S}, 1}, list_coeff::Array{struct_coefficients{T, S}, 1}, μ::Vector{T}) where {S<:Integer, T<:Real}
    
    wrap_plot_graph(0, list_w, nw)
    for k in 1:n
        BFP!(nw, list_w, list_r, list_coeff, μ)
        FFP!(nw, list_w, list_coeff)
        update_μ!(μ, 5.0)
        wrap_plot_graph(k, list_w, nw)
    end
end

function wrap_plot_graph(idx, list_w::Array{tuple_w{T, S}, 1}, nw::S) where {S<:Integer, T<:Real}
    plot_graph(
        idx, 
        [list_w[i].x[j] for j in 1:list_w[1].nx, i in 1:nw],
        [list_w[i].u[j] for j in 1:list_w[1].nu, i in 1:nw], 
        [list_w[i].s[j] for j in 1:list_w[1].ns, i in 1:nw], 
        [list_w[i].y[j] for j in 1:list_w[1].ns, i in 1:nw], 
        [0.01 * i for i in 1:nw]
    )
end

function init_μ!(μ::Vector{T}, list_w::Array{tuple_w{T, S}, 1}, nw::S, ns::S) where {S<:Integer, T<:Real}
    μ[:] .= (sum(l(list_w[i].x, list_w[i].u) for i in 1:nw)/nw/ns + lf(list_w[end].x))*zeros(ns)
end

function update_μ!(μ::Vector{T}, κ::T) where T<:Real
    μ[:] .= min.(μ/κ, μ.^1.2)
end

function init_xu!(list_w::Array{tuple_w{T, S}, 1}, nw::S, bx::Tuple{Array{T, 1}, Array{T, 1}}, bu::Tuple{Array{T, 1}, Array{T, 1}}) where {S<:Integer, T<:Real}
    bx0, bxf = bx
    bx0f = bxf .- bx0
    bu0, buf = bu
    bu0f = buf .- bu0
    for i in 1:nw
        list_w[i].x .= bx0f*i/nw .+ bx0
        list_w[i].u .= bu0f*i/nw .+ bu0
    end
end

function main()
    nw::Int64 = 600
    nx::Int64 = 3
    nu::Int64 = 1
    ns::Int64 = 7
    list_w::Array{tuple_w{Float64, Int64}, 1} = [tuple_w(nx, nu, ns, Float64) for _ in 1:nw]
    list_r::Array{tuple_r{Float64, Int64}, 1} = [tuple_r(ns, Float64) for _ in 1:nw]
    list_coeff::Array{struct_coefficients{Float64, Int64}, 1} = [struct_coefficients(nx, nu, ns, Float64) for _ in 1:nw]
    μ::Vector{Float64} = zeros(ns)

    init_xu!(list_w, nw, ([-10.0, 0.0, 0.0], [-10.0, 0.0, 0.0]), ([0.0], [0.0]))
    init_μ!(μ, list_w, nw, ns)

    loop!(10, nw, list_w, list_r, list_coeff, μ)
end


function plot_graph(index, plot_x, plot_u, plot_s, plot_y, plot_t)

    plots_x = Plots.plot(
        plot_t[:], 
        [ plot_x[1, :] plot_x[2, :]], 
        label=["rx" "ry" ],
        xlabel = "t",
        st=:scatter,
    )

    plots_u = Plots.plot(
        plot_t, 
        [ plot_u[1, :] ], 
        label=["u"],
        xlabel = "t",
        st=:scatter
    )

    plots_s = Plots.plot(
        plot_t, 
        [ plot_s[i, :] for i in 1:7], 
        xlabel = "t",
        st=:scatter
    )

    plots_y = Plots.plot(
        plot_t, 
        [ plot_y[i, :] for i in 1:7], 
        xlabel = "t",
        st=:scatter
    )
    
    plot_ref = Plots.plot(
        plots_x,
        plots_u,
        plots_s,
        plots_y,
        layout = (2, 2),
        # legend = false,
        # margin = 1Plots.cm,
    )
    png(string(index, base = 10, pad = 2))
end


@time @CPUtime main()