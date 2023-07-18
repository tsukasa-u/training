using LinearAlgebra
using ForwardDiff
using Plots
using CPUTime


c(x, u) = begin
    r, v_r, v_t = x
    u_rp, u_rm, u_tp, u_tm = u
    return [
        -r,
        -v_r,
        -v_t,
        u_rp*(u_rp - 0.01),
        u_rm*(u_rm - 0.01),
        u_tp*(u_tp - 0.01),
        u_tm*(u_tm - 0.01),
    ]
end
cx(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(x -> c(x, u), x)
cu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(u -> c(x, u), u)

f(x, u) = begin
    r, v_r, v_t = x
    u_rp, u_rm, u_tp, u_tm = u
    return [
        (v_r),
        (v_t^2.0/r - 1.0/r^2.0 + u_rp - u_rm),
        (-v_r*v_t/r + u_tp - u_tm)
    ]*0.01 + x
end
fx(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(x -> f(x, u), x)
fu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(u -> f(x, u), u)
fxx(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(x -> ForwardDiff.jacobian(x -> f(x, u), x), x)
fux(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(x -> ForwardDiff.jacobian(u -> f(x, u), u), x)
# fxu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(u -> ForwardDiff.jacobian(x -> f(x, u), x), u)
fuu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(u -> ForwardDiff.jacobian(u -> f(x, u), u), u)

l(x, u) = begin
    r, v_r, v_t = x
    u_rp, u_rm, u_tp, u_tm = u
    return u_rp + u_rm + u_tp + u_tm + (u_rp*u_rm + u_tp*u_tm)
end
lx(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.gradient(x -> l(x, u), x)
lu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.gradient(u -> l(x, u), u)
lxx(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.hessian(x -> l(x, u), x)
lux(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(x -> ForwardDiff.gradient(u -> l(x, u), u), x)
# lxu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(u -> ForwardDiff.gradient(x -> l(x, u), x), u)
luu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.hessian(u -> l(x, u), u)

lf(x) = sum((x .- [4.0, 0.0, 0.5]).^2.0)

vN(x) = begin
    return lf(x)
end
vNx(x::Vector{T}) where T<:Real = ForwardDiff.gradient(x -> vN(x), x)
vNxx(x::Vector{T}) where T<:Real = ForwardDiff.hessian(x -> vN(x), x)

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

function init_Q!(_Q::struct_Q, w::tuple_w)

    _Q.Qs[:] .= c(w.x, w.u)
    _Q.Qsx[:, :] .= cx(w.x, w.u)
    _Q.Qsu[:, :] .= cu(w.x, w.u)
    _Q.Qss[:, :] .= zeros(_Q.ns, _Q.ns)

    _fx = fx(w.x, w.u)
    _fu = fu(w.x, w.u)
    fx_T = transpose(_fx)
    fu_T = transpose(_fu)
    _fxx = reshape(fxx(w.x, w.u), (_Q.nx, _Q.nx, _Q.nx))
    _fux = reshape(fux(w.x, w.u), (_Q.nx, _Q.nu, _Q.nx))
    _fuu = reshape(fuu(w.x, w.u), (_Q.nx, _Q.nu, _Q.nu))

    _vNx = vNx(w.x)
    _vNxx = vNxx(w.x)

    _Q.Qx[:] .= lx(w.x, w.u) + fx_T*_vNx
    _Q.Qu[:] .= lu(w.x, w.u) + fu_T*_vNx

    _Q.Qxx[:, :] .= lxx(w.x, w.u) + fx_T*_vNxx*_fx + sum(_vNx[i] * _fxx[i, :, :] for i in 1:_Q.nx)
    # _Q.Qxu[:, :] .= lxu(w.x, w.u) + fx_T*_vNxx*_fu + sum(_vNx[i] * _fxu[i, :, :] for i in 1:_Q.nx)
    _Q.Qxu[:, :] .= transpose(lux(w.x, w.u)) + fx_T*_vNxx*_fu + sum(_vNx[i] * transpose(_fux[i, :, :]) for i in 1:_Q.nx)
    _Q.Quu[:, :] .= luu(w.x, w.u) + fu_T*_vNxx*_fu + sum(_vNx[i] * _fuu[i, :, :] for i in 1:_Q.nx)
end

φ(x, w::tuple_w, coeff::struct_coefficients) = w.u + coeff.α + coeff.β*(x - w.x)
ψ(x, w::tuple_w, coeff::struct_coefficients) = w.s + coeff.η + coeff.θ*(x - w.x)
ξ(x, w::tuple_w, coeff::struct_coefficients) = w.y + coeff.χ + coeff.ζ*(x - w.x)

function compute_coeff!(coeff::struct_coefficients, w::tuple_w, r::tuple_r, Q::struct_Q, μ::Vector{T}) where T<:Real
    S = diagm(w.s)
    Y = diagm(w.y)
    I = ones(w.ns, w.ns)
    zero_us = zeros(w.nu, w.ns)
    zero_su = zeros(w.ns, w.nu)
    zero_ss = zeros(w.ns, w.ns)
    zero_sx = zeros(w.ns, w.nx)

    r.rp[:] .= c(w.x, w.u) + w.y
    r.rd[:] .= S*w.y - μ
    r.rhat[:] .= S*r.rp - r.rd

    # P = [Q.Quu transpose(Q.Qsu) zero_us; Q.Qsu zero_ss I; zero_su Y S]

    # ret = -inv(P)*[Q.Qu Q.Qux; rp Q.Qsx; rd zero_sx]

    # coeff.α .= vec(ret[1:w.nx, 1])
    # coeff.η .= vec(ret[(w.nx + 1):(w.nx + w.nu), 1])
    # coeff.χ .= vec(ret[(w.nx + w.nu + 1):end, 1])

    # coeff.β .= ret[1:w.nx, 2:end]
    # coeff.θ .= ret[(w.nx + 1):(w.nx + w.nu), 2:end]
    # coeff.ζ .= ret[(w.nx + w.nu + 1):end, 2:end]

    SY_1 = S*inv(Y)
    Y_1 = inv(Y)

    R = (Q.Quu + transpose(Q.Qsu)*SY_1*Q.Qsu)
    ret = -inv(R)*[Q.Qu + transpose(Q.Qsu)*Y_1*r.rhat transpose(Q.Qxu) + transpose(Q.Qsu)*SY_1*Q.Qsx]
    
    coeff.α[:] .= vec(ret[:, 1])

    coeff.β[:, :] .= ret[:, 2:end]

    coeff.η[:] .= Y_1*(r.rhat + S*Q.Qsu*coeff.α)
    coeff.θ[:, :] .= SY_1*(Q.Qsx + Q.Qsu*coeff.β)

    coeff.χ[:] .= -r.rp - Q.Qsu*coeff.α
    coeff.ζ[:, :] .= -Q.Qsx - Q.Qsu*coeff.β

end

function update_Q!(Q::struct_Q, w::tuple_w, r::tuple_r)
    Y_1 = inv(diagm(w.y))
    S = diagm(w.s)

    Qsx_T = transpose(Q.Qsx)
    Qsu_T = transpose(Q.Qsu)

    Q.Qx[:] .+= Qsx_T*Y_1*r.rhat
    Q.Qu[:] .+= Qsu_T*Y_1*r.rhat
    Q.Qxx[:, :] .+= Qsx_T*Y_1*S*Q.Qsx
    Q.Qxu[:, :] .+= Qsx_T*Y_1*S*Q.Qsu
    Q.Quu[:, :] .+= Qsu_T*Y_1*S*Q.Qsu

    Q.Qs[:] .= c(w.x, w.u)
    Q.Qsx[:, :] .= cx(w.x, w.u)
    Q.Qsu[:, :] .= cu(w.x, w.u)
    Q.Qss[:, :] .= zeros(Q.ns, Q.ns)
end

function FFP!(nw::S, list_w::Array{tuple_w{T, S}, 1}, list_coeff::Array{struct_coefficients{T, S}, 1}) where {S<:Integer, T<:Real}
    x = list_w[1].x
    u = zeros(list_w[1].nu)
    s = zeros(list_w[1].ns)
    y = zeros(list_w[1].ns)
    for i in 1:nw
        u .= φ(x, list_w[i], list_coeff[i])
        s .= ψ(x, list_w[i], list_coeff[i])
        y .= ξ(x, list_w[i], list_coeff[i])

        list_w[i].u[:] .= u
        list_w[i].s[:] .= s
        list_w[i].y[:] .= y
        list_w[i].x[:] .= x

        x .= f(x, u)
    end
    
end

function BFP!(nw::S, list_w::Array{tuple_w{T, S}, 1}, list_r::Array{tuple_r{T, S}, 1}, list_coeff::Array{struct_coefficients{T, S}, 1}, μ::Vector{T}) where {S<:Integer, T<:Real}
    
    Q = struct_Q(list_w[1].nx, list_w[1].nu, list_w[1].ns, Float64)
    init_Q!(Q, list_w[1])
    compute_coeff!(list_coeff[1], list_w[1], list_r[1], Q, μ)
    for i in 2:nw
        update_Q!(Q, list_w[i], list_r[i])
        compute_coeff!(list_coeff[i], list_w[i], list_r[i], Q, μ)
    end
end

# function check_constraints(nw::S, list_w::Array{tuple_w{T, S}, 1}, list_r::Array{tuple_r{T, S}, 1})

# end

function loop!(n::S, nw::S, list_w::Array{tuple_w{T, S}, 1}, list_r::Array{tuple_r{T, S}, 1}, list_coeff::Array{struct_coefficients{T, S}, 1}, μ::Vector{T}) where {S<:Integer, T<:Real}
    for k in 1:n
        BFP!(nw, list_w, list_r, list_coeff, μ)
        FFP!(nw, list_w, list_coeff)
        plot_graph(
            k, 
            [list_w[i].x[j] for j in 1:list_w[1].nx, i in 1:nw],
            [list_w[i].u[j] for j in 1:list_w[1].nu, i in 1:nw], 
            [list_w[i].s[j] for j in 1:list_w[1].ns, i in 1:nw], 
            [list_w[i].y[j] for j in 1:list_w[1].ns, i in 1:nw], 
            [0.01 * i for i in 1:nw]
        )
    end
end

function main()
    nw = 5500
    nx = 3
    nu = 4
    ns = 7
    list_w = [tuple_w(nx, nu, ns, Float64) for _ in 1:nw]
    list_r = [tuple_r(ns, Float64) for _ in 1:nw]
    list_coeff = [struct_coefficients(nx, nu, ns, Float64) for _ in 1:nw]
    μ = ones(ns)
    loop!(10, nw, list_w, list_r, list_coeff, μ)
end


function plot_graph(index, plot_x, plot_u, plot_s, plot_y, plot_t)

    plots_x = Plots.plot(
        plot_t[:], 
        [ plot_x[1, :] plot_x[2, :] plot_x[3, :] ], 
        label=["r" "v_r" "v_t"],
        xlabel = "t",
        st=:scatter,
    )

    plots_u = Plots.plot(
        plot_t, 
        [ plot_u[1, :].-plot_u[2, :] plot_u[3, :].-plot_u[4, :] ], 
        label=["u_r" "u_t"],
        xlabel = "t",
        st=:scatter
    )

    plots_s = Plots.plot(
        plot_t, 
        [ plot_s[i, :] for i in 1:7], 
        label=["s$i" for i in 1:7],
        xlabel = "t",
        st=:scatter
    )

    plots_y = Plots.plot(
        plot_t, 
        [ plot_y[i, :] for i in 1:7], 
        label=["y$i" for i in 1:7],
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