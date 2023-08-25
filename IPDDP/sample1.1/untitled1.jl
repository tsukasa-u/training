using LinearAlgebra
using ForwardDiff
using Plots
using CPUTime


c(x, u) = begin
    rx, ry = x
    ux, uy= u
    return [
        ux - 1.5,
        -ux - 1.5,
        uy - 1.5,
        -uy - 1.5,
        ry - 1.0,
        -ry - 1.0,
        (1.0^2 - ((rx + 5.0)^2 + (ry + 1.0)^2)),
        (0.5^2 - ((rx + 8.0)^2 + (ry - 0.2)^2)),
        (1.5^2 - ((rx + 2.5)^2 + (ry - 1.0)^2)),
    ]
end
cx(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(x -> c(x, u), x)
cu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(u -> c(x, u), u)
cxx(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(x -> ForwardDiff.jacobian(x -> c(x, u), x), x)
cux(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(x -> ForwardDiff.jacobian(u -> c(x, u), u), x)
# fxu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(u -> ForwardDiff.jacobian(x -> c(x, u), x), u)
cuu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(u -> ForwardDiff.jacobian(u -> c(x, u), u), u)

f(x, u) = begin
    rx, ry= x
    ux, uy = u
    h = 1E-0
    v = 1E-0
    return [
        rx + h*ux,
        ry + v*uy
    ]
end
fx(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(x -> f(x, u), x)
fu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(u -> f(x, u), u)
fxx(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(x -> ForwardDiff.jacobian(x -> f(x, u), x), x)
fux(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(x -> ForwardDiff.jacobian(u -> f(x, u), u), x)
# fxu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(u -> ForwardDiff.jacobian(x -> f(x, u), x), u)
fuu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(u -> ForwardDiff.jacobian(u -> f(x, u), u), u)

l(x, u) = begin
    rx, ry = x
    ux, uy = u
    return ux^2 + uy^2
end
lx(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.gradient(x -> l(x, u), x)
lu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.gradient(u -> l(x, u), u)
lxx(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.hessian(x -> l(x, u), x)
lux(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(x -> ForwardDiff.gradient(u -> l(x, u), u), x)
# lxu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.jacobian(u -> ForwardDiff.gradient(x -> l(x, u), x), u)
luu(x::Vector{T}, u::Vector{T}) where T<:Real = ForwardDiff.hessian(u -> l(x, u), u)

lf(x) = begin
    rx, ry = x
    return (rx-0.0)^2+(ry-0.0)^2
end

vN(x) = begin
    return lf(x)
end
vNx(x::Vector{T}) where T<:Real = ForwardDiff.gradient(x -> vN(x), x)
vNxx(x::Vector{T}) where T<:Real = ForwardDiff.hessian(x -> vN(x), x)

mutable struct struct_V{T<:Real, S<:Integer}
    nx::S
    
    ΔV::T

    Vx::Matrix{T}

    Vxx::Matrix{T}

    function struct_V(_nx::S, T::DataType) where S
        self = new{T, S}()

        self.nx = _nx

        self.Vx = zeros(_nx, 1)
        self.Vxx = zeros(_nx, _nx)

        return self
    end
end

mutable struct struct_Q{T<:Real, S<:Integer}
    ns::S
    nx::S
    nu::S

    Qx::Matrix{T}
    Qu::Matrix{T}
    Qs::Matrix{T}

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

        self.Qx = zeros(_nx, 1)
        self.Qu = zeros(_nu, 1)
        self.Qs = zeros(_ns, 1)

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
        self.s = 0.1*ones(_ns)
        self.y = 0.01*ones(_ns)

        return self
    end
end

mutable struct array_w{T<:Real, S<:Integer}
    nw::S
    ns::S
    nx::S
    nu::S
    a::Vector{tuple_w{T, S}}
    b::Vector{tuple_w{T, S}}
    new::Base.RefValue{Vector{tuple_w{T, S}}}
    old::Base.RefValue{Vector{tuple_w{T, S}}}
    flag::Bool
    swap::Function

    function array_w(_nx::S, _nu::S, _ns::S, _nw::S, T::DataType) where S
        self = new{T, S}()
        
        self.flag = true

        self.nw = _nw
        self.nx = _nx
        self.nu = _nu
        self.ns = _ns

        self.a = [tuple_w(_nx, _nu, _ns, T) for _ in 1:_nw]
        self.b = [tuple_w(_nx, _nu, _ns, T) for _ in 1:_nw]
        # self.b = copy(self.a)
        
        self.new = Ref{Vector{tuple_w{T, S}}}(self.a)
        self.old = Ref{Vector{tuple_w{T, S}}}(self.b)
        
        self.swap = function()
            if self.flag
                self.new = Ref{Vector{tuple_w{T, S}}}(self.b)
                self.old = Ref{Vector{tuple_w{T, S}}}(self.a)
            else
                self.new = Ref{Vector{tuple_w{T, S}}}(self.a)
                self.old = Ref{Vector{tuple_w{T, S}}}(self.b)
            end
            self.flag = !self.flag
        end

        return self
    end
end

function Base.getindex(U::Base.RefValue{Vector{tuple_w{T, S}}}, i::Int) where {T<:Real, S<:Integer}
    1 <= i <= length(U[]) || throw(BoundsError(U, i))
    return U[][i]
end

function Base.setindex!(U::Base.RefValue{Vector{tuple_w{T, S}}}, j::T, i::Int) where {T<:Real, S<:Integer}
    1 <= i <= length(U[]) || throw(BoundsError(U, i))
    return U[][i] = j
end

Base.firstindex(U::Base.RefValue{Vector{tuple_w{T, S}}}) where {T<:Real, S<:Integer} = 1

Base.lastindex(U::Base.RefValue{Vector{tuple_w{T, S}}}) where {T<:Real, S<:Integer} = length(U[])

function Base.getindex(U::array_w{T, S}, i::Int) where {T<:Real, S<:Integer}
    1 <= i <= U.nw || throw(BoundsError(U, i))
    return U.new[][i]
end

function Base.setindex!(U::array_w{T, S}, j::T, i::Int) where {T<:Real, S<:Integer}
    1 <= i <= U.nw || throw(BoundsError(U, i))
    return U.new[][i] = j
end

Base.firstindex(U::array_w{T, S}) where {T<:Real, S<:Integer} = 1

Base.lastindex(U::array_w{T, S}) where {T<:Real, S<:Integer} = U.nw

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

    α::Matrix{T}
    η::Matrix{T}
    χ::Matrix{T}
    β::Matrix{T}
    θ::Matrix{T}
    ζ::Matrix{T}

    function struct_coefficients(_nx::S, _nu::S, _ns::S, T::DataType) where S
        self = new{T, S}()

        self.nx = _nx
        self.nu = _nu
        self.ns = _ns

        self.α = zeros(_nu, 1)
        self.η = zeros(_ns, 1)
        self.χ = zeros(_ns, 1)

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

function reshape_tensor_3(x, n::Tuple{S, S, S}) where S<:Integer
    return [x[(i-1)*n[1]+1, j] for v in 1:n[1], i in 1:n[2], j in 1:n[3]]
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
    _fxx::Array{T, 3} = reshape_tensor_3(fxx(w.x, w.u), (_Q.nx, _Q.nx, _Q.nx))
    _fux::Array{T, 3} = reshape_tensor_3(fux(w.x, w.u), (_Q.nx, _Q.nu, _Q.nx))
    _fuu::Array{T, 3} = reshape_tensor_3(fuu(w.x, w.u), (_Q.nx, _Q.nu, _Q.nu))
    _cxx::Array{T, 3} = reshape_tensor_3(cxx(w.x, w.u), (_Q.ns, _Q.nx, _Q.nx))
    _cux::Array{T, 3} = reshape_tensor_3(cux(w.x, w.u), (_Q.ns, _Q.nu, _Q.nx))
    _cuu::Array{T, 3} = reshape_tensor_3(cuu(w.x, w.u), (_Q.ns, _Q.nu, _Q.nu))
    # _fxx::Array{T, 3} = reshape(fxx(w.x, w.u), (_Q.nx, _Q.nx, _Q.nx))
    # _fux::Array{T, 3} = reshape(fux(w.x, w.u), (_Q.nx, _Q.nu, _Q.nx))
    # _fuu::Array{T, 3} = reshape(fuu(w.x, w.u), (_Q.nx, _Q.nu, _Q.nu))
    # _cxx::Array{T, 3} = reshape(cxx(w.x, w.u), (_Q.ns, _Q.nx, _Q.nx))
    # _cux::Array{T, 3} = reshape(cux(w.x, w.u), (_Q.ns, _Q.nu, _Q.nx))
    # _cuu::Array{T, 3} = reshape(cuu(w.x, w.u), (_Q.ns, _Q.nu, _Q.nu))
    
    # _fxx::Array{T, 3} = reshape(fxx(w.x, w.u), (_Q.nx, _Q.nx, _Q.nx))
    # _fux::Array{T, 3} = reshape(fux(w.x, w.u), (_Q.nu, _Q.nx, _Q.nx))
    # _fuu::Array{T, 3} = reshape(fuu(w.x, w.u), (_Q.nu, _Q.nx, _Q.nu))
    # _cxx::Array{T, 3} = reshape(cxx(w.x, w.u), (_Q.nx, _Q.ns, _Q.nx))
    # _cux::Array{T, 3} = reshape(cux(w.x, w.u), (_Q.nu, _Q.ns, _Q.nx))
    # _cuu::Array{T, 3} = reshape(cuu(w.x, w.u), (_Q.nu, _Q.ns, _Q.nu))

    # println(size(_cux), " : ", size(cux(w.x, w.u)), _cux)

    _Q.Qx[:] .= lx(w.x, w.u) + _fx_T*_V.Vx
    _Q.Qu[:] .= lu(w.x, w.u) + _fu_T*_V.Vx

    _Q.Qxx[:, :] .=           lxx(w.x, w.u)  + _fx_T*_V.Vxx*_fx + sum(_V.Vx[i] *           _fxx[i, :, :]  for i in 1:_Q.nx)
    # _Q.Qxu[:, :] .= lxu(w.x, w.u) + _fx_T*_V.Vxx*_fu + sum(_V.Vx[i] * _fxu[i, :, :] for i in 1:_Q.nx)
    _Q.Qxu[:, :] .= transpose(lux(w.x, w.u)) + _fx_T*_V.Vxx*_fu + sum(_V.Vx[i] * transpose(_fux[i, :, :]) for i in 1:_Q.nx)
    _Q.Quu[:, :] .=           luu(w.x, w.u)  + _fu_T*_V.Vxx*_fu + sum(_V.Vx[i] *           _fuu[i, :, :]  for i in 1:_Q.nx)
    
    # _Q.Qxx[:, :] .=           lxx(w.x, w.u)  + _fx_T*_V.Vxx*_fx + sum(_V.Vx[i] *           _fxx[:, i, :]  for i in 1:_Q.nx)
    # _Q.Qxu[:, :] .= transpose(lux(w.x, w.u)) + _fx_T*_V.Vxx*_fu + sum(_V.Vx[i] * transpose(_fux[:, i, :]) for i in 1:_Q.nx)
    # _Q.Quu[:, :] .=           luu(w.x, w.u)  + _fu_T*_V.Vxx*_fu + sum(_V.Vx[i] *           _fuu[:, i, :]  for i in 1:_Q.nx)


    _Q.Qx[:] .+= transpose(_Q.Qsx)*w.s
    _Q.Qu[:] .+= transpose(_Q.Qsu)*w.s
    _Q.Qxx[:, :] .+= sum(w.s[i] *           _cxx[i, :, :]  for i in 1:_Q.ns)
    _Q.Qxu[:, :] .+= sum(w.s[i] * transpose(_cux[i, :, :]) for i in 1:_Q.ns)
    _Q.Quu[:, :] .+= sum(w.s[i] *           _cuu[i, :, :]  for i in 1:_Q.ns)
    # _Q.Qxx[:, :] .+= sum(w.s[i] *           _cxx[:, i, :]  for i in 1:_Q.ns)
    # _Q.Qxu[:, :] .+= sum(w.s[i] * transpose(_cux[:, i, :]) for i in 1:_Q.ns)
    # _Q.Quu[:, :] .+= sum(w.s[i] *           _cuu[:, i, :]  for i in 1:_Q.ns)

end

φ(x, w::tuple_w, coeff::struct_coefficients, param::T) where T = w.u + param*coeff.α + coeff.β*(x - w.x)
ψ(x, w::tuple_w, coeff::struct_coefficients, param::T) where T = w.s + param*coeff.η + coeff.θ*(x - w.x)
ξ(x, w::tuple_w, coeff::struct_coefficients, param::T) where T = w.y + param*coeff.χ + coeff.ζ*(x - w.x)

# φ(x, w::tuple_w, coeff::struct_coefficients, param::T) where T = w.u + param*coeff.α + coeff.β*(x - w.x)
# ψ(x, w::tuple_w, coeff::struct_coefficients, param::T) where T = w.s + param*coeff.η + coeff.θ*(x - w.x)
# ξ(x, w::tuple_w, coeff::struct_coefficients, param::T) where T = w.y + param*coeff.χ + coeff.ζ*(x - w.x)

ReLU(x, a::T) where T<:Real = begin
    if x < -10/a 
        1E-9
    elseif x > 10/a
        x
    else
        0.5*(x + log(exp(a*x) + exp(-a*x))/a)
    end
end
function compute_coeff!(coeff::struct_coefficients, w::tuple_w, r::tuple_r, Q::struct_Q, μ::Vector{T}) where T<:Real
    _S::Matrix{T} = diagm(w.s)
    Y::Matrix{T} = diagm(w.y)
    #---------------------------------------------------
    I = ones(w.ns, w.ns)
    zero_us = zeros(w.nu, w.ns)
    zero_su = zeros(w.ns, w.nu)
    zero_ss = zeros(w.ns, w.ns)
    zero_sx = zeros(w.ns, w.nx)
    #---------------------------------------------------
    r.rp[:] .= c(w.x, w.u) + w.y
    r.rd[:] .= _S*w.y - μ
    r.rhat[:] .= _S*r.rp - r.rd

    # P = [Q.Quu transpose(Q.Qsu) zero_us; Q.Qsu zero_ss I; zero_su Y _S]
    # ret = -inv(P)*[Q.Qu Q.Qux; rp Q.Qsx; rd zero_sx]

    # coeff.α .= vec(ret[1:w.nx, 1])
    # coeff.η .= vec(ret[(w.nx + 1):(w.nx + w.nu), 1])
    # coeff.χ .= vec(ret[(w.nx + w.nu + 1):end, 1])

    # coeff.β .= ret[1:w.nx, 2:end]
    # coeff.θ .= ret[(w.nx + 1):(w.nx + w.nu), 2:end]
    # coeff.ζ .= ret[(w.nx + w.nu + 1):end, 2:end]
#---------------------------------------------------

    Y_1::Matrix{T} = inv(Y)
    SY_1::Matrix{T} = _S*Y_1

    R::Matrix{T} = (Q.Quu + transpose(Q.Qsu)*SY_1*Q.Qsu)
    ret::Matrix{T} = -inv(R)*[Q.Qu + transpose(Q.Qsu)*Y_1*r.rhat transpose(Q.Qxu) + transpose(Q.Qsu)*SY_1*Q.Qsx]
    
    coeff.α[:] .= ret[:, 1]

    coeff.β[:, :] .= ret[:, 2:end]

    coeff.η[:] .= Y_1*(r.rhat + _S*Q.Qsu*coeff.α)
    coeff.θ[:, :] .= SY_1*(Q.Qsx + Q.Qsu*coeff.β)

    coeff.χ[:] .= -r.rp - Q.Qsu*coeff.α
    coeff.ζ[:, :] .= -Q.Qsx - Q.Qsu*coeff.β
#---------------------------------------------------

    P = [Q.Quu transpose(Q.Qsu) zero_us; Q.Qsu zero_ss I; zero_su Y _S]
    println(size(P))
    println(size(Q.Qu), " : ", size(Q.Qxu), " : ", size(Q.Qsx), " : ", size(zero_sx), " : ", size(r.rp), " : ", size(r.rd))
    R = [transpose(Q.Qu) Q.Qxu; r.rp transpose(Q.Qsx); r.rd transpose(zero_sx)]
    C = [coff.α coff.β; coff.η coff.θ; coff.χ coff.ζ]

    println(P*C + R)

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

function FFP!(nw::S, list_w::array_w{T, S}, list_coeff::Array{struct_coefficients{T, S}, 1}) where {S<:Integer, T<:Real}
    # x::Vector{T} = copy(list_w[1].x)
    # u::Vector{T} = zeros(list_w[1].nu)
    # s::Vector{T} = zeros(list_w[1].ns)
    # y::Vector{T} = zeros(list_w[1].ns)

    isfailed::Bool = false
    failed_step::S = -1

    println("start")
    # println(list_w.old)

    # for i in 1:nw-1
    #     if any(list_w.new[i].s .< 0) || any(list_w.new[i].y .< 0)
    #         println("failed")
    #     end
    # end

    param_steps::Vector{T} = [2.0^i for i in 0:-1:-20]
    
    for param_step in param_steps
    # for param_step in [2.0^i for i in -3:-1:-3]
        isfailed = false
        for i in 1:nw-1
            # list_w.new[i].s[:] = ReLU.(ψ(list_w.new[i].x, list_w.old[i], list_coeff[i], param_step), 1E+0)
            # list_w.new[i].y[:] = ReLU.(ξ(list_w.new[i].x, list_w.old[i], list_coeff[i], param_step), 1E+0)
            list_w.new[i].s[:] = ψ(list_w.new[i].x, list_w.old[i], list_coeff[i], param_step)
            list_w.new[i].y[:] = ξ(list_w.new[i].x, list_w.old[i], list_coeff[i], param_step)
    
            if any(list_w.new[i].s .<= 1E-3) || any(list_w.new[i].y .<= 1E-3)
            # if any(list_w.new[i].s .< 0.01*list_w.old[i].s) || any(list_w.new[i].y .< 0.01*list_w.old[i].y) || any(list_w.new[i].s .< 0) || any(list_w.new[i].y .< 0)
            # if any(list_w.new[i].s .< 0.01*list_w.old[i].s) || any(list_w.new[i].y .< 0.01*list_w.old[i].y) || any(list_w.new[i].s .<= 1E-3) || any(list_w.new[i].y .<= 1E-3)
                isfailed = true
                failed_step = i
            #     println("step : ", i, " ", param_step, " ", list_w.new[i].x - list_w.old[i].x)
            #     println("α : ", list_coeff[i].α)
            #     println("β : ", list_coeff[i].β)
            #     println("η : ", list_coeff[i].η)
            #     println("θ : ", list_coeff[i].θ)
            #     println("χ : ", list_coeff[i].χ)
            #     println("ζ : ", list_coeff[i].ζ)  
            #     # println("failed:" , i, " : ",  param_step, " : ", list_w.new[i].s, " : ", list_w.new[i].y, " : ", list_w.old[i].s, " : ", list_w.old[i].y)
                break
            end
    
            list_w.new[i].u[:] = φ(list_w.new[i].x, list_w.old[i], list_coeff[i], param_step)
    
            # list_w[i].u[:] .= u
            # list_w[i].s[:] .= s
            # list_w[i].y[:] .= y
            # list_w[i].x[:] .= x
    
            # check_constraints(list_w[i])
    
            list_w.new[i + 1].x[:] = f(list_w.new[i].x, list_w.new[i].u)
        end

        if !isfailed
            break
        end
    end
    if isfailed
        # for i in failed_step:failed_step
        for i in failed_step:nw-1
            # println("step : ", i, " ", param_steps[end], " ", list_w.new[i].x - list_w.old[i].x, " ", ψ(list_w.new[i].x, list_w.old[i], list_coeff[i], param_steps[end]))
            list_w.new[i].s[:] = ReLU.(ψ(list_w.new[i].x, list_w.old[i], list_coeff[i], param_steps[end]), 1E+1)
            # println("s : ", list_w.new[i].s)
            list_w.new[i].y[:] = ReLU.(ξ(list_w.new[i].x, list_w.old[i], list_coeff[i], param_steps[end]), 1E+1)
    
            list_w.new[i].u[:] = φ(list_w.new[i].x, list_w.old[i], list_coeff[i], param_steps[end])
    
            list_w.new[i + 1].x[:] = f(list_w.new[i].x, list_w.new[i].u)
        end
    end
    println("isfailed: ", isfailed)
end

function BFP!(nw::S, list_w::array_w{T, S}, list_r::Array{tuple_r{T, S}, 1}, list_coeff::Array{struct_coefficients{T, S}, 1}, μ::Vector{T}) where {S<:Integer, T<:Real}

    Q::struct_Q = struct_Q(list_w[1].nx, list_w[1].nu, list_w[1].ns, T)
    V::struct_V = struct_V(list_w[1].nx, T)
    init_V!(V, list_w[nw])
    for i in nw:-1:1
        compute_Q!(Q, V, list_w[i])
        compute_coeff!(list_coeff[i], list_w[i], list_r[i], Q, μ)
        update_Q!(Q, list_w[i], list_r[i])
        update_V!(Q, V, list_coeff[i])
    end
end

function check_constraints(list_w::array_w{T, S}) where {S<:Integer, T<:Real}
    for i in 1:list_w.nw
        list_w.new[i].s[:] .= max.(1E-3::T, list_w.new[i].s)
        list_w.new[i].y[:] .= max.(1E-3::T, list_w.new[i].y)
    end
end

function loop!(n::S, nw::S, list_w::array_w{T, S}, list_r::Array{tuple_r{T, S}, 1}, list_coeff::Array{struct_coefficients{T, S}, 1}, μ::Vector{T}) where {S<:Integer, T<:Real}
    
    wrap_plot_graph(0, list_w, nw)
    for k in 1:n
        check_constraints(list_w)
        BFP!(nw, list_w, list_r, list_coeff, μ)
        list_w.swap()
        # check_constraints(list_w)
        FFP!(nw, list_w, list_coeff)
        update_μ!(μ, 5.0)
        
        tmp = [sqrt(list_coeff[i].α[1]^2 + list_coeff[i].α[2]^2) for i in 1:nw]
        vec1 = [[list_coeff[i].α[1]/tmp[i] for i in 1:nw], [list_coeff[i].α[2]/tmp[i] for i in 1:nw]]
        tmp = [sqrt(list_coeff[i].β[1]^2 + list_coeff[i].β[2]^2) for i in 1:nw]
        vec2 = [[list_coeff[i].β[1]/tmp[i] for i in 1:nw], [list_coeff[i].β[2]/tmp[i] for i in 1:nw]]
        
        wrap_plot_graph(k, list_w, nw, vec1=vec1[:], vec2=vec2[:])
        wrap_plot_graph(-k, list_w, nw, check_nan=true)
    end
end

function wrap_plot_graph(idx, list_w::array_w{T, S}, nw::S; check_nan=false, vec1::Array{Array{T, 1}, 1}=[T[],T[]], vec2::Array{Array{T, 1}, 1}=[T[],T[]]) where {S<:Integer, T<:Real}
    
    if check_nan
        plot_graph(
            idx, 
            [isnan(list_w[i].x[j]) ? j : isinf(list_w[i].x[j]) ? -j : NaN for j in 1:list_w[1].nx, i in 1:nw],
            [isnan(list_w[i].u[j]) ? j : isinf(list_w[i].u[j]) ? -j : NaN for j in 1:list_w[1].nu, i in 1:nw], 
            [isnan(list_w[i].s[j]) ? j : isinf(list_w[i].s[j]) ? -j : NaN for j in 1:list_w[1].ns, i in 1:nw], 
            [isnan(list_w[i].y[j]) ? j : isinf(list_w[i].y[j]) ? -j : NaN for j in 1:list_w[1].ns, i in 1:nw], 
            [0.01 * i for i in 1:nw],
            vec1 = vec1[:],
            vec2 = vec2[:]
        )
    else
        plot_graph(
            idx, 
            [list_w[i].x[j] for j in 1:list_w[1].nx, i in 1:nw],
            [list_w[i].u[j] for j in 1:list_w[1].nu, i in 1:nw], 
            [list_w[i].s[j] for j in 1:list_w[1].ns, i in 1:nw], 
            [list_w[i].y[j] for j in 1:list_w[1].ns, i in 1:nw], 
            [0.01 * i for i in 1:nw],
            vec1 = vec1[:],
            vec2 = vec2[:]
        )
    end
end

function init_μ!(μ::Vector{T}, list_w::array_w{T, S}, nw::S, ns::S) where {S<:Integer, T<:Real}
    μ[:] .= (sum(l(list_w[i].x, list_w[i].u) for i in 1:nw)/nw/ns + lf(list_w[end].x))*zeros(ns)
end

function update_μ!(μ::Vector{T}, κ::T) where T<:Real
    μ[:] .= min.(μ/κ, μ.^1.2)
end

function init_xu!(list_w::array_w{T, S}, bx::Tuple{Array{T, 1}, Array{T, 1}}, bu::Tuple{Array{T, 1}, Array{T, 1}}) where {S<:Integer, T<:Real}
    bx0, bxf = bx
    bx0f = bxf .- bx0
    bu0, buf = bu
    bu0f = buf .- bu0
    for i in 1:list_w.nw
        list_w.new[i].x .= bx0f*i/list_w.nw .+ bx0
        list_w.new[i].u .= bu0f*i/list_w.nw .+ bu0
        list_w.old[i].x .= bx0f*i/list_w.nw .+ bx0
        list_w.old[i].u .= bu0f*i/list_w.nw .+ bu0
    end
end

function get_list_init(_nx::S, _nu::S, _ns::S, _nw::S, T::DataType) where {S <: Integer}
    T <: Real || throw(ArgumentError("T must be a subtype of Real"))
    list_w::array_w{T, S} = array_w(_nx, _nu, _ns, _nw, T)
    list_r::Array{tuple_r{T, S}, 1} = [tuple_r(_ns, T) for _ in 1:_nw]
    list_coeff::Array{struct_coefficients{T, S}, 1} = [struct_coefficients(_nx, _nu, _ns, T) for _ in 1:_nw]
    μ::Vector{T} = zeros(_ns)
    return list_w, list_r, list_coeff, μ
end

function main()
    nw::Int64 = 20
    nx::Int64 = 2
    nu::Int64 = 2
    ns::Int64 = 9

    list_w, list_r, list_coeff, μ = get_list_init(nx, nu, ns, nw, Float64)

    init_xu!(list_w, ([-10.0, 0.0], [0.0, 0.0]), ([0.0, 0.0], [0.0, 0.0]))
    init_μ!(μ, list_w, nw, ns)

    loop!(10, nw, list_w, list_r, list_coeff, μ)
end


function plot_graph(index, plot_x, plot_u, plot_s, plot_y, plot_t; vec1::Array{Array{T, 1}, 1}=[T[],T[]], vec2::Array{Array{T, 1}, 1}=[T[],T[]]) where T<:Real

    plots_x = Plots.plot(
        plot_x[1, :], 
        plot_x[2, :], 
        st=:scatter,
        legend = :none
    )
    plot!([-5.5 + sqrt(1  )*cos(i) for i in 0:0.01:6.28], [-1.0 + sqrt(1  )*sin(i) for i in 0:0.01:6.28])
    plot!([-8.0 + sqrt(0.5)*cos(i) for i in 0:0.01:6.28], [ 0.2 + sqrt(0.5)*sin(i) for i in 0:0.01:6.28])
    plot!([-2.5 + sqrt(1.5)*cos(i) for i in 0:0.01:6.28], [ 1.0 + sqrt(1.5)*sin(i) for i in 0:0.01:6.28])

    
    if !isempty(vec1) && size(vec1[1], 1) > 0
        quiver!(
            plot_x[1, :],
            plot_x[2, :],
            quiver=(vec1[1], vec1[2])
        )
    end

    if !isempty(vec2) && size(vec2[1], 1) > 0
        quiver!(
            plot_x[1, :] + vec1[1][:],
            plot_x[2, :] + vec1[2][:],
            quiver=(vec2[1], vec2[2])
        )
    end

    plots_u = Plots.plot(
        plot_t, 
        [ plot_u[i, :] for i in 1:2], 
        ylabel = "u",
        xlabel = "t",
        st=:scatter,
        legend = :none
    )

    plots_s = Plots.plot(
        plot_t, 
        [ plot_s[i, :] for i in 1:9], 
        ylabel = "s",
        xlabel = "t",
        st=:scatter,
        legend = :none
    )

    plots_y = Plots.plot(
        plot_t, 
        [ plot_y[i, :] for i in 1:9], 
        ylabel = "y",
        xlabel = "t",
        st=:scatter,
        legend = :none
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
    savefig(string(index, base = 10, pad = 2)*".svg")
end


@time @CPUtime main()