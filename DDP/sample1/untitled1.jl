using LinearAlgebra
using ForwardDiff
using Plots
using CPUTime

f(x, u) = begin
    r, v_r, v_t = x
    u_rp, u_rm, u_tp, u_tm = u
    return [
        (v_r),
        (v_t^2.0/r - 1.0/r^2.0 + u_rp - u_rm),
        (-v_r*v_t/r + u_tp - u_tm)
    ]
end
fx(x::Vector{Real}, u::Vector{Real}) = ForwardDiff.jacobian(x -> f(x, u), x)
fu(x::Vector{Real}, u::Vector{Real}) = ForwardDiff.jacobian(u -> f(x, u), u)
fxx(x::Vector{Real}, u::Vector{Real}) = ForwardDiff.jacobian(x -> ForwardDiff.jacobian(x -> f(x, u), x), x)
fux(x::Vector{Real}, u::Vector{Real}) = ForwardDiff.jacobian(x -> ForwardDiff.jacobian(u -> f(x, u), u), x)
fuu(x::Vector{Real}, u::Vector{Real}) = ForwardDiff.jacobian(u -> ForwardDiff.jacobian(u -> f(x, u), u), u)

l(x, u) = begin
    r, v_r, v_t = x
    u_rp, u_rm, u_tp, u_tm = u
    return u_rp + u_rm + u_tp + u_tm + (u_rp*u_rm + u_tp*u_tm)
end
lx(x::Vector{Real}, u::Vector{Real}) = ForwardDiff.gradient(x -> l(x, u), x)
lu(x::Vector{Real}, u::Vector{Real}) = ForwardDiff.gradient(u -> l(x, u), u)
lxx(x::Vector{Real}, u::Vector{Real}) = ForwardDiff.hessian(x -> l(x, u), x)
lux(x::Vector{Real}, u::Vector{Real}) = ForwardDiff.jacobian(x -> ForwardDiff.gradient(u -> l(x, u), u), x)
luu(x::Vector{Real}, u::Vector{Real}) = ForwardDiff.hessian(u -> l(x, u), u)

lf(x) = sum((x .- [4.0, 0.0, 0.5]).^2.0)
# lf(x) = sum((x .- [4.0, 0.0, 0.5]).^2.0)/sum(x.^2.0)

# v(i, x, u) = sum(l(_x, _u) for (_x, _u) in zip(x[i:end], u[i:end]); init = 0.0) + lf(x[end])
v(i, x, u, n) = begin
    if i == n
        return l(x, u[i]) + lf(x)
    else
        return l(x, u[i]) + v(i + 1, f(x, u[i]), u, n)
    end
end
vx(i, x::Vector{Real}, u::Vector{Vector{Real}}, n) = ForwardDiff.gradient(x -> v(i, x, u, n), x)
vxx(i, x::Vector{Real}, u::Vector{Vector{Real}}, n) = ForwardDiff.hessian(x -> v(i, x, u, n), x)

function compute_Q(idx, x, u; __vxx = [;], __vx = [;])

    # println(__vxx)
    # println(__vx)

    nx = size(x[1], 1)
    nu = size(u[1], 1)
    n = size(x, 1)

    _lx = lx(x[idx], u[idx])
    _lu = lu(x[idx], u[idx])
    _lxx = lxx(x[idx], u[idx])
    _lux = lux(x[idx], u[idx])
    _luu = luu(x[idx], u[idx])

    # println(isassigned(__vxx))
    # println(isempty(__vxx))
    # println(size(__vx))
    if isassigned(__vxx) && isassigned(__vx)
        _vx = __vx
        _vxx = __vxx
    else
        _vx = vx(idx + 1, x[idx], u, n)
        _vxx = vxx(idx + 1, x[idx], u, n)
    end

    _fx = fx(x[idx], u[idx])
    _fu = fu(x[idx], u[idx])
    _fxx = reshape(fxx(x[idx], u[idx]), (nx, nx, nx))
    _fux = reshape(fux(x[idx], u[idx]), (nx, nu, nx))
    _fuu = reshape(fuu(x[idx], u[idx]), (nx, nu, nu))

    Qx = _lx + transpose(_fx) * _vx
    Qu = _lu + transpose(_fu) * _vx

    # Qxx = _lxx + transpose(_fx) * _vxx * _fx + sum(_vx[i] * _fxx[i, :, :] for i in 1:nx)
    # Qux = _lux + transpose(_fu) * _vxx * _fx + sum(_vx[i] * _fux[i, :, :] for i in 1:nx)
    # Quu = _luu + transpose(_fu) * _vxx * _fu + sum(_vx[i] * _fuu[i, :, :] for i in 1:nx)

    Qxx = _lxx + transpose(_fx) * _vxx * _fx + sum(_vx[i] * _fxx[i, :, :] for i in 1:nx)
    Qux = _lux + transpose(_fu) * _vxx * _fx + sum(_vx[i] * _fux[i, :, :] for i in 1:nx)
    Quu = _luu + transpose(_fu) * _vxx * _fu + sum(_vx[i] * _fuu[i, :, :] for i in 1:nx)

    # println(Qx)
    # println(Qu)
    # println(Qxx)
    # println(Qux)
    # println(Quu)

    # Qxx = _lxx + transpose(_fx) * _vxx * _fx + kron(_vx, _fxx)
    # Qux = _lux + transpose(_fu) * _vxx * _fx + kron(_vx, _fux)
    # Quu = _luu + transpose(_fu) * _vxx * _fu + kron(_vx, _fuu)

    return Qx, Qu, Qxx, Qux, Quu
end

function inverse_Quu(Quu)
    @assert det(Quu) != 0.0 "non-regular matrix"
    return inv(Quu)
end

function compute_k(Qxx_1, Qu)
    return -Qxx_1 * Qu
end

function compute_K(Quu_1, Qux)
    return -Quu_1 * Qux
end

function compute_Vx(Qx, Qu, K)
    return Qx + transpose(K) * Qu
end

function compute_Vxx(Qxx, Qux, K)
    return Qxx + transpose(K) * Qux
end

function set_constraint!(x, xc, u, uc)
    x .= min.(max.(x, xc[1]), xc[2])
    u .= min.(max.(u, uc[1]), uc[2])
end

function BFP(n, x, u)

    list_k = []
    list_K = []
    
    Qx, Qu, Qxx, Qux, Quu = compute_Q(n - 1, x, u)
    Quu_1 = inverse_Quu(Quu)
    k = compute_k(Quu_1, Qu)
    K = compute_K(Quu_1, Qux)
    push!(list_k, k)
    push!(list_K, K)

    for i in (n - 2):-1:1

        Vxx = compute_Vxx(Qxx, Qux, K)
        Vx = compute_Vx(Qx, Qu, K)
        # Vx = Qx - transpose(Qux) * Quu_1 * Qu
        # Vxx = Qxx - transpose(Qux) * Quu_1 * Qux

        Qx, Qu, Qxx, Qux, Quu = compute_Q(i, x, u; __vxx = Vxx, __vx = Vx)
        Quu_1 = inverse_Quu(Quu)
        k = compute_k(Quu_1, Qu)
        K = compute_K(Quu_1, Qux)

        push!(list_k, k)
        push!(list_K, K)
    end

    return list_k, list_K
end

function FFP(n, x, u, xc, uc, α, dt, k, K)

    _x = copy(x)
    _u = copy(u)
    for i in 1:(n - 1)
        # println(u[i])
        # println(α*k[i])
        # println(K[i] * (_x[i] - x[i]))
        _u[i] = u[i] .+ α*k[i] .+ K[i] * (_x[i] - x[i])

        set_constraint!(_x[i], xc, _u[i], uc)

        _x[i + 1] = _x[i] + f(_x[i], _u[i])*dt
        # println(_u[i])
    end
    
    return _x, _u
end

function BFFP!(n, x, u, xc, uc, α, dt)

    k, K = BFP(n, x, u)
    _x, _u = FFP(n, x, u, xc, uc, α, dt, k, K)

    ret = maximum(maximum(abs.(_x[i] - x[i]) for i in 1:n)) < 1e-6
    ret = maximum(maximum(abs.(_x[i + 1] - f(_x[i], _u[i])) for i in 1:(n - 1))) < 1e-6 && ret

    x .= _x
    u .= _u

    return ret
end

function loop()
    n = 101
    t0 = 0
    tf = 55
    dt = (tf - t0)/(n - 1)

    nx = 3
    nu = 4

    x0 = [1.0, 0.0, 1.0]
    xf = [4.0, 0.0, 0.5]
    
    x = [Real[(xf[j] - x0[j])*i/(n - 1) + x0[j] for j in 1:nx] for i in 1:n]
    u = [Real[0.001 for j in 1:nu] for i in 1:n]

    xc = [[0.01, 0.0, 0.0], [10.0, 10.0, 10.0]]
    uc = [[0.0, 0.0, 0.0, 0.0], [0.01, 0.01, 0.01, 0.01]]

    plot_graph(0, [x[i][j] for j in 1:nx, i in 1:n], [u[i][j] for j in 1:nu, i in 1:n], [dt*i for i in 1:n])
    for idx in 1:100
        ret = BFFP!(n, x, u, xc, uc, 0.8, dt)
        # println(x)
        plot_graph(idx, [x[i][j] for j in 1:nx, i in 1:n], [u[i][j] for j in 1:nu, i in 1:n], [dt*i for i in 1:n])
        if ret
            break
        end
    end
end

function plot_graph(index, plot_x, plot_u, plot_t)

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
    
    plot_ref = Plots.plot(
        plots_x,
        plots_u,
        layout = (2, 1),
        # legend = false,
        # margin = 1Plots.cm,
    )
    png(string(index, base = 10, pad = 2))
end

@time @CPUtime loop()