using JuMP
# import Ipopt
using NLopt
using Plots

using Polynomials
using SpecialPolynomials
using SpecialFunctions

#-------------------------------------------------------

function compute_η(
        n::Int, 
        dx::Array{Real, 1}, 
        f::Array{Real, 1}, 
        w::Array{Real, 1}
    )::Real
        return sum(abs(dx[k] - f[k])*w[k] for k in 1:(n+1))
end

function compute_w(
        ns::Int, 
        n::Array{Int, 1}, 
        dx::Array{Array{Real, 1}, 1}, 
        x::Array{Array{Real, 1}, 1}
    )::Real
        return maximum([maximum([max(x[i][k], dx[i][k]) for k in 1:(n[i] + 1)]) for i in 1:ns])
end

function compute_ε(
        ns::Int, 
        n::Array{Int, 1}, 
        dx::Array{Array{Real, 1}, 1}, 
        x::Array{Array{Real, 1}, 1}, 
        f::Array{Array{Real, 1}, 1}, 
        W::Array{Array{Real, 1}, 1}
    )::Tuple{Array{Real, 1}, Real}
        w::Real = compute_w(ns, n, dx, x)
        return ([compute_η(n[i], dx[i], f[i], W[n[i]])/w for i in 1:ns], w)
end

function check_ε_tolerance(ε_tol::Real, ns::Real, ε::Array{Array{Real, 1}, 1})::Real
    tmp = maximum([sum(ele[i] for i in ns) for ele in ε])
    return tmp < ε_tol ? tmp : -1
end

function select_refine_segment(ε_tol::Real, α::Real, nx::Int, ε::Array{Array{Real, 1}, 1})::Array{Int, 1}
    ret_index::Array{Int, 1} = fill(-1, nx)
    for (j, ele) in enumerate(ε)
        idx = sortperm(ele)
        sum = 0
        r = 0
        ret_index[j] = idx[1]
        for i in idx
            if sum < α*ε_tol
                sum += ele[i]
                r = i
            else
                ret_index[j] = r
                break
            end
        end
    end
    return ret_index
end

function compute_difficulty(n::Int, f::Array{Real, 2}, w::Array{Real, 1}, ε::Array{Real, 1}, τ::Array{Real, 1}, t0::Real, tf::Real)::Array{Real, 1}
    F_k = [sum(abs(ε[j]*f[j, k]/w[j]) for (j, ele) in enumerate(ε)) for k in 1:(n+1)]
    F_1 = [(F_k[k+1] - F_k[k])*2/(tf-t0)/(τ[k+1] - τ[k]) for k in 1:n]
    insert!(F_1, 1, NaN)
    F_2 = [(F_k[k] - F_k[k-1])*2/(tf-t0)/(τ[k+1] - τ[k-1]) for k in 2:n]
    insert!(F_2, 1, NaN)

    sum_F = sum(F_2[l] for l in 2:n)
    G_k = [(n - 1)*abs(F_2[k])/sum_F for k in 2:n]
    insert!(G_k, 1, NaN)
    return G_k
end

function convert_coefficient_n_m(
    n::Real,
    m::Array{Real, 1},
    nx::Real,
    nu::Real,
    xn::Array{Real, 2}, 
    un::Array{Real, 2}, 
    τn::Array{Real, 1}, 
    τm::Array{Array{Real, 1}}, 
    tk::Array{Int, 1}
)::Tuple{Array{Array{Real, 2}, 1}, Array{Array{Real, 2}, 1}}
    
    # xm = tuple(Tuple(Real[]))
    xm = Array{Real, 2}[]
    um = Array{Real, 2}[]
    # for tk_idx in 2:length(tk)
    tk_size = size(tk, 1)
    for tk_idx in 2:tk_size
        _tk_idx = tk_idx - 1
        _tk = tk[_tk_idx]
        tk_ = tk[tk_idx]
        tau = τn[_tk:tk_]

        if size(tau, 1) > 0

            tau = tau.*(2.0/(tau[end] - tau[1]))
            tau = tau .- (1 + tau[1])
            # tau = tau .- tau[1]
            # tau = tau./tau[end]
            # println(tau)

            _xm = fill(0.0, (nx, m[_tk_idx] + 1))
            _um = fill(0.0, (nu, m[_tk_idx] + 1))

            i = 1
            for k in 2:m[_tk_idx]
                while τm[_tk_idx][k] < tau[i]
                    i -= 1
                end
                while !(τm[_tk_idx][k] > tau[i] && τm[_tk_idx][k] <= tau[i + 1])
                    i += 1
                end
                ratio_ = tau[i + 1]- tau[i]
                ratio_a = (τm[_tk_idx][k] - tau[i + 1])/ratio_
                ratio_b = (τm[_tk_idx][k] - tau[i])/ratio_
                for j in 1:nx
                    _xm[j, k] = xn[j, _tk + i]*ratio_b - xn[j, _tk + i - 1]*ratio_a
                end
                for j in 1:nu
                    _um[j, k] = un[j, _tk + i]*ratio_b - un[j, _tk + i - 1]*ratio_a
                end
            end
            for j in 1:nx
                _xm[j, 1] = xn[j, _tk]
                _xm[j, m[_tk_idx] + 1] = xn[j, tk_]
            end
            for j in 1:nu
                _um[j, 1] = un[j, _tk]
                _um[j, m[_tk_idx] + 1] = un[j, tk_]
            end
            push!(xm, _xm)
            push!(um, _um)
        end
    end
    return (xm, um)
end

function refine_segment(
        n::Int, 
        nx::Int, 
        nu::Int, 
        xn::Array{Real, 2},
        un::Array{Real, 2},
        f::Array{Real, 2}, 
        w::Array{Real}, 
        ε::Array{Real, 1}, 
        τ::Array{Array{Real, 1}, 1}, 
        t0::Real, 
        tf::Real, 
        g_tolerance::Real=2
    )::Tuple{Int, Array{Int, 1}, Array{Real, 1}, Array{Real, 1}, Array{Array{Real, 2}, 1}, Array{Array{Real, 2}, 1}}
        if n <= 40
            G_k = compute_difficulty(n, f, w, ε, τ[n], t0, tf)
            (max_idx, max_val) = findmax(G_k)
            if max_val > g_tolerance  #G_k tolerance 2
                _m = max(floor(n/3), 7)
                t1 = (tf - t0)*(max_idx - 2)/n + t0
                t2 = (tf - t0)*(max_idx    )/n + t0
                (rx, ru) = convert_coefficient_n_m(n, Real[_m, _m, _m], nx, nu, xn, un, τ[n], Array{Real, 1}[τ[_m], τ[_m], τ[_m]], Int[1, max_idx - 1, max_idx + 1, n+1])
                return (3, [_m, _m, _m], Real[t0, t1, t2], Real[t1, t2, tf], rx, ru)
            else
                _m = n + 8
                (rx, ru) = convert_coefficient_n_m(n, Real[_m], nx, nu, xn, un, τ[n], Array{Real, 1}[τ[_m]], Int[1, n+1])
                return (1, [_m], Real[t0], Real[tf], rx, ru)
            end
        else
            _m = Int(floor(n/2))
            t1 = (tf - t0)*(floor(0.4*n) - 1)/n + t0
            t2 = (tf - t0)*(floor(0.7*n) - 1)/n + t0
            (rx, ru) = convert_coefficient_n_m(n, Real[_m, _m, _m], nx, nu, xn, un, τ[n], Array{Real, 1}[τ[_m], τ[_m], τ[_m]], Int[1, Int(floor(0.4*n)), Int(floor(0.7*n)), n+1])
            return (3, [_m, _m, _m], Real[t0, t1, t2], Real[t1, t2, tf], rx, ru)
        end
end

function refine(
        ns::Int, 
        nx::Int,  
        nu::Int,  
        n::Array{Int, 1}, 
        x::Array{Array{Real, 2}, 1}, 
        u::Array{Array{Real, 2}, 1}, 
        dx::Array{Array{Real, 2}, 1}, 
        fx::Array{Array{Real, 2}, 1},
        τ::Array{Array{Real, 1}, 1},
        t0::Array{Real, 1},
        tf::Array{Real, 1},
        ε_tol::Real, 
        α::Real,
        W::Array{Array{Real, 1}, 1}
    )::Tuple{Int , Array{Int, 1}, Array{Array{Real, 2}, 1}, Array{Array{Real, 2}, 1}, Array{Real, 1}, Array{Real, 1}}
        # εx::Array{Array{Real, 1}, 1} = Array{Real, 1}[]
        # wx::Array{Real, 1} = Real[]
        εx = Array{Real, 1}[]
        wx = Real[]
        for j in 1:nx
            tmp = compute_ε(
                ns, 
                n, 
                Array{Real, 1}[[dx[i][j, k] for k in 1:(n[i] + 1)] for i in 1:ns], 
                Array{Real, 1}[[x[i][j, k] for k in 1:(n[i] + 1)] for i in 1:ns], 
                Array{Real, 1}[[fx[i][j, k] for k in 1:(n[i] + 1)] for i in 1:ns], 
                W
            )
            push!(εx, tmp[1])
            push!(wx, tmp[2])
        end

        if check_ε_tolerance(ε_tol, ns, εx) > 0
            return (-ns, n, x, u, t0, tf)
        end

        x_idx = select_refine_segment(ε_tol, α, nx, εx)

        rx = copy(x)
        ru = copy(u)
        rn = copy(n)
        rt0 = copy(t0)
        rtf = copy(tf)
        offset = 0

        # println(unique(sort(x_idx)))

        for idx in unique(sort(x_idx))
            # println(isempty(rx))
            deleteat!(rx, idx + offset)
            deleteat!(ru, idx + offset)
            deleteat!(rn, idx + offset)
            deleteat!(rt0, idx + offset)
            deleteat!(rtf, idx + offset)
            # println(isempty(rx))
            (_len, _m, _t0, _tf, segx, segu) = refine_segment(n[idx], nx, nu, x[idx], u[idx], fx[idx], wx, εx[:][idx], τ, t0[idx], tf[idx])
            # println(isempty(rx), " ", isempty(segx), " ", idx + offset, " ", _len, " ", size(segx, 1))
            for i in _len:-1:1
                # print(i)
                insert!(rx, idx + offset, segx[i])
                insert!(ru, idx + offset, segu[i])
                insert!(rn, idx + offset, _m[i])
                insert!(rt0, idx + offset, _t0[i])
                insert!(rtf, idx + offset, _tf[i])
            end
            
            # println(isempty(rx), " ", isempty(segx), " ", idx + offset, " ", _len, " ", size(segx, 1))

            offset = offset + _len - 1
        end

        # println(isempty(rx))
        return (ns + offset, rn, rx, ru, rt0, rtf)  # copied or returned pointer?
end

function Real_binomial(n::Real, k::Real)
    return gamma(1 + n)/gamma(1 + k)/gamma(1 + n - k)
end

function compute_W(n, α, β, τ::Array{Real, 1})
    
    P = Jacobi{α, β}
    Pn = basis(P, n)
    p_n = Pn.(τ)

    # fill(pi/n, n)
    _w = fill(0.0, n + 1)

    a_1 = 2^(α + β + 1)*gamma(α + 2)*gamma(β + 2)/gamma(α + β + 3)/n^2*Real_binomial(n + α, n - 1)*Real_binomial(n + β, n - 1)/Real_binomial(n + α + β + 1, n - 1)
    a_2 = (4*(n + α)*(n + β) + (α - β)^2)/(2*n + α + β)^2

    for i in 2:n
        _w[i] = a_1/(a_2 - τ[i]^2)*(1 - τ[i]^2)/p_n[i]^2
    end
    _w[1] = 2^(α + β + 1)*gamma(α + 2)*gamma(β + 1)/gamma(α + β + 3)*Real_binomial(n + α, n - 1)/Real_binomial(n + β, n - 1)/Real_binomial(n + α + β + 1, n - 1)
    _w[n + 1] = _w[1]
    return _w
end

function compute_tau(n)
    -cos.(Vector(0:n)*π/n)
end

function compute_D(n, alpha, beta, tau)
    P = Jacobi{alpha, beta}
    Pn = basis(P, n)
    p_n = Pn.(tau)

    tmp = fill(0.0, (n + 1, n + 1))
    tmp[1, 1] = (alpha - n*(n + alpha + beta + 1.0))/2.0/(beta + 2.0)
    tmp[n + 1, n + 1] = (n*(n + alpha + beta + 1.0) - beta)/2.0/(alpha + 2.0)
    for k = 2:n
        for j = (k + 1):(n + 1)
            tmp[k, j] = p_n[k]/p_n[j]/(tau[k] - tau[j])
            tmp[j, k] = p_n[j]/p_n[k]/(tau[j] - tau[k])
        end
    end

    for k = 2:n
        tmp[k, k] = ((alpha + beta)*tau[k] + alpha - beta)/2.0/(1.0 - tau[k]^2.0)
    end
        
    for j = 2:n
        tmp[1, j] = -p_n[1]/p_n[j]/(beta + 1.0)/(1.0 + tau[j])
    end

    for j = 2:n
        tmp[n + 1, j] = p_n[n + 1]/p_n[j]/(alpha + 1.0)/(1.0 - tau[j])
    end
    
    for k = 2:n
        tmp[k, 1] = p_n[k]/p_n[1]*(beta + 1.0)/(1.0 + tau[k])
    end
    
    for k = 2:n
        tmp[k, n + 1] = -p_n[k]/p_n[n + 1]*(alpha + 1.0)/(1.0 - tau[k]);
    end
    
    tmp[1, n + 1] = -(alpha + 1.0)/2.0/(beta + 1.0)*p_n[1]/p_n[n + 1];
    tmp[n + 1, 1] = (beta + 1.0)/2.0/(alpha + 1.0)*p_n[n + 1]/p_n[1];
    tmp
end

#-------------------------------------------------------
#-----
function solve_NLP(
    ns::Int, 
    nx::Int,
    nu::Int,
    n::Array{Int, 1}, 
    x::Array{Array{Real, 2}, 1}, 
    constraints_x::Array{Tuple{Real, Real}, 1},
    u::Array{Array{Real, 2}, 1}, 
    constraints_u::Array{Tuple{Real, Real}, 1},
    bound_x::Array{Tuple{Real, Real}, 1},
    bound_u::Array{Tuple{Real, Real}, 1},
    constraints_eq::Array{Tuple{Function, Array{Real, 1}}},
    constraints_le::Array{Tuple{Function, Array{Real, 1}}},
    f::Array{Function},
    na::Int,
    a::Array{Real},
    obj_f::Function,
    na_obj::Int,
    obj_a::Array{Real, 1},
    τ::Array{Array{Real, 1}},
    t0::Array{Real, 1},
    tf::Array{Real, 1},
    D::Array{Array{Real, 2}},
    W::Array{Array{Real, 1}}
)::Tuple{Array{Array{Real, 2}, 1}, Array{Array{Real, 2}, 1}, Array{Array{Real, 2}, 1}, Array{Array{Real, 2}, 1}}

    rocket = JuMP.Model(NLopt.Optimizer)
    JuMP.set_optimizer_attribute(rocket, "algorithm", :LD_SLSQP)

    println(constraints_x)
    println(constraints_u)

    ε_bound = 1e-2
    JuMP.@variables(
        rocket, 
        begin
            constraints_x[j][1] + ε_bound ≥ _x[i = 1:ns, j = 1:nx, k = 1:(n[i] + 1)] ≥ constraints_x[j][2] - ε_bound , (start = x[i][j, k])
            constraints_u[j][1] + ε_bound ≥ _u[i = 1:ns, j = 1:nu, k = 1:(n[i] + 1)] ≥ constraints_u[j][2] - ε_bound , (start = u[i][j, k])
        end
    )

    JuMP.@NLexpressions(
        rocket,
        begin
            T_2[i = 1:ns], (tf[i] - t0[i])/2
            t[i = 1:ns, k = 1:(n[i] + 1)], (T_2[i]*τ[n[i]][k] + (tf[i] + t0[i])/2)
            Dx[i = 1:ns, j = 1:nx, k = 1:(n[i] + 1)], sum(D[n[i]][k, l]*_x[i, j, l] for l in 1:(n[i] + 1))
            # wDx[i = 1:ns, j = 1:nx], sum(Dx[i, j, k]*W[n[i]][k] for k in 1:(n[i] + 1))
        end
    )

    f_sym_obj = Symbol("f_0")
    JuMP.register(rocket, f_sym_obj, nx + nu + na_obj + 4, obj_f, autodiff = true)
    obj_expr = JuMP.add_nonlinear_expression(
        rocket, 
        :($(f_sym_obj)(
            $(nx), 
            $(nu), 
            $(na_obj), 
            $(t[ns, n[ns] + 1]), 
            $([_x[ns, j, n[ns] + 1] for j in 1:nx] ...), 
            $([_u[ns, j, n[ns] + 1] for j in 1:nu] ...), 
            $(obj_a...)
        ))
    )
    @NLobjective(rocket, Max, obj_expr)
    
    f_sym = [Symbol("f_$(j)") for (j, ele_f) in enumerate(f)]

    for (j, ele_f) in enumerate(f)
        JuMP.register(rocket, f_sym[j], nx + nu + na + 4, ele_f, autodiff = true)
    end
        
    fx = JuMP.Containers.SparseAxisArray(
        Dict(
            ((i, j, k), JuMP.add_nonlinear_expression(
                rocket, 
                :($(f_sym[j])(
                    $(nx), 
                    $(nu), 
                    $(na), 
                    $(t[i, k]), 
                    $([_x[i, l, k] for l in 1:nx] ...), 
                    $([_u[i, l, k] for l in 1:nu] ...),
                    $(a...)
                ))
            ))
            for i in 1:ns for (j, ele_f) in enumerate(f) for k in 1:(n[i] + 1)
        )
    )
    # println(fx)

    JuMP.@NLconstraint(rocket, con1[i = 1:ns, j = 1:nx, k = 1:(n[i] + 1)], Dx[i, j, k] == T_2[i]*fx[i, j, k])

    if (ns >= 1)
        JuMP.@NLconstraint(rocket,con3[j = 1:nx],  _x[ns, j, n[ns] + 1] - _x[ns, j, 1] == wDx[ns, j])
        if (ns > 1)
            JuMP.@NLconstraint(rocket, con2[i = 2:ns, j = 1:nx], _x[i, j, 1] - _x[i - 1, j, 1] == wDx[i, j])
            JuMP.@NLconstraint(rocket, con4[i = 2:ns, j = 1:nx], _x[i, j, 1] == _x[i - 1, j, n[i - 1] + 1])
            JuMP.@NLconstraint(rocket, con5[i = 2:ns, j = 1:nu], _u[i, j, 1] == _u[i - 1, j, n[i - 1] + 1])
        end
    end

    for (j, ele_f) in enumerate(constraints_eq)
        f_sym_eq  = Symbol("f_eq_$(j)")
        tmp = size(ele_f[2], 1)
        JuMP.register(rocket, f_sym_eq, tmp + 2*nx + 2*nu + 3, ele_f[1], autodiff = true)
        if tmp > 0
            JuMP.add_nonlinear_constraint(
                rocket, 
                :($(f_sym_eq)(
                    $(nx), 
                    $(nu), 
                    $(tmp), 
                    $([_x[1, j, 1] for j in 1:nx] ...), 
                    $([_x[ns, j, n[ns] + 1] for j in 1:nx] ...), 
                    $([_u[1, j, 1] for j in 1:nu] ...), 
                    $([_u[ns, j, n[ns] + 1] for j in 1:nu] ...), 
                    $(ele_f[2]...)
                ) == 0)
            )
        else
            JuMP.add_nonlinear_constraint(
                rocket, 
                :($(f_sym_eq)(
                    $(nx), 
                    $(nu), 
                    $(tmp), 
                    $([_x[1, j, 1] for j in 1:nx] ...), 
                    $([_x[ns, j, n[ns] + 1] for j in 1:nx] ...), 
                    $([_u[1, j, 1] for j in 1:nu] ...), 
                    $([_u[ns, j, n[ns] + 1] for j in 1:nu] ...)
                ) == 0)
            )
        end
    end

    for (j, ele_f) in enumerate(constraints_le)
        f_sym_le  = Symbol("f_le_$(j)")
        tmp = size(ele_f[2], 1)
        JuMP.register(rocket, f_sym_le, tmp + 2*nx + 2*nu + 3, ele_f[1], autodiff = true)
        if tmp > 0
            
            JuMP.add_nonlinear_constraint(
                rocket, 
                :($(f_sym_le)(
                    $(nx), 
                    $(nu), 
                    $(tmp), 
                    $([_x[1, j, 1] for j in 1:nx] ...), 
                    $([_x[ns, j, n[ns] + 1] for j in 1:nx] ...), 
                    $([_u[1, j, 1] for j in 1:nu] ...), 
                    $([_u[ns, j, n[ns] + 1] for j in 1:nu] ...),
                    $(ele_f[2]...)
                ) <= 0)
            )
            # JuMP.add_nonlinear_constraint(rocket, :($(f_sym_le)($(nx), $(nu), $(tmp), $(_x[0, :, 0]...), $(_x[ns, :, n[ns] + 1]...), $(_u[0, :, 0]...), $(_u[ns, :, n[ns] + 1]...), $(ele_f[2]...)) <= 0))
        else
            JuMP.add_nonlinear_constraint(
                rocket, 
                :($(f_sym_le)(
                    $(nx), 
                    $(nu), 
                    $(tmp), 
                    $([_x[1, j, 1] for j in 1:nx] ...), 
                    $([_x[ns, j, n[ns] + 1] for j in 1:nx] ...), 
                    $([_u[1, j, 1] for j in 1:nu] ...), 
                    $([_u[ns, j, n[ns] + 1] for j in 1:nu] ...)
                ) <= 0)
            )
            # JuMP.add_nonlinear_constraint(rocket, :($(f_sym_le)($(nx), $(nu), $(tmp), $(_x[1, :, 1]...), $(_x[ns, :, n[ns] + 1]...), $(_u[1, :, 1]...), $(_u[ns, :, n[ns] + 1]...)) <= 0))
        end
    end

    for (i, ele) in enumerate(bound_x)
        if !isnan(ele[1]) 
            JuMP.fix(_x[1, i, 1], ele[1]; force = true)
        end
        if !isnan(ele[2]) 
            JuMP.fix(_x[ns, i, n[ns] + 1], ele[2]; force = true)
        end
    end

    for (i, ele) in enumerate(bound_u)
        if !isnan(ele[1]) 
            JuMP.fix(_u[1, i, 1], ele[1]; force = true)
        end
        if !isnan(ele[2]) 
            JuMP.fix(_u[ns, i, n[ns] + 1], ele[2]; force = true)
        end
    end

    println("Solving...")
    JuMP.optimize!(rocket)
    print(solution_summary(rocket))

    rx = JuMP.value.(_x)
    ru = JuMP.value.(_u)
    rDx = JuMP.value.(Dx)
    rfx = JuMP.value.(fx)

    _rx = [[rx[i, j, k] for j in 1:nx, k in 1:(n[i] + 1)] for i in 1:ns]
    _ru = [[ru[i, j, k] for j in 1:nu, k in 1:(n[i] + 1)] for i in 1:ns]
    _rDx = [[rDx[i, j, k] for j in 1:nx, k in 1:(n[i] + 1)] for i in 1:ns]
    _rfx = [[rfx[i, j, k] for j in 1:nx, k in 1:(n[i] + 1)] for i in 1:ns]

    delete_x = JuMP.all_variables(rocket)
    JuMP.delete(rocket, delete_x)
    JuMP.unregister(rocket, :delete_x)
    
    delete_con = JuMP.all_nonlinear_constraints(rocket)
    for ele in delete_con
        JuMP.delete(rocket, ele)
        JuMP.unregister(rocket, :ele)
    end
    
    return (
        _rx,
        _ru,
        _rDx,
        _rfx,
    )
end

function main()

    n_max = 50
    n_min = 1
    τ = Array{Real, 1}[compute_tau(n) for n in n_min:n_max]
    W = Array{Real, 1}[compute_W(n, -0.5, -0.5, τ[n]) for n in n_min:n_max]
    D = Array{Real, 2}[compute_D(n, -0.5, -0.5, τ[n]) for n in n_min:n_max]

    ns = 1
    nx = 4
    nu = 1
    na = 2
    # n = Int[7 for i in 1:ns]
    n = Int[12 for i in 1:ns]
    
    x = [
        Real[(j == 1 || j == 4) ? 1.0 : 0.0 for j in 1:nx, k in 1:(n[i] + 1)] for i in 1:ns
    ]

    cx = Tuple{Real, Real}[(10.0, 0.0) for j in 1:nx]
    u = [Real[0.001 for j in 1:nu, k in 1:(n[i] + 1)] for i in 1:ns]
    cu = Tuple{Real, Real}[(10.0, 0.0) for j in 1:nu]
    bx = Tuple{Real, Real}[
        (1.0, NaN),
        (0.0, NaN),
        (0.0, 0.0),
        (1.0, NaN)
    ]
    bu = Tuple{Real, Real}[
        (0.001, NaN)
    ]
    c_eq = Tuple{Function, Array{Real, 1}}[
        (
            ((nx, nu, na, r0, θ0, v_r0, v_θ0, rf, θf, v_rf, v_θf, γ0, γf) -> v_θf - 1/sqrt(rf)), []
        )
    ]
    c_le = Tuple{Function, Array{Real, 1}}[]
    f = Function[
        ((nx, nu, na, t, r, θ, v_r, v_θ, γ, T, m_p) -> (v_r)),
        ((nx, nu, na, t, r, θ, v_r, v_θ, γ, T, m_p) -> (v_θ/r)),
        ((nx, nu, na, t, r, θ, v_r, v_θ, γ, T, m_p) -> (v_θ^2.0/r - r^(-2.0) + T/(1 - m_p*t)*sin(γ))),
        ((nx, nu, na, t, r, θ, v_r, v_θ, γ, T, m_p) -> (-v_r*v_θ/r + T/(1 - m_p*t)*cos(γ)))
    ]
    a = Real[0.1405, 0.0749]
    f_obj = (nx, nu, na, t, x...) -> x[1]
    na_obj = 0
    a_obj = Real[]
    t0 = Real[0.0]
    tf = Real[3.32]

    ε_tol = 2

    while (true)
        (x, u, dx, fx) = solve_NLP(ns, nx, nu, na, n, x, cx, u, cu, bx, bu, c_eq, c_le, f, a, f_obj, na_obj, a_obj, τ, t0, tf, D, W)
        
        (ns, n, x, u, t0, tf) = refine(ns, nx, nu, n, x, u, dx, fx, τ, t0, tf, ε_tol, 0.8, W)
        if ns < 0
            ns *= -1
            break
        end
    end
    
    # println(fx[1][1, :])
    # println(dx[1][1, :])
    # println([dx[1][1, i] - (tf[1]-t0[1])/2*fx[1][1, i] for i in 1:(n[1] + 1)])

    plot_x = cat(x..., dims = 2)
    plot_u = cat(u..., dims = 2)
    plot_t = cat([((τ[n[i]]).*((tf[i] - t0[i])/2)).+((tf[i] + t0[i])/2) for i in 1:ns]..., dims = 1)

    minmax_plot_x = ((x) -> (min(x...), max(x...)))(cat(plot_x..., dims = 1))
    minmax_plot_u = ((x) -> (min(x...), max(x...)))(cat(plot_u..., dims = 1))

    open("file.txt","w") do out
        println(out, "idx,t,x1,x2,x3,x4,u1")
        for (i, ele_t) in enumerate(plot_t)
            println(
                out,
                i,
                ",",
                ele_t,
                ",",
                plot_x[1, i],
                ",",
                plot_x[2, i],
                ",",
                plot_x[3, i],
                ",",
                plot_x[4, i],
                ",",
                plot_u[1, i]
            )
        end
    end

    plots_x = Plots.plot(
        plot_t[:], 
        [ plot_x[1, :] plot_x[2, :] plot_x[3, :] plot_x[4, :] ], 
        label=["r" "θ" "u" "v"],
        xlabel = "t",
        st=:scatter
    )
    for ti in cat(t0, [tf[end]], dims = 1)
        plot!(
            plots_x,
            [ti, ti],
            [minmax_plot_x[1], minmax_plot_x[2]],
            label="",
            ls=:dash,
            lc=:gray
        )
    end

    plots_u = Plots.plot(
        plot_t, 
        plot_u[1, :], 
        label="γ",
        xlabel = "t",
        st=:scatter
    )
    for ti in cat(t0, [tf[end]], dims = 1)
        plot!(
            plots_u,
            [ti, ti],
            [minmax_plot_u[1], minmax_plot_u[2]],
            label="",
            ls=:dash,
            lc=:gray
        )
    end
    
    Plots.plot(
        plots_x,
        plots_u,
        layout = (2, 1),
        # legend = false,
        # margin = 1Plots.cm,
    )
end

main()