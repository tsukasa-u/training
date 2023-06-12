using JuMP
import Ipopt
# using NLopt
using Plots

using Polynomials
using SpecialPolynomials
using SpecialFunctions

using ForwardDiff

using Roots

index = 0

#-------------------------------------------------------

function compute_η(
        n::Int, 
        dx::Array{Real, 1}, 
        f::Array{Real, 1}, 
        w::Array{Real, 1}
    )::Real
        return sum(abs(dx[k] - f[k])*w[k] for k in 1:n)
end

function compute_w(
        ns::Int, 
        n::Array{Int, 1}, 
        x::Array{Array{Real, 1}, 1},
        dx::Array{Array{Real, 1}, 1}
    )::Real
        # return maximum([maximum([max(abs(x[i][k]), abs(dx[i][k])) for k in 1:(n[i] + 1)]) for i in 1:ns])
        return maximum([
            max(
                maximum([
                    abs(x[i][k]) 
                    for k in 1:(n[i] + 1)
                ]), 
                maximum([
                    abs(dx[i][k])
                    for k in 1:n[i]
                ])
            ) 
            for i in 1:ns
        ])
end

function compute_ε(
        ns::Int, 
        n::Array{Int, 1}, 
        dx::Array{Array{Real, 1}, 1}, 
        x::Array{Array{Real, 1}, 1}, 
        f::Array{Array{Real, 1}, 1}, 
        W::Array{Array{Real, 1}, 1}
    )::Tuple{Array{Real, 1}, Real}
        w::Real = compute_w(ns, n, x, dx)
        return ([compute_η(n[i], dx[i], f[i], W[n[i]])/w for i in 1:ns], w)
end

function check_ε_tolerance(ε_tol::Real, ns::Real, ε::Array{Array{Real, 1}, 1})::Real
    tmp = maximum([sum(ele[i] for i in ns) for ele in ε])
    return tmp < ε_tol ? tmp : -1
end

function select_refine_segment(ε_tol::Real, α::Real, nx::Int, ε::Array{Array{Real, 1}, 1})::Array{Int, 1}
    ret_index::Array{Int, 1} = []
    for (j, ele) in enumerate(ε)
        idx = sortperm(ele)
        println(idx)
        sum = 0
        # r = 0
        # push!(ret_index, idx[end])
        union!(ret_index, idx[end])
        for i in idx
            sum += ele[i]
            if sum > α*ε_tol
                # push!(ret_index, i)
                union!(ret_index, i)
                # ret_index[j] = r
                # break
            end
        end
    end
    return ret_index
end

function compute_difficulty(n::Int, f::Array{Real, 2}, w::Array{Real, 1}, ε::Array{Real, 1}, τ::Array{Real, 1}, t0::Real, tf::Real)::Array{Real, 1}

    F_k = [sum(ε[j]*f[j, k]/w[j] for (j, ele) in enumerate(ε)) for k in 1:(n+1)]
    # F_1 = [(F_k[k+1] - F_k[k])*2/(tf - t0)/(τ[k+1] - τ[k]) for k in 1:n]
    # insert!(F_1, 1, NaN)
    F_2 = [(F_k[k] - F_k[k-1])*2/(tf - t0)/(τ[k+1] - τ[k-1]) for k in 2:n]
    insert!(F_2, 1, NaN)

    sum_F = sum(abs(F_2[l]) for l in 2:n)
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
    
    xm = Array{Real, 2}[]
    um = Array{Real, 2}[]
    tk_size = size(tk, 1)
    for tk_idx in 2:tk_size
        _tk_idx = tk_idx - 1
        _tk = tk[_tk_idx]
        tk_ = tk[tk_idx]
        tau = τn[_tk:tk_]

        if size(tau, 1) > 0

            tau = tau.*(2.0/(tau[end] - tau[1]))
            tau = tau .- (1 + tau[1])

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
                    # _um[j, k] = un[j, _tk + i]*ratio_b - un[j, _tk + i - 1]*ratio_a
                    _um[j, k] = un[j, _tk + i]
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
        W::Array{Array{Real, 1}, 1};
        index = -1
    )::Tuple{Int , Array{Int, 1}, Array{Array{Real, 2}, 1}, Array{Array{Real, 2}, 1}, Array{Real, 1}, Array{Real, 1}, Array{Array{Real, 1}, 1}, Array{Real, 1}}
        # εx::Array{Array{Real, 1}, 1} = Array{Real, 1}[]
        # wx::Array{Real, 1} = Real[]
        εx = Array{Real, 1}[]
        wx = Real[]
        t = Array{Real, 1}[[(tf[i] - t0[i])/2*τ[n[i]][k] + (tf[i] + t0[i])/2 for k in 1:(n[i] + 1)] for i in 1:ns]
        for j in 1:nx
            tmp = compute_ε(
                ns, 
                n, 
                Array{Real, 1}[[dx[i][j, k] for k in 1:n[i]] for i in 1:ns], 
                Array{Real, 1}[[x[i][j, k] for k in 1:(n[i] + 1)] for i in 1:ns], 
                Array{Real, 1}[[fx[i][j, k] for k in 1:(n[i] + 1)] for i in 1:ns], 
                W,
            )
            push!(εx, tmp[1])
            push!(wx, tmp[2])
        end

        option(index, ns, n, x, u, τ, tf, t0, ε = (true, nx, εx, wx, fx))
        for j in 1:nx
            for i in 1:ns
                print(" ", εx[j][i])
            end
            println("")
        end

        if check_ε_tolerance(ε_tol, ns, εx) > 0
            return (-ns, n, x, u, t0, tf, εx, wx)
        end

        x_idx = select_refine_segment(ε_tol, α, nx, εx)
        println(x_idx)

        rx = copy(x)
        ru = copy(u)
        rn = copy(n)
        rt0 = copy(t0)
        rtf = copy(tf)
        offset = 0

        for idx in unique(sort(x_idx))
            deleteat!(rx, idx + offset)
            deleteat!(ru, idx + offset)
            deleteat!(rn, idx + offset)
            deleteat!(rt0, idx + offset)
            deleteat!(rtf, idx + offset)
            
            (_len, _m, _t0, _tf, segx, segu) = refine_segment(n[idx], nx, nu, x[idx], u[idx], fx[idx], wx, Real[εx[j][idx] for j in 1:nx], τ, t0[idx], tf[idx])
            
            for i in _len:-1:1
                insert!(rx, idx + offset, segx[i])
                insert!(ru, idx + offset, segu[i])
                insert!(rn, idx + offset, _m[i])
                insert!(rt0, idx + offset, _t0[i])
                insert!(rtf, idx + offset, _tf[i])
            end

            offset = offset + _len - 1
        end

        return (ns + offset, rn, rx, ru, rt0, rtf, εx, wx)  # copied or returned pointer?
end

function fix_constraints(
    ns::Int,
    nx::Int,
    nu::Int,
    n::Array{Int, 1},
    x::Array{Array{Real, 2}, 1}, 
    u::Array{Array{Real, 2}, 1}, 
    constraints_x::Array{Tuple{Real, Real}, 1},
    constraints_u::Array{Tuple{Real, Real}, 1}
)::Tuple{Array{Array{Real, 2}, 1}, Array{Array{Real, 2}, 1}}

    rx = copy(x)
    ru = copy(u)

    for j in 1:nx
        x_u = constraints_x[j][1]
        x_d = constraints_x[j][2]
        for i in 1:ns
            for k in 1:(n[i] + 1)
                if x[i][j, k] >= x_u
                    rx[i][j, k] = x_u
                end
                if x[i][j, k] <= x_d
                    rx[i][j, k] = x_d
                end
            end
        end
    end

    for j in 1:nu
        u_u = constraints_u[j][1]
        u_d = constraints_u[j][2]
        for i in 1:ns
            for k in 1:(n[i] + 1)
                if u[i][j, k] >= u_u
                    ru[i][j, k] = u_u
                end
                if u[i][j, k] <= u_d
                    ru[i][j, k] = u_d
                end
            end
        end
    end

    return (rx, ru)
end

function Real_binomial(n::Real, k::Real)
    return gamma(1 + n)/gamma(1 + k)/gamma(1 + n - k)
end

function compute_tau(n, alpha, beta, max_steps = 500, tolerance = 1E-12)

    if n == 1
        return [-1.0, 1.0]
    end

    # k1 = 1
    # P1 = Jacobi{alpha + k1, beta + k1}
    # Pn1 = gamma(alpha + beta + n + 1 + k1)/2^k1/gamma(alpha + beta + n + 1)*basis(P1, n - k1)

    # k2 = 2
    # P1 = Jacobi{alpha + k2, beta + k2}
    # Pn2 = gamma(alpha + beta + n + 1 + k2)/2^k1/gamma(alpha + beta + n + 1)*basis(P1, n - k2)

    # inix = updatex = -1.0
    # array_x = Real[-1.0]
    # for i in -(max_steps - 1):1:(max_steps -1)
    #     inix = i/max_steps
    #     for j in 1:50

    #         updatex =  inix - Pn1(inix) / Pn2(inix)
    #         if abs(updatex - inix) < tolerance
    #             if ((i - 0.5)/max_steps <= updatex && (i + 0.5)/max_steps > updatex)
    #                 push!(array_x, updatex)
    #                 sort!(array_x)
    #             end
    #             break
    #         end
    #         inix = updatex
    #     end
    # end
    # push!(array_x, 1.0)
    # return array_x
    P_n = basis(Jacobi{alpha, beta}, n)
    dPn_dτ(x) = ForwardDiff.derivative(P_n, x)
    ret = cat([-1.0], find_zeros(dPn_dτ, -1.0,  1.0)..., [1.0]; dims = 1)
    @assert (length(ret) == n + 1) "$(n + 1) op $(length(ret))"
    return ret
end

function base_L(x, k, τ::Array{Real, 1})
    a1 = cat(τ[1:(k - 1)], τ[(k + 1):end]; dims = 1)
    tk = τ[k]
    b1 = (x .- a1)./(tk .- a1)
    return reduce(*, b1)
end

dL_dτ(x, k, τ::Array{Real, 1}) = ForwardDiff.derivative(x -> base_L(x, k, τ), x)

# function dL_dτ(n, i, k, τ::Array{Real, 1})
#     if n == 1
#         return 1.0
#     end
#     taui = copy(τ)
#     deleteat!(taui, i)
#     ret = 0.0
#     for l in 1:(n + 1)
#         tauk = copy(τ)
#         deleteat!(tauk, unique(sort([i, l])))
#         ret += reduce(*, τ[k] .- tauk)
#     end
#     return ret/reduce(*, τ[i] .- taui)
#     # return sum(reduce(*, deleteat!(τ, [k, l])) for l in 1:(n + 1))/bi
# end

# function compute_W(n, τ::Array{Real, 1})
    
# end

# function compute_W(n, α, β, τ::Array{Real, 1})
    
#     P = Jacobi{α, β}
#     Pn = basis(P, n)
#     p_n = Pn.(τ)

#     # fill(pi/n, n)
#     _w = fill(0.0, n + 1)

#     a_1 = 2^(α + β + 1)*gamma(α + 2)*gamma(β + 2)/gamma(α + β + 3)/n^2*Real_binomial(n + α, n - 1)*Real_binomial(n + β, n - 1)/Real_binomial(n + α + β + 1, n - 1)
#     a_2 = (4*(n + α)*(n + β) + (α - β)^2)/(2*n + α + β)^2

#     for i in 2:n
#         _w[i] = a_1/(a_2 - τ[i]^2)*(1 - τ[i]^2)/p_n[i]^2
#     end
#     _w[1] = 2^(α + β + 1)*gamma(α + 2)*gamma(β + 1)/gamma(α + β + 3)*Real_binomial(n + α, n - 1)/Real_binomial(n + β, n - 1)/Real_binomial(n + α + β + 1, n - 1)
#     _w[n + 1] = _w[1]
#     return _w
# end

function compute_W(n, D::Array{Real, 2})
    return vec(inv(D[:, 2:end])[end, :])
end

# function compute_W(n, α, β, τ)
#     return [(1 - x)^α*(1 + x)^β for x in τ]
# end

function compute_D(n::Int, tau::Array{Real, 1})
    
    D = zeros(n, n + 1)
    for i in 1:n
        for j in 1:(n + 1)
            D[i, j] = dL_dτ(tau[i], j, tau)
        end
    end
    return D
    # return [dL_dτ(n, i, j, tau) for i in 1:n, j in 1:(n + 1)]
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
    bound_eq::Array{Tuple{Function, Array{Real, 1}}},
    bound_le::Array{Tuple{Function, Array{Real, 1}}},
    f::Array{Function},
    na::Int,
    a::Array{Real},
    obj_f::Function,
    na_obj::Int,
    obj_a::Array{Real, 1},
    int_f::Function,
    na_int::Int,
    int_a::Array{Real},
    τ::Array{Array{Real, 1}},
    t0::Array{Real, 1},
    tf::Array{Real, 1},
    constraints_t::Tuple{Real, Real},
    bound_t::Tuple{Real, Real},
    D::Array{Array{Real, 2}},
    W::Array{Array{Real, 1}},
    ε_bound::Real
)::Tuple{Array{Array{Real, 2}, 1}, Array{Array{Real, 2}, 1}, Array{Array{Real, 2}, 1}, Array{Array{Real, 2}, 1}, Array{Real, 1}, Array{Real, 1}}

    # rocket = JuMP.Model(NLopt.Optimizer)
    # JuMP.set_attribute(rocket, "algorithm", :LD_SLSQP)
    rocket = JuMP.Model(Ipopt.Optimizer)
    JuMP.set_attribute(rocket, "tol", 1e-12)

    t0f = cat(t0, [tf[end]], dims = 1)

    # ε_bound = 1e-2
    JuMP.@variables(
        rocket, 
        begin
            constraints_t[1] + ε_bound ≥ _t0f[i = 1:(ns + 1)] ≥ constraints_t[2] - ε_bound, (start = t0f[i])
            constraints_x[j][1] + ε_bound ≥ _x[i = 1:ns, j = 1:nx, k = 1:(n[i] + 1)] ≥ constraints_x[j][2] - ε_bound , (start = x[i][j, k])
            constraints_u[j][1] + ε_bound ≥ _u[i = 1:ns, j = 1:nu, k = 1:(n[i] + 1)] ≥ constraints_u[j][2] - ε_bound , (start = u[i][j, k])
        end
    )

    JuMP.@NLconstraint(rocket, [i = 1:ns], _t0f[i + 1] >= _t0f[i] + ε_bound)
    for i in 2:ns
        JuMP.fix(_t0f[i], t0f[i]; force = true)
    end
    if !isnan(bound_t[1])
        JuMP.fix(_t0f[1], bound_t[1]; force = true)
    end
    if !isnan(bound_t[2])
        JuMP.fix(_t0f[ns + 1], bound_t[2]; force = true)
    end

    JuMP.@NLexpressions(
        rocket,
        begin
            T_2[i = 1:ns], (_t0f[i + 1] - _t0f[i])/2
            t[i = 1:ns, k = 1:(n[i] + 1)], (T_2[i]*τ[n[i]][k] + (_t0f[i + 1] + _t0f[i])/2)
            Dx[i = 1:ns, j = 1:nx, k = 1:n[i]], sum(D[n[i]][k, l]*_x[i, j, l] for l in 1:(n[i] + 1))
            # xDw[i = 1:ns, j = 1:nx], sum(Dx[i, j, k]*W[n[i]][k] for k in 1:(n[i] + 1))
            # xDw[i = 1:ns, j = 1:nx], sum(Dx[i, j, k]*W[n[i]][k] for k in 1:n[i])
            # xDw[i = 1:ns, j = 1:nx], sum(_x[i, j, l]*sum(D[n[i]][k, l]*W[n[i]][k] for k in 1:(n[i] + 1)) for l in 1:(n[i] + 1))
        end
    )
    
    f_sym_obj = Symbol("f_obj")
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

    f_sym_int = Symbol("f_int")
    JuMP.register(rocket, f_sym_int, nx + nu + na_int + 4, int_f, autodiff = true)
    int_expr = JuMP.Containers.SparseAxisArray(
        Dict(
            ((i, k), JuMP.add_nonlinear_expression(
                rocket, 
                :($(f_sym_int)(
                    $(nx), 
                    $(nu), 
                    $(na_int), 
                    $(t[i, k]), 
                    $([_x[i, j, k] for j in 1:nx] ...), 
                    $([_u[i, j, k] for j in 1:nu] ...), 
                    $(int_a...)
                ))
            ))
            for i in 1:ns for k in 1:(n[i] + 1)
        )
    )

    # JuMP.@NLexpression(rocket, wf, sum(W[n[i]][k]*int_expr[i, k] for i in 1:ns for k in 1:(n[i] + 1)))
    # JuMP.@NLexpression(rocket, wf, sum(T_2[i]*sum(W[n[i]][k]*int_expr[i, k] for k in 1:(n[i] + 1)) for i in 1:ns ))
    # JuMP.@NLexpression(rocket, wf, sum(T_2[i]*W[n[i]][k]*int_expr[i, k] for i in 1:ns for k in 1:(n[i] + 1)))
    JuMP.@NLexpression(rocket, wf, sum(T_2[i]*W[n[i]][k]*int_expr[i, k] for i in 1:ns for k in 1:n[i]))
    JuMP.@NLexpression(rocket, wur, sum(((_u[i, 1, k] - _u[i, 2, k]) - (_u[i, 1, k + 1] - _u[i, 2, k + 1]))^2 for i in 1:ns for k in 1:n[i]))
    JuMP.@NLexpression(rocket, wut, sum(((_u[i, 3, k] - _u[i, 4, k]) - (_u[i, 3, k + 1] - _u[i, 4, k + 1]))^2 for i in 1:ns for k in 1:n[i]))

    @NLobjective(rocket, Min, obj_expr + wf + wur + wut)
    
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

    # JuMP.@NLconstraint(rocket, [i = 1:ns, j = 1:nx, k = 1:(n[i] + 1)], Dx[i, j, k] == T_2[i]*fx[i, j, k])
    JuMP.@NLconstraint(rocket, [i = 1:ns, j = 1:nx, k = 1:n[i]], Dx[i, j, k] == T_2[i]*fx[i, j, k])

    if (ns >= 1)
        # JuMP.@NLconstraint(rocket,[j = 1:nx],  _x[ns, j, n[ns] + 1] - _x[ns, j, 1] == xDw[ns, j])
        if (ns > 1)
            # JuMP.@NLconstraint(rocket, [i = 2:ns, j = 1:nx], _x[i, j, 1] - _x[i - 1, j, 1] == xDw[i, j])
            JuMP.@NLconstraint(rocket, [i = 2:ns, j = 1:nx], _x[i, j, 1] == _x[i - 1, j, n[i - 1] + 1])
            # JuMP.@NLconstraint(rocket, [i = 2:ns, j = 1:nx], Dx[i, j, 1] == Dx[i - 1, j, n[i - 1]])

            JuMP.@NLconstraint(rocket, [i = 2:ns, j = 1:nu], _u[i, j, 1] == _u[i - 1, j, n[i - 1] + 1])
            JuMP.@NLconstraint(rocket, [i = 2:ns, j = 1:nu], _u[i, j, 1] == _u[i - 1, j, n[i - 1]])
        end
    end

    for (j, ele_f) in enumerate(constraints_eq)
        f_sym_eq  = Symbol("f_eq_$(j)")
        tmp = size(ele_f[2], 1)
        JuMP.register(rocket, f_sym_eq, tmp + nx + nu + 4, ele_f[1], autodiff = true)
        for i in 1:ns
            for k in 1:(n[i] + 1)
                JuMP.add_nonlinear_constraint(
                    rocket, 
                    :($(f_sym_eq)(
                        $(nx), 
                        $(nu), 
                        $(tmp), 
                        $(t[i, k]), 
                        $([_x[i, j, k] for j in 1:nx] ...), 
                        $([_u[i, j, k] for j in 1:nu] ...), 
                        $(ele_f[2]...)
                    ) == 0)
                )
            end
        end
    end

    for (j, ele_f) in enumerate(constraints_le)
        f_sym_le  = Symbol("f_le_$(j)")
        tmp = size(ele_f[2], 1)
        JuMP.register(rocket, f_sym_le, tmp + nx + nu + 4, ele_f[1], autodiff = true)
        for i in 1:ns
            for k in 1:(n[i] + 1)
                JuMP.add_nonlinear_constraint(
                    rocket, 
                    :($(f_sym_le)(
                        $(nx), 
                        $(nu), 
                        $(tmp), 
                        $(t[i, k]), 
                        $([_x[i, j, k] for j in 1:nx] ...), 
                        $([_u[i, j, k] for j in 1:nu] ...), 
                        $(ele_f[2]...)
                    ) <= 0)
                )
            end
        end
    end

    for (j, ele_f) in enumerate(bound_eq)
        f_sym_eq  = Symbol("f_eq_$(j)")
        tmp = size(ele_f[2], 1)
        JuMP.register(rocket, f_sym_eq, tmp + 2*nx + 2*nu + 5, ele_f[1], autodiff = true)
        JuMP.add_nonlinear_constraint(
            rocket, 
            :($(f_sym_eq)(
                $(nx), 
                $(nu), 
                $(tmp), 
                $(_t0f[1]),
                $(_t0f[ns + 1]),
                $([_x[1, j, 1] for j in 1:nx] ...), 
                $([_x[ns, j, n[ns] + 1] for j in 1:nx] ...), 
                $([_u[1, j, 1] for j in 1:nu] ...), 
                $([_u[ns, j, n[ns] + 1] for j in 1:nu] ...), 
                $(ele_f[2]...)
            ) == 0)
        )
    end

    for (j, ele_f) in enumerate(bound_le)
        f_sym_le  = Symbol("f_le_$(j)")
        tmp = size(ele_f[2], 1)
        JuMP.register(rocket, f_sym_le, tmp + 2*nx + 2*nu + 5, ele_f[1], autodiff = true)
        JuMP.add_nonlinear_constraint(
            rocket, 
            :($(f_sym_le)(
                $(nx), 
                $(nu), 
                $(tmp),
                $(_t0f[1]),
                $(_t0f[ns + 1]),
                $([_x[1, j, 1] for j in 1:nx] ...), 
                $([_x[ns, j, n[ns] + 1] for j in 1:nx] ...), 
                $([_u[1, j, 1] for j in 1:nu] ...), 
                $([_u[ns, j, n[ns] + 1] for j in 1:nu] ...),
                $(ele_f[2]...)
            ) <= 0)
        )
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
    # print(solution_summary(rocket))
    println(JuMP.objective_value(rocket))

    rx = JuMP.value.(_x)
    ru = JuMP.value.(_u)
    rDx = JuMP.value.(Dx)
    rfx = JuMP.value.(fx)
    rt0f = JuMP.value.(_t0f)

    _rx = [[rx[i, j, k] for j in 1:nx, k in 1:(n[i] + 1)] for i in 1:ns]
    _ru = [[ru[i, j, k] for j in 1:nu, k in 1:(n[i] + 1)] for i in 1:ns]
    _rDx = [[rDx[i, j, k] for j in 1:nx, k in 1:n[i]] for i in 1:ns]
    _rfx = [[rfx[i, j, k] for j in 1:nx, k in 1:(n[i] + 1)] for i in 1:ns]
    _rt0 = [rt0f[i] for i in 1:ns]
    _rtf = [rt0f[i] for i in 2:(ns + 1)]

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
        _rt0,
        _rtf
    )
end
function main()

    n_max = 110
    n_min = 5
    alpha = 0.0
    beta = 0.0
    τ = Array{Real, 1}[n >= n_min ? compute_tau(n, alpha, beta) : [] for n in 1:n_max]
    D = Array{Real, 2}[n >= n_min ? compute_D(n, τ[n]) : [0 0; 0 0;] for n in 1:n_max]
    W = Array{Real, 1}[n >= n_min ? compute_W(n, D[n]) : [] for n in 1:n_max]
    # W = Array{Real, 1}[n >= n_min ? compute_W(n, alpha, beta, τ[n]) : [] for n in 1:n_max]

    ns = 1
    nx = 3
    nu = 4
    na = 0
    n = Int[100 for i in 1:ns]
    
    x = [
        Real[
            if j == 1
                3.0*(k - 1)/n[i] + 1.0
            elseif j == 2
                0.0
            else
                -0.5*(k - 1)/n[i] + 1.0
            end
            for j in 1:nx, k in 1:(n[i] + 1)
        ] for i in 1:ns
    ]

    cx = Tuple{Real, Real}[(10.0, 0.0) for j in 1:nx]
    u = [Real[0.001 for j in 1:nu, k in 1:(n[i] + 1)] for i in 1:ns]
    cu = Tuple{Real, Real}[(0.01, 0.0) for j in 1:nu]
    ct = (55.0, 0.0)
    bt = (0.0, NaN)
    bx = Tuple{Real, Real}[
        (1.0, 4.0),
        (0.0, 0.0),
        (1.0, 0.5)
    ]
    bu = Tuple{Real, Real}[
        (NaN, NaN),
        (NaN, NaN),
        (NaN, NaN),
        (NaN, NaN)
    ]
    b_eq = Tuple{Function, Array{Real, 1}}[]
    b_le = Tuple{Function, Array{Real, 1}}[]
    c_eq = Tuple{Function, Array{Real, 1}}[]
    c_le = Tuple{Function, Array{Real, 1}}[
        # ((nx, nu, na, t, r, v_r, v_t, u_rp, u_rm, u_tp, u_tm) -> u_rp + u_rm - 0.01, []),
        # ((nx, nu, na, t, r, v_r, v_t, u_rp, u_rm, u_tp, u_tm) -> u_tp + u_tm - 0.01, []),
        # ((nx, nu, na, t, r, v_r, v_t, u_rp, u_rm, u_tp, u_tm) -> - u_rp * u_rm, []),
        # ((nx, nu, na, t, r, v_r, v_t, u_rp, u_rm, u_tp, u_tm) -> - u_tp * u_tm, [])
    ]
    f = Function[
        ((nx, nu, na, t, r, v_r, v_t, u_rp, u_rm, u_tp, u_tm) -> (v_r)),
        ((nx, nu, na, t, r, v_r, v_t, u_rp, u_rm, u_tp, u_tm) -> (v_t^2.0/r - 1.0/r^2.0 + u_rp - u_rm)),
        ((nx, nu, na, t, r, v_r, v_t, u_rp, u_rm, u_tp, u_tm) -> (-v_r*v_t/r + u_tp - u_tm))
    ]
    a = Real[]
    f_obj = (nx, nu, na, t, r, v_r, v_t, u_rp, u_rm, u_tp, u_tm) -> 0.0
    na_obj = 0
    a_obj = Real[]
    f_int = (nx, nu, na, t, r, v_r, v_t, u_rp, u_rm, u_tp, u_tm) -> u_rp + u_rm + u_tp + u_tm
    na_int = 0
    a_int = Real[]
    t0 = Real[0.0]
    tf = Real[55.0]

    ε_tol = 0.02

    while (true)
    # for _ in 1:1
        global index += 1
        
        (x, u) = fix_constraints(ns, nx, nu, n, x, u, cx, cu)
        
        option(index, ns > 0 ? ns : -ns, n, x, u, τ, tf, t0)

        (x, u, dx, fx, t0, tf) = solve_NLP(ns, nx, nu, n, x, cx, u, cu, bx, bu, c_eq, c_le, b_eq, b_le, f, na, a, f_obj, na_obj, a_obj, f_int, na_int, a_int, τ, t0, tf, ct, bt, D, W, 1e-12)
        global index += 1
        # option(index, ns > 0 ? ns : -ns, n, x, u, τ, tf, t0)
        
        (ns, n, x, u, t0, tf, εx, wx) = refine(ns, nx, nu, n, x, u, dx, fx, τ, t0, tf, ε_tol, 0.8, W, index = index)
        
        if ns < 0
            ns *= -1
            global index += 1
            option(index, ns, n, x, u, τ, tf, t0)
            break
        end
    end
    
end

function option(index, ns, n, x, u, τ, tf, t0; ε = (false, 0, NaN, NaN, NaN))

    plot_x = cat(x..., dims = 2)
    plot_u = cat(u..., dims = 2)
    plot_t = cat([((τ[n[i]]).*((tf[i] - t0[i])/2)).+((tf[i] + t0[i])/2) for i in 1:ns]..., dims = 1)
    (isε, nx, εx, wx, fx) = ε
    plot_ε = []
    if isε
        plot_ε = cat([compute_difficulty(n[i], fx[i], wx, Real[εx[j][i] for j in 1:nx], τ[n[i]], t0[i], tf[i]) for i in 1:ns]..., dims = 1)
    end

    # save_file(plot_x, plot_u, plot_t)
    plot_graph(index, plot_x, plot_u, plot_ε, plot_t, tf, t0)
end

function save_file(plot_x, plot_u, plot_t)

    open("file.txt","w") do out
        println(out, "idx,t,x1,x2,x3,u1,u2")
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
                plot_u[1, i],
                ",",
                plot_u[2, i]
            )
        end
    end
end

function plot_graph(index, plot_x, plot_u, plot_ε, plot_t, tf, t0)
    minmax_plot_x = ((x) -> (min(x...), max(x...)))(cat(plot_x..., dims = 1))
    # minmax_plot_u = ((x) -> (min(x...), max(x...)))(cat(plot_u..., dims = 1))
    minmax_plot_u = (-0.01, 0.01)

    plots_x = Plots.plot(
        plot_t[:], 
        [ plot_x[1, :] plot_x[2, :] plot_x[3, :] ], 
        label=["r" "v_r" "v_t"],
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
        [ plot_u[1, :].-plot_u[2, :] plot_u[3, :].-plot_u[4, :] ], 
        label=["u_r" "u_t"],
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
    
    if isempty(plot_ε)
        
        plot_ref = Plots.plot(
            plots_x,
            plots_u,
            layout = (2, 1),
            # legend = false,
            # margin = 1Plots.cm,
        )
    else
        plots_ε = Plots.plot(
            plot_t, 
            [ plot_ε[:] ], 
            label=["difficulty" ],
            xlabel = "t",
            st=:scatter
        )
        
        plot_ref = Plots.plot(
            plots_x,
            plots_u,
            plots_ε,
            layout = (3, 1),
            # legend = false,
            # margin = 1Plots.cm,
        )
    end
    png(string(index, base = 10, pad = 2))
end

main()