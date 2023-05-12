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
        tk::Array{Real, 1}
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
                    for j in 1:nx
                        _xm[j, k] = (xn[j, _tk + i] - xn[j, _tk + i - 1])*(τm[_tk_idx][k] - tau[i]) + tau[i]
                    end
                    for j in 1:nu
                        _xm[j, k] = (un[j, _tk + i] - un[j, _tk + i - 1])*(τm[_tk_idx][k] - tau[i]) + tau[i]
                    end
                end
                for j in 1:nx
                    _xm[j, 1] = xn[j, 1]
                    _xm[j, m[_tk_idx] + 1] = xn[j, n + 1]
                end
                for j in 1:nu
                    _um[j, 1] = un[j, 1]
                    _um[j, m[_tk_idx] + 1] = un[j, n + 1]
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
                t1 = (tf - t0)*(max_idx - 2) + t0
                t2 = (tf - t0)*(max_idx    ) + t0
                (rx, ru) = convert_coefficient_n_m(n, Real[_m, _m, _m], nx, nu, xn, un, τ[n], Array{Real, 1}[τ[_m], τ[_m], τ[_m]], Real[1, max_idx - 1, max_idx + 1, n+1])
                return (3, [_m, _m, _m], Real[t0, t1, t2], Real[t1, t2, tf], rx, ru)
            else
                _m = n + 8
                (rx, ru) = convert_coefficient_n_m(n, Real[_m], nx, nu, xn, un, τ[n], Array{Real, 1}[τ[_m]], Real[1, n+1])
                return (1, [_m], Real[t0], Real[tf], rx, ru)
            end
        else
            _m = floor(n/2)
            t1 = (tf - t0)*(floor(0.4*n) - 1) + t0
            t2 = (tf - t0)*(floor(0.7*n) - 1) + t0
            (rx, ru) = convert_coefficient_n_m(n, Real[_m, _m, _m], nx, nu, xn, un, τ[n], Array{Real, 1}[τ[_m], τ[_m], τ[_m]], Real[1, floor(0.4*n), floor(0.7*n), n+1])
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
                Array{Real, 1}[[dx[i][i, k] for k in 1:(n[i] + 1)] for i in 1:ns], 
                Array{Real, 1}[[x[i][i, k] for k in 1:(n[i] + 1)] for i in 1:ns], 
                Array{Real, 1}[[fx[i][i, k] for k in 1:(n[i] + 1)] for i in 1:ns], 
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

        for idx in unique(sort(x_idx))
            deleteat!(rx, idx + offset)
            deleteat!(ru, idx + offset)
            deleteat!(rn, idx + offset)
            deleteat!(rt0, idx + offset)
            deleteat!(rtf, idx + offset)
            (_len, _m, _t0, _tf, segx, segu) = refine_segment(n[idx], nx, nu, x[idx], u[idx], fx[idx], wx, εx[:][idx], τ, t0[idx], tf[idx])
            for i in _len:1
                insert!(rx, idx + offset, segx[i])
                insert!(ru, idx + offset, segu[i])
                insert!(rn, idx + offset, _m[i])
                insert!(rt0, idx + offset, _t0[i])
                insert!(rtf, idx + offset, _tf[i])
            end
            offset = offset + _len - 1
        end

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
        na::Int,
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

        JuMP.@variables(
            rocket, 
            begin
                constraints_x[j][1] ≥ _x[i = 1:ns, j = 1:nx, k = 1:(n[i] + 1)] ≥ constraints_x[j][2] , (start = x[i][j, k])
                constraints_u[j][1] ≥ _u[i = 1:ns, j = 1:nu, k = 1:(n[i] + 1)] ≥ constraints_u[j][2] , (start = u[i][j, k])
                _a[1:na]
                _nx
                _nu
                _na 
                _na_obj
                _a_obj[1:na_obj]
            end
        )

        for i in 1:na
            JuMP.fix(_a[i], a[i]; force = true)
        end
        
        for i in 1:na_obj
            JuMP.fix(_a_obj[i], obj_a[i]; force = true)
        end

        JuMP.fix(_nx, nx; force = true)
        JuMP.fix(_nu, nu; force = true)
        JuMP.fix(_na, na; force = true)
        JuMP.fix(_na_obj, na; force = true)

        JuMP.@NLexpressions(
            rocket,
            begin
                T_2[i = 1:ns], (tf[i] - t0[i])/2
                t[i = 1:ns, k = 1:(n[i] + 1)], (tf[i] - t0[i])/2*τ[n[i]][k] + (tf[i] + t0[i])/2
            end
        )

        # JuMP.Containers.SparseAxisArray{VariableRef, 3, Tuple{Int64, Int64, Int64}}
        # println("---------------------------------------------------------------------------")
        # println(typeof(_x))

        f_sym_obj = Symbol("f_0")
        JuMP.register(rocket, f_sym_obj, nx + nu + na_obj + 4, obj_f, autodiff = true)
        obj_expr = (
            na_obj > 0 
            ? JuMP.add_nonlinear_expression(rocket, :($(f_sym_obj)($(_nx), $(_nu), $(_na_obj), $(t[ns, n[ns] + 1]), $(_x[ns, :, n[ns] + 1]...), $(_u[ns, :, n[ns] + 1]...), $(_a_obj...))))
            : JuMP.add_nonlinear_expression(rocket, :($(f_sym_obj)($(_nx), $(_nu), $(_na_obj), $(t[ns, n[ns] + 1]), $(_x[ns, :, n[ns] + 1]...), $(_u[ns, :, n[ns] + 1]...))))
        )
        @NLobjective(rocket, Max, obj_expr)
        # @NLobjective(rocket, Max, obj_f(_nx, _nu, t[end, end], _x[ns, :, n[ns] + 1]...))
        # println(JuMP.objective_function(rocket))
        
        f_sym = [Symbol("f_$(j)") for (j, ele_f) in enumerate(f)]

        for (j, ele_f) in enumerate(f)
            JuMP.register(rocket, f_sym[j], nx + nu + na + 4, ele_f, autodiff = true)
        end
            
        fx = JuMP.Containers.SparseAxisArray(
            na > 0
            ? Dict(
                ((i, j, k), JuMP.add_nonlinear_expression(rocket, :($(f_sym[j])($(_nx), $(_nu), $(_na), $(t[i, k]), $(_x[i, :, k]...), $(_u[i, :, k]...), $(_a...)))))
                for i in 1:ns for (j, ele_f) in enumerate(f) for k in 1:(n[i] + 1)
            )
            : Dict(
                ((i, j, k), JuMP.add_nonlinear_expression(rocket, :($(f_sym[j])($(_nx), $(_nu), $(_na), $(t[i, k]), $(_x[i, :, k]...), $(_u[i, :, k]...)))))
                for i in 1:ns for (j, ele_f) in enumerate(f) for k in 1:(n[i] + 1)
            )
        )
        # println("---------------------------------------------------------------------------")
        # println(typeof(fx))

        JuMP.@NLexpressions(
            rocket,
            begin
                Dx[i = 1:ns, j = 1:nx, k = 1:(n[i] + 1)], sum(D[n[i]][k, l]*_x[i, j, l] for l in 1:(n[i] + 1))
                wDx[i = 1:ns, j = 1:nx], sum(Dx[i, j, k]*W[n[i]][k] for k in 1:(n[i] + 1))
                # fx[i = 1:ns][j = 1:nx, k = 1:(n[i] + 1)], f[j](_x[i][:, k], _u[i][:, k], _a)
                # T_2[i = 1:ns], (tf[i] - t0[i])/2
                # t[i = 1:ns, k = 1:(n[i] + 1)], (tf[i] - t0[i])/2*τ[n[i]][k] + (tf[i] + t0[i])/2
            end
        )
        # println("---------------------------------------------------------------------------")
        # println(typeof(Dx))

        JuMP.@NLconstraint(rocket, con1[i = 1:ns, j = 1:nx, k = 1:(n[i] + 1)], Dx[i, j, k] == T_2[i]*fx[i, j, k])

        if (ns >= 1)
            JuMP.@NLconstraint(rocket,con2[[1], j = 1:nx],  _x[ns, j, n[ns] + 1] - _x[ns, j, 1] == wDx[ns, j])
            if (ns > 1)
                JuMP.@NLconstraint(rocket, con2[i = 2:ns, j = 1:nx], _x[i, j, 1] - _x[i - 1, j, 1] == wDx[i, j])
            end
        end

        for (j, ele_f) in enumerate(constraints_eq)
            f_sym_eq  = Symbol("f_eq_$(j)")
            tmp = size(ele_f[2], 1)
            JuMP.register(rocket, f_sym_eq, tmp + 2*nx + 2*nu + 3, ele_f[1], autodiff = true)
            if tmp > 0
                JuMP.add_nonlinear_constraint(rocket, :($(f_sym_eq)($(_nx), $(_nu), $(tmp), $(_x[1, :, 1]...), $(_x[ns, :, n[ns] + 1]...), $(_u[1, :, 1]...), $(_u[ns, :, n[ns] + 1]...), $(ele_f[2]...)) == 0))
            else
                a = JuMP.add_nonlinear_constraint(rocket, :($(f_sym_eq)($(_nx), $(_nu), $(tmp), $(_x[1, 1:1:nx, 1]...), $(_x[ns, 1:1:nx, n[ns] + 1]...), $(_u[1, :, 1]...), $(_u[ns, :, n[ns] + 1]...)) == 0))
                println(:($(_x[1,:,1] ...)))
                println(:($([_x[1,j,1] for j in 1:nx] ...)))
                println(a)
            end
            # @NLconstraint(rocket, conx[j], ele_f[1](_x, _u, ele_f[2]) == 0)
        end

        for (j, ele_f) in enumerate(constraints_le)
            f_sym_le  = Symbol("f_le_$(j)")
            tmp = size(ele_f[2], 1)
            JuMP.register(rocket, f_sym_le, tmp + 2*nx + 2*nu + 3, ele_f[1], autodiff = true)
            if tmp > 0
                JuMP.add_nonlinear_constraint(rocket, :($(f_sym_le)($(_nx), $(_nu), $(tmp), $(_x[0, :, 0]...), $(_x[ns, :, n[ns] + 1]...), $(_u[0, :, 0]...), $(_u[ns, :, n[ns] + 1]...), $(ele_f[2]...)) <= 0))
            else
                JuMP.add_nonlinear_constraint(rocket, :($(f_sym_le)($(_nx), $(_nu), $(tmp), $(_x[0, :, 0]...), $(_x[ns, :, n[ns] + 1]...), $(_u[0, :, 0]...), $(_u[ns, :, n[ns] + 1]...)) <= 0))
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
        # print(solution_summary(rocket))

        # println("got ", objective_value(rocket))
        # println(JuMP.value(r_f))

        rx = JuMP.value.(_x)
        ru = JuMP.value.(_u)
        rDx = JuMP.value.(Dx)
        rfx = JuMP.value.(fx)

        _rx = [[rx[i, j, k] for j in 1:nx, k in 1:(n[i] + 1)] for i in 1:ns]
        _ru = [[ru[i, j, k] for j in 1:nu, k in 1:(n[i] + 1)] for i in 1:ns]
        _rDx = [[rDx[i, j, k] for j in 1:nx, k in 1:(n[i] + 1)] for i in 1:ns]
        _rfx = [[rfx[i, j, k] for j in 1:nx, k in 1:(n[i] + 1)] for i in 1:ns]

        delete_x = all_variables(rocket)
        delete(rocket, delete_x)
        unregister(model, :delete_x)

        delete_con = all_nonlinear_constraints(rocket)
        delete(rocket, delete_con)
        unregister(model, :delete_con)
        
        return (
            _rx,
            _ru,
            _rDx,
            _rfx,
        )
end


# function add_plot(y, ylabel, t)
#     return Plots.plot!(
#         value.(t)[:],
#         value.(y)[:];
#         xlabel = "t",
#         ylabel = ylabel,
#         st=:scatter
#     )
# end

# function my_plot(y, ylabel, t)
#     return Plots.plot(
#         value.(t)[:],
#         value.(y)[:];
#         xlabel = "t",
#         ylabel = ylabel,
#         st=:scatter
#     )
# end

# Plots.plot(
#     Plots.plot(
#         value.(t)[:], 
#         [ value.(r)[:] value.(θ)[:] value.(u)[:] value.(v)[:] ], 
#         label=["r" "θ" "u" "v"],
#         xlabel = "t",
#         st=:scatter
#     ),
#     my_plot(γ, "γ", t),
#     layout = (2, 1),
#     # legend = false,
#     # margin = 1Plots.cm,
# )

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
    n = Int[7 for i in 1:ns]
    x = [Real[0.0 for j in 1:nx, k in 1:(n[i] + 1)] for i in 1:ns]
    cx = Tuple{Real, Real}[(1.0, -1.0) for j in 1:nx]
    u = [Real[0.0 for j in 1:nu, k in 1:(n[i] + 1)] for i in 1:ns]
    cu = Tuple{Real, Real}[(1.0, -1.0) for j in 1:nu]
    bx = Tuple{Real, Real}[(NaN, NaN) for _ in 1:nx]
    bu = Tuple{Real, Real}[(NaN, NaN) for _ in 1:nu]
    c_eq = Tuple{Function, Array{Real, 1}}[
        (
            ((nx, nu, na, r0, θ0, v_r0, v_θ0, rf, θf, v_rf, v_θf, γ0, γf) -> v_rf - 1/sqrt(rf)), []
        )
    ]
    c_le = Tuple{Function, Array{Real, 1}}[]
    f = Function[
        ((nx, nu, na, t, r, θ, v_r, v_θ, γ, T, m_p) -> v_r),
        ((nx, nu, na, t, r, θ, v_r, v_θ, γ, T, m_p) -> v_θ/r),
        ((nx, nu, na, t, r, θ, v_r, v_θ, γ, T, m_p) -> v_θ^2.0/r - r^(-2.0) + T/(1 - m_p*t)*sin(γ)),
        ((nx, nu, na, t, r, θ, v_r, v_θ, γ, T, m_p) -> -v_r*v_θ/r + T/(1 - m_p*t)*cos(γ))
    ]
    a = Real[0.1405, 0.0749]
    f_obj = (nx, nu, na, t, x...) -> x[1]
    na_obj = 0
    a_obj = Real[]
    t0 = Real[0.0]
    tf = Real[3.32]

    ε_tol = 2


    while (true)
        println(x)
        (x, u, dx, fx) = solve_NLP(ns, nx, nu, na, n, x, cx, u, cu, bx, bu, c_eq, c_le, f, a, f_obj, na_obj, a_obj, τ, t0, tf, D, W)
        
        (ns, n, x, u, t0, tf) = refine(ns, nx, nu, n, x, u, dx, fx, τ, t0, tf, ε_tol, 0.8, W)
        if ns < 0
            println(x)
            break
        end
    end
end

main()