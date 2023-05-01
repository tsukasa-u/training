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
        w::Array{Real ,1}
    )::Real
        return sum(abs(dx[k] - f[k])*w[k] for k in 1:(n+1))
end

function compute_w(
        ns::Int, 
        n::Array{Int, 1}, 
        dx::Array{Real, 2}, 
        x::Array{Real, 2}
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
    )::(Array{Real, 1}, Real)
        w::Real = compute_w(ns, n, dx, x)
        return ([compute_η(n[i], dx[i], f[i], W[n[i]])/w for i in 1:ns], w)
end

function check_ε_tolerance(ε_tol::Real, ns::Real, ε::Array{Array{Real, 1}, 1})::Real
    tmp = maximum([sum(ele[i] for i in ns) for ele in ε])
    return tmp < ε_tol ? tmp : -1
end

function select_refine_segment(ε_tol::Real, α::Real, ns::Int, ε::Array{Array{Real, 1}, 1})::Array{Int, 1}
    ret_index::Array{Int, 1} = fill(-1, ns)
    for (j, ele) in enumerate(ε)
        idx = sortperm(ele)
        sum = 0
        r = 0
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
    return insert(G_k, 1, NaN)
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
    )::(Array{Array{Real, 2}, 1}, Array{Array{Real, 2}, 1})
        
        # xm = tuple(Tuple(Real[]))
        xm = Array{Real, 2}[]
        um = Array{Real, 2}[]
        # for tk_idx in 2:length(tk)
        tk_size = size(tk, 1)
        for tk_idx in 2:tk_size
            tau = τn[tk[tk_idx-1]:tk[tk_idx]]

            if size(tau, 1) > 0

                tau = tau - tau[0]
                tau = tau*tau[end]

                _xm = fill(0.0, (nx, m + 1))
                _um = fill(0.0, (nx, m + 1))

                i = 1
                for k in 2:m
                    if τm[k] > tau[i] && τm[k] <= tau[i + 1]
                        for (j, ele_n) in enumerate(xn)
                            _xm[j][k] = (ele_n[i + 1] - ele_n[i])*(τm[k] - tau[i]) + tau[i]
                        end
                        for (j, ele_n) in enumerate(un)
                            _um[j][k] = (ele_n[i + 1] - ele_n[i])*(τm[k] - tau[i]) + tau[i]
                        end
                    else
                        i += 1
                        continue
                    end
                end
                for (j, ele_n) in enumerate(xn)
                    _xm[j][1] = ele_n[1]
                    _xm[j][m + 1] = ele_n[n + 1]
                end
                for (j, ele_n) in enumerate(un)
                    _um[j][1] = ele_n[1]
                    _um[j][m + 1] = ele_n[n + 1]
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
        ε::Array{real}, 
        τ::Array{Array{Real, 1}, 1}, 
        t0::Real, 
        tf::Real, 
        g_tolerance::Real=2
    )::(Real, Array{Real, 1}, Array{Real, 1}, Array{Real, 1}, Array{Array{Real, 2}, 1})
        if n <= 40
            G_k = compute_difficulty(n, f, w, ε, τ, t0, tf)
            (max_idx, max_val) = findMax(G_k)
            if max_val > g_tolerance  #G_k tolerance 2
                _m = max(floor(n/3), 7)
                t1 = (tf - t0)*(max_idx - 2) + t0
                t2 = (tf - t0)*(max_idx    ) + t0
                (rx, ru) = convert_coefficient_n_m(n, [_m, _m, _m], nx, nu, xn, un, τ[n], [τ[_m], τ[_m], τ[_m]], [1, max_idx - 1, max_idx + 1, n+1])
                return (3, [_m, _m, _m], [t0, t1, t2], [t1, t2, tf], rx, ru)
            else
                _m = n + 8
                (rx, ru) = convert_coefficient_n_m(n, [_m], nx, nu, xn, un, τ[n], [τ[_m]], [1, n+1])
                return (1, [_m], [t0], [tf], rx, ru)
            end
        else
            _m = floor(n/2)
            t1 = (tf - t0)*(floor(0.4*n) - 1) + t0
            t2 = (tf - t0)*(floor(0.7*n) - 1) + t0
            (rx, ru) = convert_coefficient_n_m(n, [_m, _m, _m], nx, nu, xn, un, τ[n], [τ[_m], τ[_m], τ[_m]], [1, floor(0.4*n), floor(0.7*n), n+1])
            return (3, [_m, _m, _m], [t0, t1, t2], [t1, t2, tf], rx, ru)
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
    )::(Real, Array{Int, 1}, Array{Array{Real, 2}, 1}, Array{Array{Real, 2}, 1}, Array{Real, 1}, Array{Real, 1})
        # εx::Array{Array{Real, 1}, 1} = Array{Real, 1}[]
        # wx::Array{Real, 1} = Real[]
        εx = Array{Real, 1}[]
        wx = Real[]
        for i in 1:nx
            tmp = compute_ε(ns, n, dx[:][i, :], x[:][i, :], fx[:][i, :], W)
            push!(εx, tmp[0])
            push!(wx, tmp[1])
        end

        if check_ε_tolerance(ε_tol, ns, εx) > 0
            return (-ns, n, x, t0, tf)
        end

        x_idx = select_refine_segment(ε, α, ns, εx)

        rx = copy(x)
        ru = copy(u)
        rn = copy(n)
        rt0 = copy(t0)
        rtf = copy(tf)
        offset = 0

        for idx in unique(sort(x_idx))
            delete!(rx, idx + offset)
            delete!(ru, idx + offset)
            delete!(rn, idx + offset)
            delete!(rt0, idx + offset)
            delete!(rtf, idx + offset)
            (_len, _m, _t0, _tf, segx, segu) = refine_segment(n[idx], nx, nu, x[idx], u[idx], fx[idx], wx, εx, τ, t0[idx], tf[idx])
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

function compute_W(n)
    fill(pi/n, n)
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
        constraints_x::Array{Array{Tuple{Real, Real}, 2}, 1},
        u::Array{Array{Real, 2}, 1}, 
        constraints_u::Array{Array{Tuple{Real, Real}, 2}, 1},
        bound_x::Array{Tuple{Real, Real}, 1},
        bound_u::Array{Tuple{Real, Real}, 1},
        constraints_eq::Array{Tuple{Function, Array{Real, 1}}},
        constraints_le::Array{Tuple{Function, Array{Real, 1}}},
        f::Array{Function},
        τ::Array{Array{Real, 1}},
        t0::Array{Real, 1},
        tf::Array{Real, 1},
        D::Array{Array{Real, 2}},
        W::Array{Array{Real, 1}}
    )::(Array{{Array{Real, 2}}, 1}, Array{Array{Real, 2}, 1}, Array{Array{Real, 2}, 1}, Array{Array{Real, 2}, 1})
        # n = 64    # Time steps
        # T = 0.1405
        # m_p = 0.0749

        # alpha = -0.5
        # beta = -0.5

        # tau = compute_tau(n)
        # D = compute_D(n, alpha, beta, tau)
        # t_f = 3.32
        # t_0 = 0

        #-----

        rocket = Model(NLopt.Optimizer)
        set_optimizer_attribute(rocket, "algorithm", :LD_SLSQP)

        @variables(
            rocket, 
            begin
                constraints_x[i][j, k][1] ≥ _x[i = 1:ns][j = 1:nx, k = 1:(n[i] + 1)] ≥ constraints_x[i][j, k][0] , (start = x[i][j, k])
                constraints_u[i][j, k][1] ≥ _u[i = 1:ns][j = 1:nu, k = 1:(n[i] + 1)] ≥ constraints_u[i][j, k][0] , (start = u[i][j, k])
            end
        )

        @NLobjective(rocket, Max, _x[ns][0, n + 1])

        # fx = Array{NonlinearExpression, 2}[]

        # for (j, ele_f) in enumerate(f)
        #     f_sym  = Symbol("f_$(j)")
        #     register(rocket, f_sym, 3, ele_f, autodiff = true)
        #     # register(rocket, :f_$i, 3, ele_f, autodiff = true)

        #     _fx = [[add_nonlinear_expression(rocket, :($(f_sym)($(_x[i][:, k]), $(_u[i][:, k]), $(_a)))) for k = 1:(n + 1)] for i = 1:ns]
        #     # for i = 1:ns
        #     #         for k = 1:(n + 1)
        #     #             add_nonlinear_expression(rocket, :($(f_sym)($(_x[i][:, k]), $(_u[i][:, k]), $(_a))))
        #     #         end
        #     # end
        # end
        # fx = (Array{NonlinearExpression, 2})[
        fx = [
            [
                (begin
                    f_sym  = Symbol("f_$(j)")
                    register(rocket, f_sym, 3, ele_f, autodiff = true)
                    return add_nonlinear_expression(rocket, :($(f_sym)($(_x[i][:, k]), $(_u[i][:, k]), $(_a))))
                end)
                for (j, ele_f) in enumerate(f), k = 1:(n + 1)
            ] for i = 1:ns
        ]

        @NLexpressions(
            rocket,
            begin
                Dx[i = 1:ns][j = 1:nx, k = 1:(n[i] + 1)], sum(D[n[i]][k, l]*_x[i][j, l] for l in 1:(n[i] + 1))
                wDx[i = 1:ns][j = 1:nx], sum(Dx[i][j, k]*W[k] for k in 1:(n[i] + 1))
                # fx[i = 1:ns][j = 1:nx, k = 1:(n[i] + 1)], f[j](_x[i][:, k], _u[i][:, k], _a)
                T_2[i = 1:ns], (tf[i] - t0[i])/2
                t[i = 1:ns][k = 1:(n[i] + 1)], (tf[i] - t0[i])/2*τ[n[i]][k] + (tf[i] + t0[i])/2
            end
        );

        @NLconstraint(rocket, con1[i = 1:ns][j = 1:nx, k = 1:(n[i] + 1)], Dx[i][j, k] == T_2[i]*fx[i][j, k])

        if (ns >= 1)
            @NLconstraint(rocket,con2[1][j = 1:nx],  _x[ns][j, n[i] + 1] - _x[ns][j, 1] == wDx[ns][j])
            if (ns > 1)
                @NLconstraint(rocket, con2[i = 2:ns][j = 1:nx], _x[i][j, 1] - _x[i - 1][j, 1] == wDx[i][j])
            end
        end

        for (j, ele_f) in enumerate(constraints_eq)
            f_sym  = Symbol("f_$(j)")
            register(rocket, f_sym, 3, ele_f, autodiff = true)
            add_nonlinear_constraint(rocket, :($(f_sym)($(_x[i][:, k]), $(_u[i][:, k]), $(ele_f[2]))))
            # @NLconstraint(rocket, conx[j], ele_f[1](_x, _u, ele_f[2]) == 0)
        end

        for (j, ele_f) in enumerate(constraints_le)
            f_sym  = Symbol("f_$(j)")
            register(rocket, f_sym, 3, ele_f, autodiff = true)
            add_nonlinear_constraint(rocket, :($(f_sym)($(_x[i][:, k]), $(_u[i][:, k]), $(ele_f[2]))))
            # @NLconstraint(rocket, conu[j], ele_f[1](_x, _u, ele_f[2]) == 0)
        end

        for (i, ele) in enumerate(bound_x)
            if !isnan(ele[1]) 
                fix(_x[1][i, 1], ele[1]; force = true)
            end
            if !isnan(ele[2]) 
                fix(_x[ns][i, n[ns] + 1], ele[2]; force = true)
            end
        end

        for (i, ele) in enumerate(bound_u)
            if !isnan(ele[1]) 
                fix(_u[1][i, 1], ele[1]; force = true)
            end
            if !isnan(ele[2]) 
                fix(_u[ns][i, n[ns] + 1], ele[2]; force = true)
            end
        end

        # println("Solving...")
        JuMP.optimize!(rocket)
        # print(solution_summary(rocket))

        # println("got ", objective_value(rocket))
        # println(JuMP.value(r_f))

        rx = value.(_x)
        ru = value.(_u)
        rDx = value.(Dx)
        rfx = value.(fx)

        # delete_x = all_variables(rocket)
        # delete(rocket, delete_x)
        # unregister(model, :delete_x)

        # delete_con = all_nonlinear_constraints(rocket)
        # delete(rocket, delete_con)
        # unregister(model, :delete_con)
        
        return (rx, ru, rDx, rfx)
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
    W = Real[compute_W(n) n in n_min:n_max]
    τ = Real[compute_tau(n) n in n_min:n_max]
    D = Real[compute_D(n, -0.5, -0.5, τ[n]) n in n_min:n_max]

    while (true)
        (x, u, dx, fx) = solve_NLP(ns, nx, nu, n, x, cx, u, cu, bx, bu, c_eq, c_le, f, τ, t0, tf, D, W)
        
        (ns, n, x, t0, tf) = refine(ns, nx, nu, n, x, u, dx, fx, τ, t0, tf, ε_tol, 0.8, W)
        if ns<0
            break
        end
    end
end

main()