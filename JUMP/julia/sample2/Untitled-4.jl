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

function compute_difficulty(n::Int, f::Array{Real, 2}, w::Tuple{Real}, ε::Tuple{Real}, τ::Array{Real, 1}, t0::Real, tf::Real)::Array{Real, 1}
    F_k = [sum(abs(ε[j]*f[j, k]/w[j]) for (j, ele) in enumerate(ε)) for k in 1:(n+1)]
    F_1 = [(F_k[k+1] - F_k[k])*2/(tf-t0)/(τ[k+1] - τ[k]) for k in 1:n]
    insert!(F_1, 1, -1)
    F_2 = [(F_k[k] - F_k[k-1])*2/(tf-t0)/(τ[k+1] - τ[k-1]) for k in 2:n]
    insert!(F_2, 1, -1)

    sum_F = sum(F_2[l] for l in 2:n)
    G_k = [(n - 1)*abs(F_2[k])/sum_F for k in 2:n]
    return insert(G_k, 1, -1)
end

function convert_coefficient_n_m(
        n::Real,
        m::Tuple{Real},
        nx::Real,
        nu::Real,
        xn::Array{Real, 2}, 
        un::Array{Real, 2}, 
        τn::Array{Real}, 
        τm::Tuple{Array{Real, 1}}, 
        tk::Tuple{Real}
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
        w::Tuple{Real}, 
        ε::Tuple{real}, 
        τ::Tuple{Array{Real, 1}}, 
        t0::Real, 
        tf::Real, 
        g_tolerance::Real=2
    )::(Real, Tuple{Real}, Tuple{Real}, Tuple{Real}, Array{Array{Real, 2}, 1})
        if n <= 40
            G_k = compute_difficulty(n, f, w, ε, τ, t0, tf)
            (max_idx, max_val) = findMax(G_k)
            if max_val > g_tolerance  #G_k tolerance 2
                _m = max(floor(n/3), 7)
                t1 = (tf - t0)*(max_idx - 2) + t0
                t2 = (tf - t0)*(max_idx    ) + t0
                (rx, ru) = convert_coefficient_n_m(n, (_m, _m, _m), nx, nu, xn, un, τ[n], (τ[_m], τ[_m], τ[_m]), (1, max_idx - 1, max_idx + 1, n+1))
                return (3, (_m, _m, _m), (t0, t1, t2), (t1, t2, tf), rx, ru)
            else
                _m = n + 8
                (rx, ru) = convert_coefficient_n_m(n, (_m), nx, nu, xn, un, τ[n], (τ[_m]), (1, n+1))
                return (1, (_m), (t0), (tf), rx, ru)
            end
        else
            _m = floor(n/2)
            t1 = (tf - t0)*(floor(0.4*n) - 1) + t0
            t2 = (tf - t0)*(floor(0.7*n) - 1) + t0
            (rx, ru) = convert_coefficient_n_m(n, (_m, _m, _m), nx, nu, xn, un, τ[n], (τ[_m], τ[_m], τ[_m]), (1, floor(0.4*n), floor(0.7*n), n+1))
            return (3, (_m, _m, _m), (t0, t1, t2), (t1, t2, tf), rx, ru)
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
        τ::Tuple{Array{Real, 1}},
        t0::Array{Real, 1},
        tf::Array{Real, 1},
        ε_tol::Real, 
        α::Real,
        W::Tuple{Array{Real, 1}}
    )::(Real, Array{Int, 1}, Array{Array{Real, 2}, 1}, Array{Array{Real, 2}, 1}, Array{Real, 1}, Array{Real, 1})
        εx::Array{Array{Real, 1}, 1} = Array{Real, 1}[]
        wx::Array{Real, 1} = Real[]
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
        x::Tuple{Array{Real, 2}}, 
        constraints_x::Tuple{Array{Real, 2}},
        u::Tuple{Array{Real, 2}}, 
        constraints_u::Tuple{Array{Real, 2}},
        f::Tuple{Function}
        τ::Tuple{Array{Real, 1}},
        t0::Array{Real, 1},
        tf::Array{Real, 1},
        D::Tuple{Array{Real, 2}},
        W::Tuple{Array{Real, 1}}
    )::(Tuple{Array{Real, 2}}, Tuple{Array{Real, 2}}, Tuple{Array{Real, 2}})
        # n = 64    # Time steps
        T = 0.1405
        m_p = 0.0749

        # alpha = -0.5
        # beta = -0.5

        # tau = compute_tau(n)
        # D = compute_D(n, alpha, beta, tau)
        # t_f = 3.32
        # t_0 = 0

        #-----

        rocket = Model(NLopt.Optimizer)
        set_optimizer_attribute(rocket, "algorithm", :LD_SLSQP)

        @variables( rocket, constraints_x[1][i][j, k] ≥ _x[i = 1:ns][j = 1:nx, k = 1:(n[i] + 1)] ≥ constraints_x[0][i][j, k] )
        @variables( rocket, constraints_u[1][i][j, k] ≥ _u[i = 1:ns][j = 1:nu, k = 1:(n[i] + 1)] ≥ constraints_u[0][i][j, k] )

        @NLobjective(rocket, Max, _x[ns][0, n + 1])

        @NLexpressions(
            rocket,
            begin
                Dx[i = 1:ns][j = 1:nx, k = 1:(n[i] + 1)], sum(D[n[i]][k, l]*_x[i][j, l] for l in 1:(n[i] + 1))
                wDx[i = 1:ns][j = 1:nx], sum(Dx[i][j, k]*W[k] for k in 1:(n[i] + 1))
                fx[i = 1:ns][j = 1:nx, k = 1:(n + 1)], f[j](_x[i][:, k], _u[i][:, k], _a)
                T_2[i = 1:ns], (tf[i] - t0[i])/2
                t[i = 1:ns][k = 1:(n[i] + 1)], (tf[i] - t0[i])/2*τ[n[i]][k] + (tf[i] + t0[i])/2
            end
        );

        @NLconstraint(rocket, Dx[i = 1:ns][j = 1:nx, k = 1:(n[i] + 1)] == T_2[i]*fx[i][j, k])

        if (ns >= 1)
            @NLconstraint(rocket, _x[ns][j = 1:nx, n[i] + 1] - _x[ns][j = 1:nx, 1] == wDx[ns][j])
            if (ns > 1)
                @NLconstraint(rocket, _x[i = 1:ns][j = 1:nx, 1] - _x[i = 1:(ns - 1)][j = 1:nx, 1] == wDx[i][j])
            end
        end

        @NLconstraint(rocket, v[n + 1] == 1/sqrt(r[n + 1]))

        fix(r[1], 1.0; force = true)
        fix(θ[1], 0.0; force = true)
        fix(u[1], 0.0; force = true)
        fix(v[1], 1.0; force = true)
        fix(γ[1], 0.001; force = true)
        fix(u[n + 1], 0.0; force = true)

        println("Solving...")
        JuMP.optimize!(rocket)
        print(solution_summary(rocket))

        println("got ", objective_value(rocket))
        println(JuMP.value(r_f))
end


function add_plot(y, ylabel, t)
    return Plots.plot!(
        value.(t)[:],
        value.(y)[:];
        xlabel = "t",
        ylabel = ylabel,
        st=:scatter
    )
end

function my_plot(y, ylabel, t)
    return Plots.plot(
        value.(t)[:],
        value.(y)[:];
        xlabel = "t",
        ylabel = ylabel,
        st=:scatter
    )
end

Plots.plot(
    Plots.plot(
        value.(t)[:], 
        [ value.(r)[:] value.(θ)[:] value.(u)[:] value.(v)[:] ], 
        label=["r" "θ" "u" "v"],
        xlabel = "t",
        st=:scatter
    ),
    my_plot(γ, "γ", t),
    layout = (2, 1),
    # legend = false,
    # margin = 1Plots.cm,
)

function main() {

    W = tuple(compute_W(n) in 1:50)
    τ = tuple(compute_tau(n) in 1:50)
    D = tuple(compute_D(n, -0.5, -0.5, τ[n]) in 1:50)

    while (true)
        (x, dx, fx) = solve_NLP(ns, nx, nu, n, x, bx, u, bu, τ, t0, tf, D, W)
        
        (ns, n, x, t0, tf) = refine(ns, nx, nu, n, x, u, dx, fx, τ, t0, tf, ε_tol, 0.8, W)
        if ns<0
            break
        end
    end
}

main()