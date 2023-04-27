using JuMP
# import Ipopt
using NLopt
using Plots

using Polynomials
using SpecialPolynomials
using SpecialFunctions

# println(NLopt.DEFAULT_OPTIONS)

# rocket = Model(NLopt.Optimizer)
# set_optimizer_attribute(rocket, "algorithm", :LD_SLSQP)
# set_silent(rocket)

# println(NLopt.DEFAULT_OPTIONS)

#-------------------------------------------------------

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
n = 128    # Time steps
γ_max = 10
γ_min = 0
T = 0.1405
m_p = 0.0749

alpha = -0.5
beta = -0.5

tau = compute_tau(n)
D = compute_D(n, alpha, beta, tau)
t_f = 3.32
t_0 = 0
using NLopt
using Polynomials
using SpecialPolynomials
using SpecialFunctions
using Plots
#-------------------------------------------------------

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
n = 64    # Time steps
γ_max = 10
γ_min = 0
T = 0.1405
m_p = 0.0749

alpha = -0.5
beta = -0.5

tau = compute_tau(n)
D = compute_D(n, alpha, beta, tau)
t_f = 3.32
t_0 = 0

#-----

rocket = Model(NLopt.Optimizer)
set_optimizer_attribute(rocket, "algorithm", :LN_COBYLA)


@variables(
    rocket,
    begin
        # t_f ≥ 0, (start = 3.32)               # Time step
        # t_0 ≥ 0               # Time step
        # # State variables
        1e2 ≥ u[1:(n + 1)] ≥ -1e2            # Velocity
        1e2 ≥ v[1:(n + 1)] ≥ -1e2            # Velocity
        1e1 ≥ r[1:(n + 1)] ≥ -1e-2, (start = 1.0)
        1e2 ≥ θ[1:(n + 1)] ≥ -1e2  
        # Control variables
        γ_min ≤ γ[1:(n + 1)] ≤ γ_max, (start = 0.001)    # Thrust
        # v[n + 1] = 1/sqrt(r[n + 1])

        # r[i = 1:(n+1)] ≥ 0, (start = i+1)
        # v[1:(n+1)] ≥ 0
    end
)

@NLobjective(rocket, Max, r[n + 1])
# @NLobjective(model, Min, v[n + 1])

# for i in 1:(n+1)
#     JuMP.@NLconstraint(model, v[i] == (10*i*r[i]-i)^2 + i)
# end

# for i in 1:n
#     JuMP.@NLconstraint(model, (v[i+1] - v[i])/(r[i+1] - r[i]) == 1.0)
# end


@NLexpressions(
    rocket,
    begin
        Dr[k = 1:(n + 1)], sum(D[k, j]*r[j] for j in 1:(n + 1))
        Dθ[k = 1:(n + 1)], sum(D[k, j]*θ[j] for j in 1:(n + 1))
        Du[k = 1:(n + 1)], sum(D[k, j]*u[j] for j in 1:(n + 1))
        Dv[k = 1:(n + 1)], sum(D[k, j]*v[j] for j in 1:(n + 1))
        r_f, 1/sqrt(r[n + 1])
        T_2, (t_f - t_0)/2
        t[k = 1:(n + 1)], (t_f - t_0)/2*tau[k] + (t_f + t_0)/2
    end
);


for k in 1:(n + 1)

    @NLconstraint(rocket, Dr[k] == T_2*u[k])
    @NLconstraint(rocket, Dθ[k] == T_2*v[k]/r[k])
    @NLconstraint(
        rocket,
        Du[k] == T_2*(
            v[k]^2.0/r[k] - r[k]^(-2.0) + T/(1 - m_p*t[k])*sin(γ[k])
        )
    )
    @NLconstraint(
        rocket,
        Dv[k] == T_2*(
            -u[k]*v[k]/r[k] + T/(1 - m_p*t[k])*cos(γ[k])
        )
    )
end

# for k in 1:n
#     @NLconstraint(rocket, γ[k + 1] - γ[k] >= 0.1)
# end

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
