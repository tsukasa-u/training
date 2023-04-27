using JuMP
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
n = 32    # Time steps
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

model = Model(NLopt.Optimizer)
set_optimizer_attribute(model, "algorithm", :LD_SLSQP)


@variables(
    model,
    begin
        # t_f ≥ 0, (start = 3.32)               # Time step
        # t_0 ≥ 0               # Time step
        # # State variables
        u[1:(n + 1)] ≥ -1e4            # Velocity
        v[1:(n + 1)] ≥ -1e4            # Velocity
        r[1:(n + 1)] ≥ -1e4
        θ[1:(n + 1)] ≥ -1e4  
        # Control variables
        γ_min ≤ γ[1:(n + 1)] ≤ γ_max, (start = 0.001)    # Thrust
        # v[n + 1] = 1/sqrt(r[n + 1])

        x[i = 1:(n+1)] ≥ 0, (start = i+1)
        y[1:(n+1)] ≥ 0
    end
)

@NLobjective(model, Min, y[n + 1])

for i in 1:(n+1)
    JuMP.@NLconstraint(model, y[i] == (i*x[i]-i)^2 + i)
end

for i in 1:n
    JuMP.@NLconstraint(model, (y[i+1] - y[i])/(x[i+1] - x[i]) == 1.0)
end

# for i in 2:3
#     JuMP.@NLconstraint(model, ((y[i+1] - y[i])/(x[i+1] - x[i]) + (y[i] - y[i-1])/(x[i] - x[i-1]))/(x[i+1] - 2*x[i] + x[i-1])*2 == 0.0)
# end

# JuMP.set_start_value(x, [1.0 for i in 1::(n+1)])
# JuMP.set_start_value(y, [1.0 for i in 1::(n+1)])

println("Solving...")
JuMP.optimize!(model)
print(solution_summary(model))

println("got ", objective_value(model), " at ", [value(x[1]), value(y[1])])

function my_plot(y, ylabel, t)
    return Plots.plot(
        value.(t)[:],
        value.(y)[:];
        xlabel = "x",
        ylabel = ylabel,
        st=:scatter
    )
end

Plots.plot(
    my_plot(y, "y", x),
    layout = (1, 1),
    legend = false,
    margin = 1Plots.cm,
)
