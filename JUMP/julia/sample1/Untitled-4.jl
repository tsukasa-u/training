using NLopt
using Plots

using Polynomials
using SpecialPolynomials
using SpecialFunctions


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
n = 16    # Time steps
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

function myfunc(x::Vector, grad::Vector)
    m::Int = length(grad)/5
    if length(grad) > 0
        grad .= 0.0
        grad[m] = -1.0
    end
    println(x)
    return -x[m]
end

function myconstraint_r(x::Vector, grad::Vector, k, m, D)
    if length(grad) > 0
        grad[1:m] .= D[k, :]
        grad[(m + 1):5*m] .= 0.0
        grad[2*m + k] = -1.0
    end
    sum(D[k, j]*x[j] for j in 1:m) - x[2*m + k]
end

function myconstraint_θ(x::Vector, grad::Vector, k, m, D)
    if length(grad) > 0
        grad[1:m] .= 0.0
        grad[(  m + 1):2*m] .= D[k, :]
        grad[(2*m + 1):5*m] .= 0.0
        grad[      k] -= -x[3*m + k]/x[k]^2.0
        grad[3*m + k] -= 1.0/x[k]
    end
    sum(D[k, j]*x[m + j] for j in 1:m) - x[3*m + k]/x[k]
end

function myconstraint_u(x::Vector, grad::Vector, k, m, D)
    if length(grad) > 0
        grad[1:2*m] .= 0.0
        grad[(2*m + 1):3*m] .= D[k, :]
        grad[(3*m + 1):5*m] .= 0.0
        grad[      k] -= -x[3*m + k]^2.0/x[k]^2.0 + 2.0/x[k]^3.0
        grad[3*m + k] -= 2.0*x[3*m + k]/x[k]
        grad[4*m + k] -= T*cos(γ[k])/(1.0 - m_p*t)
    end
    sum(D[k, j]*x[2*m + j] for j in 1:m) - (x[3*m + k]^2.0/x[k] - 1.0/x[k]^2.0 + T*sin(γ[k])/(1.0 - m_p*t))
end

function myconstraint_u(x::Vector, grad::Vector, k, m, D)
    if length(grad) > 0
        grad[1:3*m] .= 0.0
        grad[(3*m + 1):4*m] .= D[k, :]
        grad[(4*m + 1):5*m] .= 0.0
        grad[      k] -= -x[2*m + k]*x[3*m + k]/x[k]^2.0
        grad[2*m + k] -= -x[3*m + k]/x[k]
        grad[3*m + k] -= -x[2*m + k]/x[k]
        grad[4*m + k] -= -T*sin(γ[k])/(1.0 - m_p*t)
    end
    sum(D[k, j]*x[3*m + j] for j in 1:m) - (x[2*m + k]*x[3*m + k]/x[k] + T*cos(γ[k])/(1.0 - m_p*t))
end

function myconstraint_b(x::Vector, grad::Vector, k::Int, b::Real)
    if length(grad) > 0
        grad .= 0.0
        grad[k] = 1.0
    end
    x[k] - b
end

function myconstraint_f(x::Vector, grad::Vector, k::Int, m::Int)
    if length(grad) > 0
        grad .= 0.0
        grad[3*m + m] = 1.0
        grad[      m] -= -0.5*x[k]^(-1.5)
    end
    return x[3*m + m] - x[m]^(-0.5)
end


opt = Opt(:LD_SLSQP, 5*(n + 1))
# opt.lower_bounds .= -Inf
opt.lower_bounds .= 1e-3
opt.xtol_rel = 1e-4

opt.min_objective = myfunc
for k in 1:(n + 1)
    equality_constraint!(opt, (x,g) -> myconstraint_r(x, g, k, n+1, D), 1e-8)
    # equality_constraint!(opt, (x,g) -> myconstraint_θ(x, g, k, n+1, D), 1e-8)
    # equality_constraint!(opt, (x,g) -> myconstraint_u(x, g, k, n+1, D), 1e-8)
    # equality_constraint!(opt, (x,g) -> myconstraint_v(x, g, k, n+1, D), 1e-8)
end

equality_constraint!(opt, (x,g) -> myconstraint_b(x, g,             1, 1.0), 1e-8)
equality_constraint!(opt, (x,g) -> myconstraint_b(x, g,    n + 1  + 1, 0.0), 1e-8)
equality_constraint!(opt, (x,g) -> myconstraint_b(x, g, 2*(n + 1) + 1, 0.0), 1e-8)
equality_constraint!(opt, (x,g) -> myconstraint_b(x, g, 3*(n + 1) + 1, 1.0), 1e-8)

equality_constraint!(opt, (x,g) -> myconstraint_b(x, g, 2*(n + 1) + 1, 0.0), 1e-8)
equality_constraint!(opt, (x,g) -> myconstraint_f(x, g, n + 1, n + 1), 1e-8)


# inequality_constraint!(opt, (x,g) -> myconstraint(x,g,-1,1), 1e-8)
x = fill(0.0001, 5*(n + 1))
(minf,minx,ret) = optimize(opt, x)
numevals = opt.numevals # the number of function evaluations
# println("got $minf at $minx after $numevals iterations (returned $ret)")
println("returned $ret")