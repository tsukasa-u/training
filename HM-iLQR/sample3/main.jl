using RobotZoo
using RobotDynamics
using LinearAlgebra
using Plots

include("HM-iLQR.jl")
include("M-iLQR.jl")
include("func.jl")
using .HM_iLQR
using .M_iLQR
using .func


function main()
    
    β = 1.0

    # model = RobotZoo.Cartpole()
    # n,m = RobotDynamics.dims(model)

    funcs = func.FUNC()
    funcs._h = 0.05
    # funcs.f = (x, u) -> func.dynamics_rk4(x, u, (x, u) -> RobotZoo.dynamics(model, x, u), funcs._h)
    
    funcs.f = (x, u) -> func.dynamics_rk4(
        x, 
        u, 
        (x, u) -> begin
            r, v_r, v_θ = x
            u_r_p, u_r_m, u_θ_p, u_θ_m = u
            return [
                v_r,
                v_θ^2.0/r - 1.0/r^2.0 + (u_r_p - u_r_m)/β,
                -v_r*v_θ/r + (u_θ_p - u_θ_m)/β
            ]
        end, 
        funcs._h
    )

    n = 3
    m = 4
    Q = Diagonal(0.05*ones(n))
    R = Diagonal(ones(m))
    R_ = 0.055*ones(m)
    Qn = Array(100.0*I(n))
    xgoal = [4.0, 0.0, 0.5]
    # funcs.l = (x, u) -> 0.5*((x-xgoal)'*Q*(x-xgoal) + u'*R*u)
    funcs.lf = (x) -> 0.50*((x-xgoal)'*Qn*(x-xgoal))
    funcs.lf__ = (x) -> 0.0 + 0.50*((x-xgoal)'*Qn*(x-xgoal))/100.0
    # funcs.l = (x, u) -> 0.5*((x-xgoal)'*Q*(x-xgoal) + R_'*u/β)
    funcs.l = (x, u) -> R_'*u/β
    funcs.ge = (x, u) -> begin
        r, v_r, v_θ = x
        u_r_p, u_r_m, u_θ_p, u_θ_m = u
        return [
            u_r_p - u_r_m > 0.0 ? u_r_m : u_r_p,
            u_θ_p - u_θ_m > 0.0 ? u_θ_m : u_θ_p,
        ]
    end
    funcs.gef = (x) -> [0.0]
    funcs.gl = (x, u) -> begin
        r, v_r, v_θ = x
        u_r_p, u_r_m, u_θ_p, u_θ_m = u
        return [u_r_p - β*0.01, u_θ_p - β*0.01, u_r_m - β*0.01, u_θ_m - β*0.01]*100.0
    end
    funcs.glf = (x) -> [0.0]

    N = 100
    M = 1
    MaxIter = 1000
    ϵ_v = 1e-6
    d_max = 1e-6
    
    _, _, L = M_iLQR.init(N, M)
    X_init = func.Marray(L, M, N+1, (n,), [ele for _ in 1:N+1, ele in [1.0, 0.0, 1.0]])
    X_init.limit_lower = (0.0, 0.0, 0.0)
    X_init.limit_upper = (Inf64, Inf64, Inf64)
    U_init = func.Marray(L, M, N, (m,), 0.0*ones(N,m))
    U_init.limit_lower = (0.0, 0.0, 0.0, 0.0)
    U_init.limit_upper = (Inf64, Inf64, Inf64, Inf64)

    X, U, K = HM_iLQR.RunHM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs, 0.01)

    return 

end

main()