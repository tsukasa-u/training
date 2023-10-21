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
                v_θ^2/r - 1/r^2 + u_r_p - u_r_m,
                -v_r*v_θ/r + u_θ_p - u_θ_m
            ]
        end, 
        funcs._h
    )

    n = 3
    m = 4
    Q = Diagonal(0.1*ones(n))
    R = Diagonal(ones(m))
    R_ = ones(m)
    Qn = Array(100.0*I(n))
    xgoal = [4, 0.5, 0]
    # funcs.l = (x, u) -> 0.5*((x-xgoal)'*Q*(x-xgoal) + u'*R*u)
    funcs.lf = (x) -> 0.5*((x-xgoal)'*Qn*(x-xgoal))
    funcs.l = (x, u) -> R_'*u
    # funcs.lf = (x) -> 0.0
    funcs.ge = (x, u) -> begin
        r, v_r, v_θ = x
        u_r_p, u_r_m, u_θ_p, u_θ_m = u
        return [ u_r_p*u_r_m, u_θ_p*u_θ_m]
    end
    funcs.gef = (x) -> [0.0]
    funcs.gl = (x, u) -> begin
        r, v_r, v_θ = x
        u_r_p, u_r_m, u_θ_p, u_θ_m = u
        return [-r, u_r_p - 0.01, u_r_m - 0.01, u_θ_p - 0.01, u_θ_m - 0.01, -u_r_p, -u_r_m, -u_θ_p, -u_θ_m]
    end
    funcs.glf = (x) -> [0.0]

    N = 100
    M = 1
    MaxIter = 200
    ϵ_v = 1e-3
    d_max = 1e-3
    
    _, _, L = M_iLQR.init(N, M)
    X_init = func.Marray(L, M, N+1, (n,), [ele for _ in 1:N+1, ele in [1.0, 1.0, 0.0]])
    U_init = func.Marray(L, M, N, (m,), 0.0*ones(N,m))

    X, U, K = HM_iLQR.RunHM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs)

    return 

end

main()