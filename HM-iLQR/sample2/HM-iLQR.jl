# include("func.jl")
include("ALM-iLQR.jl")
include("M-iLQR.jl")
include("RLBM-iLQR.jl")


module HM_iLQR
    using LinearAlgebra
    using Plots
    using Main.ALM_iLQR
    using Main.RLBM_iLQR
    using Main.func

    function RunHM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs)
        X_al, U_al, K_al = ALM_iLQR.RunALM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs)
        X_sol, U_sol, K_sol = RLBM_iLQR.RunRLBM_iLQR(X_al, U_al, K_al, N, M, MaxIter, ϵ_v, d_max, funcs)
        return X_sol, U_sol, K_sol
    end
end

using RobotZoo
using RobotDynamics
using LinearAlgebra
using Plots

function main()
    model = RobotZoo.Cartpole()
    n,m = RobotDynamics.dims(model)

    funcs = func.FUNC()
    funcs.a = model
    funcs._h = 0.05

    Q = Diagonal([1.0*ones(2); 1.0*ones(2)])
    R = 0.1*ones(1, 1)
    Qn = Array(100.0*I(n))
    xgoal = [0, pi, 0, 0]
    funcs.l = (x, u) -> 0.5*((x-xgoal)'*Q*(x-xgoal) + u'*R*u)
    funcs.lf = (x) -> 0.5*((x-xgoal)'*Qn*(x-xgoal))
    funcs.ge = (x, u) -> [0.0]
    funcs.gl = (x, u) -> [0.0]

    N = 100
    M = 1
    MaxIter = 1000
    ϵ_v = 1e-6
    d_max = 1e-6
    
    _, _, L = M_iLQR.init(N, M)
    X_init = func.Marray(L, M, N+1, (n,))
    U_init = func.Marray(L, M, N, (m,), 0.001*ones(N,m))

    X, U, K = HM_iLQR.RunHM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs)

end

main()