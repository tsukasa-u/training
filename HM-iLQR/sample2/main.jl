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

    model = RobotZoo.Cartpole()
    n,m = RobotDynamics.dims(model)

    funcs = func.FUNC()
    funcs._h = 0.05
    funcs.f = (x, u) -> func.dynamics_rk4(x, u, (x, u) -> RobotZoo.dynamics(model, x, u), funcs._h)

    Q = Diagonal([1.0*ones(2); 1.0*ones(2)])
    R = 0.1*ones(1, 1)
    Qn = Array(100.0*I(n))
    xgoal = [0, pi, 0, 0]
    funcs.l = (x, u) -> 0.5*((x-xgoal)'*Q*(x-xgoal) + u'*R*u)
    funcs.lf = (x) -> 0.5*((x-xgoal)'*Qn*(x-xgoal))
    funcs.ge = (x, u) -> [0.0]
    funcs.glf = (x) -> [0.0]
    funcs.gl = (x, u) -> [0.0]
    funcs.glf = (x) -> [0.0]

    N = 100
    M = 1
    MaxIter = 1000
    ϵ_v = 1e-6
    d_max = 1e-6
    
    _, _, L = M_iLQR.init(N, M)
    X_init = func.Marray(L, M, N+1, (n,))
    U_init = func.Marray(L, M, N, (m,), 0.001*ones(N,m))

    X, U, K = HM_iLQR.RunHM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs)

    return 

end

main()