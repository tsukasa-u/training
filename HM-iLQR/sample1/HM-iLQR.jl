
module ALM_iLQR
    using LinearAlgebra
    using Plots
    using Main.ALM_iLQR
    using Main.RLBM_iLQR
    include("func.jl")
    using .func

    function RunHM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs)
        X_al, U_al, K_al = M_iLQR.RunM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs)
        X_sol, U_sol, K_sol = RLBM_iLQR.RunRLBM_iLQR(X_al, U_al, N, M, MaxIter, ϵ_v, d_max, funcs)
        return X_sol, U_sol, K_sol
    end
end