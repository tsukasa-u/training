module HM_iLQR
    using LinearAlgebra
    using Plots
    include("ALM-iLQR.jl")
    using .ALM_iLQR
    include("RLBM-iLQR.jl")
    using .RLBM_iLQR

    function RunHM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs, α)
        X_al, U_al, K_al = ALM_iLQR.RunALM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs, α)
        X_sol, U_sol, K_sol = RLBM_iLQR.RunRLBM_iLQR(X_al, U_al, K_al, N, M, MaxIter, ϵ_v, d_max, funcs)
        return X_sol, U_sol, K_sol
    end
end
