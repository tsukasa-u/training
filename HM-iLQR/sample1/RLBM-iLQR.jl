include("M-iLQR.jl")
include("func.jl")
using .func

module ALM_iLQR
    using LinearAlgebra
    using Plots
    using Main.M_iLQR
    include("func.jl")
    using .func

    function compute_J2(X, U, ψ, δ, funcs)
        a = (X, U) -> begin
            funcs.B(X, U, ψ, δ)
        end
        return M_iLQR.compute_J(X, U, funcs) + Main.func.sumMarray(a, X, U) + Main.func.endMarray(a, X, U)
    end

    M_iLQR.wrapCompute_J(X, U, ψ, δ, funcs) = compute_J2(X, U, ψ, δ, funcs)

    function backwardPass!(X, U, Vx, Vxx, d, funcs)
        k = Main.func.Marray(X._L, X.M, X.N, (U.n[1],))
        K = Main.func.Marray(X._L, X.M, X.N, (U.n[1], X.n[1]))
        ΔV = Main.func.Marray(X._L, X.M, X.N, ())
        for i in U.M:-1:1
            for j in U._L:-1:1
                Vx_ = Vx[i, j+1, :] + Vxx[i, j+1, :, :] * d[i, j]
                Qx = funcs.lx(X[i, j, :], U[i, j, :]) + funcs.fx(X[i, j, :], U[i, j, :])' * Vx_
                Qu = funcs.lu(X[i, j, :], U[i, j, :]) + funcs.fu(X[i, j, :], U[i, j, :])' * Vx_
                Qxx = funcs.lxx(X[i, j, :], U[i, j, :]) + funcs.fx(X[i, j, :], U[i, j, :])' * Vxx[i, j+1, :, :] * funcs.fx(X[i, j, :], U[i, j, :])
                Quu = funcs.luu(X[i, j, :], U[i, j, :]) + funcs.fu(X[i, j, :], U[i, j, :])' * Vxx[i, j+1, :, :] * funcs.fu(X[i, j, :], U[i, j, :])
                Qux = funcs.lux(X[i, j, :], U[i, j, :]) + funcs.fu(X[i, j, :], U[i, j, :])' * Vxx[i, j+1, :, :] * funcs.fx(X[i, j, :], U[i, j, :])
                
                k[i, j, :] = -Quu \ Qu
                K[i, j, :, :] = -Quu \ Qux
                ΔV[i, j] = 0.5 * Qu' * k[i, j, :]
                Vx[i, j, :] = Qx + vec(Qu' * K[i, j, :, :])
                Vxx[i, j, :, :] = Qxx + Qux' * K[i, j, :, :]
            end
        end
        return k, K, ΔV
    end

    M_iLQR.wrapbackwardPass!(X, U, Vx, Vxx, d, funcs) = backwardPass!(X, U, Vx, Vxx, d, funcs)

    function RunM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs)
        
        _J, J_, L = M_iLQR.init(N, M)

        X, U = M_iLQR.initFP(X_init, U_init, funcs)

        d = M_iLQR.conpute_d(X, U, funcs)
        _J = M_iLQR.wrapCompute_J(X, U, ψ, δ, funcs)
        println("J = ", _J)

        Vx, Vxx = M_iLQR.initV(X, funcs)

        anim = Animation()

        for idx in 1:MaxIter
            println("start step ", idx, ", J = ", _J)
            M_iLQR.updateV!(X, Vx, Vxx, funcs)
            k, K, ΔV = M_iLQR.wrapbackwardPass!(X, U, Vx, Vxx, d, funcs)

            X, U = M_iLQR.forwardPass(X, U, k, K, d, funcs)

            J_ = M_iLQR.wrapCompute_J(X, U, ψ, δ, funcs)
            ΔJ = J_ - _J
            
            d = M_iLQR.conpute_d(X, U, funcs)

            frame(
                anim, 
                plot(X[:, :, 1], X[:, :, 2], label="$(idx) : ΔJ = $(_J)")
                # plot(Main.func.getMarray(X)[:, 1], Main.func.getMarray(X)[:, 2], label="$(idx) : ΔJ = $(_J)")
                # plot(range(start = 0.0, step = 0.01, length = U._L), Main.func.getMarray(U)[:, 1], label="$(idx) : ΔJ = $(_J)")
            )

            if (abs(ΔJ) < ϵ_v && norm(d, 2) < d_max) || idx == MaxIter
            # if idx == MaxIter
                gif(anim, "ALM-iLQR_anim_fps15.gif", fps = 15)
                return X, U, K
            end

            _J = J_
        end
        return X, U, K
    end

    M_iLQR.wrapRunM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs) = RunM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs)

    function update_ψ_δ!(ψ, δ, ω₁, ω₂, δ_min)
        @assert ψ > 0 "ψ must be positive"
        @assert ω₁ < 1.0 "ω₁ must be less than 1.0"
        @assert δ_min > 0 "δ_min must be positive"
        @assert ω₂ < 1.0 "ω₂ must be less than 1.0"
        ψ = ω₁*ψ
        δ = max(δ_min, ω₂*δ)
    end

    function init()
        return 1.0, 0.01, 0.5, 0.5, 1.0, 1E-6
    end

    function RunBLBM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs)

        δ, δ_min, ω_1, ω_2, ψ, ϵ_g = init()
        
        X, U, K = M_iLQR.wrapRunM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs)
        update_ψ_δ!(ψ, δ, ω_1, ω_2, δ_min)

        n = size(funcs.g(X[1, 1, :], U[1, 1, :]))
        if n[1] > 0
            while max(norm.(funcs.g.(X, U), 2)[:, :]) > ϵ_g
                X, U, K = M_iLQR.wrapRunM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs)
                update_ψ_δ!(ψ, δ, ω_1, ω_2, δ_min)
            end
        end

        return X, U, K
    end
end

# using RobotZoo
# using RobotDynamics
# using LinearAlgebra
# using Plots

# function main()
#     model = RobotZoo.Cartpole()
#     n,m = RobotDynamics.dims(model)

#     funcs = func.FUNC()
#     funcs.a = model
#     funcs._h = 0.05

#     Q = Diagonal([1.0*ones(2); 1.0*ones(2)])
#     R = 0.1*ones(1, 1)
#     Qn = Array(100.0*I(n))
#     xgoal = [0, pi, 0, 0]
#     funcs.l = (x, u) -> 0.5*((x-xgoal)'*Q*(x-xgoal) + u'*R*u)
#     funcs.lf = (x) -> 0.5*((x-xgoal)'*Qn*(x-xgoal))
#     funcs.g = (x, u) -> []

#     N = 100
#     M = 1
#     MaxIter = 1000
#     ϵ_v = 1e-6
#     d_max = 1e-6
    
#     _, _, L = M_iLQR.init(N, M)
#     X_init = func.Marray(L, M, N+1, (n,))
#     U_init = func.Marray(L, M, N, (m,), 0.001*ones(N,m))

#     X, U, K = ALM_iLQR.RunALM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs)

# end

# main()