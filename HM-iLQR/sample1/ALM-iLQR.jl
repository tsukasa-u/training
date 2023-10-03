include("M-iLQR.jl")
include("func.jl")
using .func

module ALM_iLQR
    using LinearAlgebra
    using Plots
    using Main.M_iLQR
    include("func.jl")
    using .func

    function compute_J1(X, U, λ, μ, funcs)
        # return Main.func.sumMarray(a, X, U, λ, μ, h) + Main.func.endMarray(a, X, λ, μ, h)
        # println(λ.n)
        if λ.n[1] == 0
            return M_iLQR.compute_J(X, U, funcs)
        else
            a = (X, U, λ, μ) -> begin
                λ'*funcs.h(X, U) + 0.5*μ*norm(funcs.h(X, U))^2
            end
            return M_iLQR.compute_J(X, U, funcs) + Main.func.sumMarray(a, X, U, λ, μ) + Main.func.endMarray(a, X, U, λ, μ)
        end
    end

    M_iLQR.wrapCompute_J(X, U, λ, μ, funcs) = compute_J1(X, U, λ, μ, funcs)

    # function compute_h(M, L, funcs)
    #     return [max(0, g(X[i, j, :], U[i, j, :])) for i in 1:M, j in 1:L, g in funcs.g]
    # end

    function backwardPass!(X, U, λ, μ, Vx, Vxx, d, funcs)
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
                
                if λ.n[1] > 0
                    Qx += funcs.hx(X[i, j, :], U[i, j, :])'*(λ[i, j, :] + μ[i, j, :]*funcs.h(X[i, j, :], U[i, j, :]))
                    Qu += funcs.hu(X[i, j, :], U[i, j, :])'*(λ[i, j, :] + μ[i, j, :]*funcs.h(X[i, j, :], U[i, j, :]))
                    Qxx += μ[i, j, :]*funcs.hx(X[i, j, :], U[i, j, :])'*funcs.hx(X[i, j, :], U[i, j, :])
                    Quu += μ[i, j, :]*funcs.hu(X[i, j, :], U[i, j, :])'*funcs.hu(X[i, j, :], U[i, j, :])
                    Qux += μ[i, j, :]*funcs.hu(X[i, j, :], U[i, j, :])'*funcs.hx(X[i, j, :], U[i, j, :])    
                end
                
                k[i, j, :] = -Quu \ Qu
                K[i, j, :, :] = -Quu \ Qux
                ΔV[i, j] = 0.5 * Qu' * k[i, j, :]
                Vx[i, j, :] = Qx + vec(Qu' * K[i, j, :, :])
                Vxx[i, j, :, :] = Qxx + Qux' * K[i, j, :, :]
            end
        end
        return k, K, ΔV
    end

    M_iLQR.wrapbackwardPass!(X, U, λ, μ, Vx, Vxx, d, funcs) = backwardPass!(X, U, λ, μ, Vx, Vxx, d, funcs)

    function RunM_iLQR(X_init, U_init, μ, λ, N, M, MaxIter, ϵ_v, d_max, funcs)
        
        _J, J_, L = M_iLQR.init(N, M)

        X, U = M_iLQR.initFP(X_init, U_init, funcs)

        d = M_iLQR.conpute_d(X, U, funcs)
        _J = M_iLQR.wrapCompute_J(X, U, λ, μ, funcs)
        println("J = ", _J)

        Vx, Vxx = M_iLQR.initV(X, funcs)

        anim = Animation()

        for idx in 1:MaxIter
            println("start step ", idx, ", J = ", _J)
            M_iLQR.updateV!(X, Vx, Vxx, funcs)
            k, K, ΔV = M_iLQR.wrapbackwardPass!(X, U, λ, μ, Vx, Vxx, d, funcs)

            X, U = M_iLQR.forwardPass(X, U, k, K, d, funcs)

            J_ = M_iLQR.wrapCompute_J(X, U, λ, μ, funcs)
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

    M_iLQR.wrapRunM_iLQR(X_init, U_init, μ, λ, N, M, MaxIter, ϵ_v, d_max, funcs) = RunM_iLQR(X_init, U_init, μ, λ, N, M, MaxIter, ϵ_v, d_max, funcs)

    function update_λ!(λ, μ, X, U, funcs, φ)
        @assert φ > 1.0 "φ must be greater than 1.0" 
        println(λ.n)
        if λ.n[1] > 0
            println(λ.n)
            copyto!(λ, max.(0, λ + μ.*funcs.g.(X, U)))
            copyto!(μ, φ.*μ)
        end
    end

    function init(_L, M, N, n)
        return Main.func.onesMarray(0_L, M, N, n), Main.func.onesMarray(_L, M, N, n), 2.0, 1E-6
    end

    function RunALM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs)

        λ, μ, φ, ϵ_g = init(X_init._L, X_init.M, X_init.N, size(funcs.g(X_init[1, 1, :], U_init[1, 1, :])))
        
        X, U, K = M_iLQR.wrapRunM_iLQR(X_init, U_init, μ, λ, N, M, MaxIter, ϵ_v, d_max, funcs)
        update_λ!(λ, μ, X, U, funcs, φ)

        if λ.n[1] > 0
            while max(norm.(funcs.g.(X, U), 2)[:, :]) > ϵ_g
                X, U, K = M_iLQR.wrapRunM_iLQR(X_init, U_init, μ, λ, N, M, MaxIter, ϵ_v, d_max, funcs)
                update_λ!(λ, μ, X, U, funcs, φ)
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