
module M_iLQR
    using LinearAlgebra
    using Plots
    include("func.jl")
    using .func

    function init(N, M)
        return 0.0, 0.0, Int64((N - 1.0)/M + 1.0)
    end

    function initFP(X_init, U_init, funcs)
        X = copy(X_init)
        U = copy(U_init)
        for i in 1:U.M
            X[i, 1, :] = X_init[i, 1, :]
            for j in 1:U._L
                X[i, j+1, :] = funcs.f(X[i, j, :], U_init[i, j, :])
            end
        end
        
        plot(X[:, :, 1], label="0")
        png("M-iLQR_init.png")
        return X, U
    end

    function conpute_d(X, U, funcs)
        return [funcs.f(X[i, j, :], U[i, j, :]) - X[i, j+1, :] for i in 1:U.M, j in 1:U._L]
    end

    function compute_J(X, U, funcs)
        return Main.func.sumMarray(funcs.l, X, U) + Main.func.endMarray(funcs.lf, X)
    end

    wrapCompute_J(X, U, funcs) = compute_J(X, U, funcs)

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

    wrapbackwardPass!(X, U, Vx, Vxx, d, funcs) = backwardPass!(X, U, Vx, Vxx, d, funcs)

    function forwardPass(_X, _U, k, K, d, funcs)
    # function forwardPass(_X, _U, k, K, d, funcs, α, _dx)
        X = copy(_X)
        U = copy(_U)
        # dx = copy(_dx)
        for i in 1:U.M
            for j in (1:U._L)
                __dx = X[i, j, :] - _X[i, j, :]
                __du = k[i, j, :] + K[i, j, :, :] * __dx
                U[i, j, :] = _U[i, j, :] + __du
                X[i, j+1, :] = _X[i, j+1, :] + funcs.fx(_X[i, j, :], _U[i, j, :]) * __dx + funcs.fu(_X[i, j, :], _U[i, j, :]) * __du + d[i, j]

            end
        end
        return X, U
    end

    function initV(X, funcs)
        vx = Main.func.Marray(X._L, X.M, X.N, (X.n[1],))
        Main.func.setEndMarray!(vx, funcs.lfx(Main.func.getEndMarray(X, [:])), [:])
        vxx = Main.func.Marray(X._L, X.M, X.N, (X.n[1], X.n[1]))
        Main.func.setEndMarray!(vxx, funcs.lfxx(Main.func.getEndMarray(X, [:])), [:, :])
        return vx, vxx
    end

    function updateV!(X, vx, vxx, funcs)
        Main.func.setEndMarray!(vx, funcs.lfx(Main.func.getEndMarray(X, [:])), [:])
        Main.func.setEndMarray!(vxx, funcs.lfxx(Main.func.getEndMarray(X, [:])), [:, :])
    end

    function compute_max_g(X, U, funcs)
        # return max(norm.(funcs.gl.(X, U))[:, :]...)
        return max(max(norm.(funcs.h.(X, U))[:, :]...), max(norm.(funcs.hf.(X))[:, :]...))
    end

    function RunM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs)
        
        _J, J_, L = init(N, M)

        X, U = initFP(X_init, U_init, funcs)

        d = conpute_d(X, U, funcs)
        _J = wrapCompute_J(X, U, funcs)
        println("J = ", _J)

        Vx, Vxx = initV(X, funcs)

        anim = Animation()

        for idx in 1:MaxIter
            println("start step ", idx, ", J = ", _J)
            updateV!(X, Vx, Vxx, funcs)
            k, K, ΔV = wrapbackwardPass!(X, U, Vx, Vxx, d, funcs)

            X, U = forwardPass(X, U, k, K, d, funcs)

            J_ = wrapCompute_J(X, U, funcs)
            ΔJ = J_ - _J
            
            d = conpute_d(X, U, funcs)

            frame(
                anim, 
                plot(X[:, :, 1], X[:, :, 2], label="$(idx) : ΔJ = $(_J)")
                # plot(Main.func.getMarray(X)[:, 1], Main.func.getMarray(X)[:, 2], label="$(idx) : ΔJ = $(_J)")
                # plot(range(start = 0.0, step = 0.01, length = U._L), Main.func.getMarray(U)[:, 1], label="$(idx) : ΔJ = $(_J)")
            )

            if (abs(ΔJ) < ϵ_v && norm(d, 2) < d_max) || idx == MaxIter
            # if idx == MaxIter
                gif(anim, "M-iLQR_anim_fps15.gif", fps = 15)
                return X, U, K
            end

            _J = J_
        end
        return X, U, K
    end

    wrapRunM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs) = RunM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs)
end