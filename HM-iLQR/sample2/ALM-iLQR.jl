
module ALM_iLQR
    using LinearAlgebra
    using Plots
    include("M-iLQR.jl")
    using .M_iLQR
    include("func.jl")
    using .func

    function compute_J1(X, U, μ, μf, λ, λf, funcs)
        if λ.n[1] == 0
            return M_iLQR.compute_J(X, U, funcs)
        else
            a = (X, U, λ, μ) -> begin
                h = funcs.h(X, U)
                return λ'*h + 0.5*μ*norm(h)^2
            end
            b = (X) -> begin
                hf = funcs.hf(X)
                return λf'*hf + 0.5*μf*norm(hf)^2
            end
            return M_iLQR.compute_J(X, U, funcs) + Main.func.sumMarray(a, X, U, λ, μ) + Main.func.endMarray(b, X)
        end
    end

    M_iLQR.wrapCompute_J(X, U, μ, μf, λ, λf, funcs) = compute_J1(X, U, μ, μf, λ, λf, funcs)

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
                    Qx += funcs.hx(X[i, j, :], U[i, j, :])'*(λ[i, j, :] + μ[i, j]*funcs.h(X[i, j, :], U[i, j, :]))
                    Qu += funcs.hu(X[i, j, :], U[i, j, :])'*(λ[i, j, :] + μ[i, j]*funcs.h(X[i, j, :], U[i, j, :]))
                    Qxx += μ[i, j]*funcs.hx(X[i, j, :], U[i, j, :])'*funcs.hx(X[i, j, :], U[i, j, :])
                    Quu += μ[i, j]*funcs.hu(X[i, j, :], U[i, j, :])'*funcs.hu(X[i, j, :], U[i, j, :])
                    Qux += μ[i, j]*funcs.hu(X[i, j, :], U[i, j, :])'*funcs.hx(X[i, j, :], U[i, j, :])    
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

    function RunM_iLQR(X_init, U_init, μ, μf, λ, λf, N, M, MaxIter, ϵ_v, d_max, funcs)
        
        _J, J_, L = M_iLQR.init(N, M)

        X, U = M_iLQR.initFP(X_init, U_init, funcs)

        d = M_iLQR.conpute_d(X, U, funcs)
        _J = M_iLQR.wrapCompute_J(X, U, μ, μf, λ, λf, funcs)
        println("J = ", _J)

        Vx, Vxx = M_iLQR.initV(X, funcs)

        anim = Animation()

        for idx in 1:MaxIter
            println("start step ", idx, ", J = ", _J)
            M_iLQR.updateV!(X, Vx, Vxx, funcs)
            k, K, ΔV = M_iLQR.wrapbackwardPass!(X, U, λ, μ, Vx, Vxx, d, funcs)

            X, U = M_iLQR.forwardPass(X, U, k, K, d, funcs)

            J_ = M_iLQR.wrapCompute_J(X, U, μ, μf, λ, λf, funcs)
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

    M_iLQR.wrapRunM_iLQR(X_init, U_init, μ, μf, λ, λf, N, M, MaxIter, ϵ_v, d_max, funcs) = RunM_iLQR(X_init, U_init, μ, μf, λ, λf, N, M, MaxIter, ϵ_v, d_max, funcs)

    function update_λ_μ!(λ, λf, μ, μf, X, U, funcs, φ)
        @assert φ > 1.0 "φ must be greater than 1.0"
        if λ.n[1] > 0
            λ = ((λ, μ, X, U) -> max.(0, λ + μ*[funcs.gl(X, U); funcs.ge(X, U)])).(λ, μ, X, U)
            
            λf = max.(0, λf + μf*[Main.func.endMarray(funcs.glf, X); Main.func.endMarray(funcs.gef, X)])
            
            μ = φ*μ
            μf = φ*μf
        end
    end

    function init(_L, M, N, n, _n)
        return Main.func.onesMarray(_L, M, N, n), ones(_n...), Main.func.onesMarray(_L, M, N, ()), 1.0, 2.0, 1E-6
    end

    function RunALM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs)

        λ, λf, μ, μf, φ, ϵ_g = init(U_init._L, U_init.M, U_init.N, size(funcs.h(X_init[1, 1, :], U_init[1, 1, :])), size(funcs.hf(X_init[1, 1, :])))
        
        X, U, K = M_iLQR.wrapRunM_iLQR(X_init, U_init, μ, μf, λ, λf, N, M, MaxIter, ϵ_v, d_max, funcs)
        update_λ_μ!(λ, λf, μ, μf, X, U, funcs, φ)

        # if λ.n[1] > 0
            while M_iLQR.compute_max_g(X, U, funcs) > ϵ_g
            # while max(norm.(funcs.g.(X, U))[:, :]...) > ϵ_g
                X, U, K = M_iLQR.wrapRunM_iLQR(X_init, U_init, μ, μf, λ, λf, N, M, MaxIter, ϵ_v, d_max, funcs)
                update_λ_μ!(λ, λf, μ, μf, X, U, funcs, φ)
            end
        # end

        return X, U, K
    end
end