include("func.jl")

module M_iLQR

    function init(N, M)
        return 0.0, 0,.0 (N - 1.0)/M + 1.0
    end

    function initFP(M::S, N::S, X_init, U_init, funcs::FUNC) where S <: Int
        X = zeros(size(X_init))
        U = copy(U_init)
        L::S = (N - 1)/M + 1
        for i in 1:M
            X[i, 1, :] = X_init[i, 1, :]
            for j in 1:L
                X[i, j+1, :] = funcs.f(X[i, j, :], U_init[i, j, :])
            end
        end
        return X, U
    end

    function conpute_d(X, U, funcs::FUNC)
        return [funcs.f(X[i, j, :], U[i, j, :]) - X[i, j+1, :] for i in 1:M, j in 1:L]
    end

    function computeJ(X, U, funcs::FUNC)
        return [sum(funcs.l(X[i, j, :], U[i, j, :]) for j in 1:L) + funcs.lf(X[i, L+1, :]) for i in 1:M]
    end

    wrap_BP = () -> ()
    wrap_FP = () -> ()

    function backwardPass!(Vx, Vxx)
        Qx = funcs.lx(X[i, j, :], u[i, j, :]) + funcs.fx(x[i, j, :], u[i, j, :])' * Vx[i, j+1, :]
        Qu = funcs.lu(X[i, j, :], u[i, j, :]) + funcs.fu(x[i, j, :], u[i, j, :])' * Vx[i, j+1, :]
        Qxx = funcs.lxx(X[i, j, :], u[i, j, :]) + funcs.fx(x[i, j, :], u[i, j, :])' * Vxx[i, j+1, :] * funcs.fx(x[i, j, :], u[i, j, :])
        Quu = funcs.luu(X[i, j, :], u[i, j, :]) + funcs.fu(x[i, j, :], u[i, j, :])' * Vxx[i, j+1, :] * funcs.fu(x[i, j, :], u[i, j, :])
        Qux = funcs.lux(X[i, j, :], u[i, j, :]) + funcs.fu(x[i, j, :], u[i, j, :])' * Vxx[i, j+1, :] * funcs.fx(x[i, j, :], u[i, j, :])
        k[i, j, :] = -Quu \ Qu
        K[i, j, :] = -Quu \ Qux
        ΔV[i, j, :] = -0.5 * Qu[i, j, :]*k[i, j, :]
        Vx[i, j, :] = Qx[i, j, :] + Qu[i, j, :] * K[i, j, :]
        Vxx[i, j, :] = Qxx[i, j, :] + Qux[i, j, :]' * K[i, j, :]
        return k, K, ΔV
    end

    function forwardPass!(X, U, k, K, funcs::FUNC)
        X_ = copy(X[:, :, :])
        for i in 1:M
            for j in 1:L
                U[i, j, :] = U[i, j, :] + k[i, j, :] + K[i, j, :] * (X[i, j, :]- X_[i, j, :])
                X[i, j+1, :] = funcs.f(X[i, j, :], U[i, j, :])
            end
        end
        return X, U
    end

    function initV(M, L, X, U, funcs::FUNC)
        return [funcs.fx(X[i, j, :], U[i, j, :]) for i in 1:M, j in 1:L], [funcs.lfxx(X[i, j, :], U[i, j, :]) for i in 1:M, j in 1:L]
    end

    function RunM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs::FUNC)
        
        _J, J_, L = init(N, M)

        X, U = initFP(M, N, X_init, U_init, funcs)

        _d = conpute_d(X, U, f)
        _J = computeJ(X, U, funcs)

        Vx, Vxx = initV(M, L, X, U, funcs)

        for _ in 1:MaxIter
            k, K, ΔV = backwardPass!(Vx, Vxx)
            X, U = forwardPass!(X, U, k, K, f)

            d_ = conpute_d(X, U, f)
            J_ = computeJ(X, U, funcs)
            ΔJ = J_ - _J

            if ΔJ < ϵ_v && norm(d_, 2) < d_max
                return X, U, K
            end

            _J = J_
        end
    end
    return X, U, K
end