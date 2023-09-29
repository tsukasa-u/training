include("func.jl")
using .func


module M_iLQR

    function init(N, M)
        return 0.0, 0.0, Int64((N - 1.0)/M + 1.0)
    end

    function initFP(M::S, N::S, X_init, U_init, funcs) where S <: Int
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

    function conpute_d(M, L, X, U, funcs)
        return [funcs.f(X[i, j, :], U[i, j, :]) - X[i, j+1, :] for i in 1:M, j in 1:L]
    end

    function computeJ(M, L, N, X, U, funcs)
        return sum(funcs.l(X[i, j, :], U[i, j, :]) for j in 1:L,  i in 1:M) + funcs.lf(X[i, N+1, :])
    end

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

    function forwardPass!(X, U, k, K, d, funcs)
        X_ = copy(X[:, :, :])
        for i in 1:M
            for j in 1:L
                dx = X[i, j, :]- X_[i, j, :]
                du = k[i, j, :] + K[i, j, :] * dx
                U[i, j, :] = U[i, j, :] + du
                X[i, j+1, :] = X[i, j + 1, :] + funcs.fx(X[i, j, :], U[i, j, :]) * dx + funcs.fu(X[i, j, :], U[i, j, :]) * du + d[i, j]
            end
        end
        return X, U
    end

    function initV(M, L, X, U, funcs)
        return [funcs.fx(X[i, j, :], U[i, j, :]) for i in 1:M, j in 1:L], [funcs.lfxx(X[i, j, :], U[i, j, :]) for i in 1:M, j in 1:L]
    end

    function RunM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs)
        
        _J, J_, L = init(N, M)

        X, U = initFP(M, N, X_init, U_init, funcs)

        d = conpute_d(M, N, X, U, funcs)
        _J = computeJ(M, L, N, X, U, funcs)

        Vx, Vxx = initV(M, L, X, U, funcs)

        for _ in 1:MaxIter
            k, K, ΔV = backwardPass!(Vx, Vxx)
            X, U = forwardPass!(X, U, k, K, d, funcs)

            d = conpute_d(M, N, X, U, f)
            J_ = computeJ(M, L, N, X, U, funcs)
            ΔJ = J_ - _J

            if ΔJ < ϵ_v && norm(d, 2) < d_max
                return X, U, K
            end

            _J = J_
        end
        return X, U, K
    end
end

using RobotZoo
using RobotDynamics
using Setfield

function main()
    model = RobotZoo.Cartpole()
    n,m = RobotDynamics.dims(model)

    funcs = func.FUNC()
    funcs.a = model
    funcs.h = 0.01

    N = 100
    M = 1
    MaxIter = 100
    ϵ_v = 1e-6
    d_max = 1e-6

    X_init = zeros(M, N+1, n)
    U_init = zeros(M, N, m)

    # for i in 1:M
    #     X_init[i, 1, :] = [0.0, 0.0, 0.0, 0.0]
    #     for j in 1:N
    #         U_init[i, j, :] = [0.0]
    #     end
    # end

    M_iLQR.RunM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs)
end

main()