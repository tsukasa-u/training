include("func.jl")
# using .func


module M_iLQR
    using LinearAlgebra
    using Plots

    function init(N, M)
        return 0.0, 0.0, Int64((N - 1.0)/M + 1.0)
    end

    function initFP(X_init, U_init, funcs)
    # function initFP(M, N, X_init, U_init, funcs)
        X = copy(X_init)
        U = copy(U_init)
        for i in 1:U.M
            X[i, 1, :] = X_init[i, 1, :]
            for j in 1:U._L
                X[i, j+1, :] = funcs.f(X[i, j, :], U_init[i, j, :])
            end
        end
        # δx = X - X_init
        # return X, U, δx
        return X, U
    end

    function conpute_d(X, U, funcs)
    # function conpute_d(M, L, X, U, funcs)
        return [funcs.f(X[i, j, :], U[i, j, :]) - X[i, j+1, :] for i in 1:U.M, j in 1:U._L]
    end

    function computeJ(X, U, funcs)
    # function computeJ(M, L, N, X, U, funcs)
        # return sum(funcs.l(X[i, j, :], U[i, j, :]) for j in 1:L,  i in 1:M) + funcs.lf(X[i, N+1, :])
        return Main.func.sumMarray(funcs.l, X, U) + Main.func.endMarray(funcs.lf, X)
    end

    function backwardPass!(X, U, Vx, Vxx, d, funcs)
        k = Main.func.Marray(X._L, X.M, X.N, [U.n[1]])
        K = Main.func.Marray(X._L, X.M, X.N, [U.n[1], X.n[1]])
        ΔV = Main.func.Marray(X._L, X.M, X.N, [])
        for i in U.M:-1:1
            for j in U._L:-1:1
                # Qx = funcs.lx(X[i, j, :], U[i, j, :]) + funcs.fx(X[i, j, :], U[i, j, :])' * Vx[i, j+1, :]
                # Qu = funcs.lu(X[i, j, :], U[i, j, :]) + funcs.fu(X[i, j, :], U[i, j, :])' * Vx[i, j+1, :]
                Vx_ = Vx[i, j+1, :] + Vxx[i, j+1, :, :] * d[i, j]
                Qx = funcs.lx(X[i, j, :], U[i, j, :]) + funcs.fx(X[i, j, :], U[i, j, :])' * Vx_
                Qu = funcs.lu(X[i, j, :], U[i, j, :]) + funcs.fu(X[i, j, :], U[i, j, :])' * Vx_

                Qxx = funcs.lxx(X[i, j, :], U[i, j, :]) + funcs.fx(X[i, j, :], U[i, j, :])' * Vxx[i, j+1, :, :] * funcs.fx(X[i, j, :], U[i, j, :])
                Quu = funcs.luu(X[i, j, :], U[i, j, :]) + funcs.fu(X[i, j, :], U[i, j, :])' * Vxx[i, j+1, :, :] * funcs.fu(X[i, j, :], U[i, j, :])
                # println(size(funcs.fu(X[i, j, :], U[i, j, :])'))
                # println(size(Vxx[i, j+1, :, :]))
                # println(size(funcs.fx(X[i, j, :], U[i, j, :])))
                # println(size(funcs.lux(X[i, j, :], U[i, j, :])))
                # println(Vxx[i,j,:,:])
                Qux = funcs.lux(X[i, j, :], U[i, j, :]) + funcs.fu(X[i, j, :], U[i, j, :])' * Vxx[i, j+1, :, :] * funcs.fx(X[i, j, :], U[i, j, :])
                k[i, j, :] = -Quu \ Qu
                K[i, j, :, :] = -Quu \ Qux
                ΔV[i, j] = 0.5 * Qu' * k[i, j, :]
                # ΔV[i, j] = -0.5 * Qu' / Quu * Qu
                # println(size(Qx))
                # println(size(Qu'))
                # println(size(K[i, j, :, :]))
                # println(size(vec(Qu' * K[i, j, :, :])))
                Vx[i, j, :] = Qx + vec(Qu' * K[i, j, :, :])
                Vxx[i, j, :, :] = Qxx + Qux' * K[i, j, :, :]
                # println(Vxx[i,j,:,:])
                # println(Vxx[i, j, :, :] - Qxx - Qux' * K[i, j, :, :])
                # Vx[i, j, :] = Qx + vec(Qu' / Quu * Qux)
                # Vxx[i, j, :, :] = Qxx + Qux' / Quu * Qux
            end
        end
        return k, K, ΔV
    end

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

                # du = α * (k[i, j, :] + K[i, j, :, :] * _dx[i, j, :])
                # U[i, j, :] = _U[i, j, :] + du
                # dx[i, j+1, :] = funcs.fx(_X[i, j, :], _U[i, j, :]) * _dx[i, j, :] + funcs.fu(_X[i, j, :], _U[i, j, :]) * du + d[i, j]
                # X[i, j+1, :] = _X[i, j+1, :] + dx[i, j+1, :]

            end
        end
        # return X, U, dx
        return X, U
    end

    function initV(X, funcs)
    # function initV(M, L, X, U, funcs)
        # return [funcs.fx(X[i, j, :], U[i, j, :]) for i in 1:M, j in 1:L], [funcs.lfxx(X[i, j, :]) for i in 1:M, j in 1:L]
        vx = Main.func.Marray(X._L, X.M, X.N, [X.n[1]])
        Main.func.setEndMarray!(vx, funcs.lfx(Main.func.getEndMarray(X, [:])), [:])
        vxx = Main.func.Marray(X._L, X.M, X.N, [X.n[1], X.n[1]])
        Main.func.setEndMarray!(vxx, funcs.lfxx(Main.func.getEndMarray(X, [:])), [:, :])
        return vx, vxx
    end

    function updateV!(X, vx, vxx, funcs)
        Main.func.setEndMarray!(vx, funcs.lfx(Main.func.getEndMarray(X, [:])), [:])
        Main.func.setEndMarray!(vxx, funcs.lfxx(Main.func.getEndMarray(X, [:])), [:, :])
    end

    function RunM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs)
        
        _J, J_, L = init(N, M)

        X, U = initFP(X_init, U_init, funcs)
        # X, U, δx = initFP(M, N, X_init, U_init, funcs)

        # d = conpute_d(M, N, X, U, funcs)
        d = conpute_d(X, U, funcs)
        _J = computeJ(X, U, funcs)
        # _J = computeJ(M, L, N, X, U, funcs)
        println("J = ", _J)

        Vx, Vxx = initV(X, funcs)

        anim = Animation()

        for idx in 1:MaxIter
            println("start step ", idx, ", J = ", _J)
            updateV!(X, Vx, Vxx, funcs)
            k, K, ΔV = backwardPass!(X, U, Vx, Vxx, d, funcs)

            # α = 0.0
            X_, U_= forwardPass(X, U, k, K, d, funcs)
            # X, U, δx = forwardPass(X, U, k, K, d, funcs, α, δx)

            # for i in 1:30
            #     d = conpute_d(M, N, X_, U_, funcs)
            #     J_ = computeJ(M, L, N, X_, U_, funcs)
            #     if J_ < _J
            #         println("finate J = ", J_, ", α = ", α)
            #         break
            #     end
            #     println("J = ", J_, ", α = ", α)
            #     α = α * 0.8
            #     X_, U_, δx_ = forwardPass(X, U, k, K, d, funcs, α, δx)
            # end

            X = copy(X_)
            U = copy(U_)
            # δx = copy(δx_)

            J_ = computeJ(X, U, funcs)
            # J_ = computeJ(M, L, N, X, U, funcs)
            ΔJ = J_ - _J
            
            d = conpute_d(X, U, funcs)
            # d = conpute_d(M, N, X, U, funcs)

            frame(
                anim, 
                plot(Main.func.getMarray(X)[:, 1], Main.func.getMarray(X)[:, 2], label="$(idx) : ΔJ = $(_J)")
                # plot(range(start = 0.0, step = 0.01, length = U._L), Main.func.getMarray(U)[:, 1], label="$(idx) : ΔJ = $(_J)")
            )

            if (abs(ΔJ) < ϵ_v && norm(d, 2) < d_max) || idx == MaxIter
            # if idx == MaxIter
                gif(anim, "anim_fps15.gif", fps = 15)
                return X, U, K
            end

            _J = J_
        end
        return X, U, K
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
    funcs.h = 0.05

    Q = Diagonal([1.0*ones(2); 1.0*ones(2)])
    R = 0.1*ones(1, 1)
    Qn = Array(100.0*I(n))
    xgoal = [0, pi, 0, 0]
    funcs.l = (x, u) -> 0.5*((x-xgoal)'*Q*(x-xgoal) + u'*R*u)
    funcs.lf = (x) -> 0.5*((x-xgoal)'*Qn*(x-xgoal))

    N = 100
    M = 1
    MaxIter = 1000
    ϵ_v = 1e-6
    d_max = 1e-6

    # X_init = zeros(M, N+1, n)
    # U_init = zeros(M, N, m)
    
    _, _, L = M_iLQR.init(N, M)
    X_init = func.Marray(L, M, N+1, [n])
    U_init = func.Marray(L, M, N, [m], 0.001*ones(N,m))

    # for i in 1:M
    #     X_init[i, 1, :] = [0.0, 0.0, 0.0, 0.0]
    #     for j in 1:N
    #         U_init[i, j, :] = [0.0]
    #     end
    # end

    X, U, K = M_iLQR.RunM_iLQR(X_init, U_init, N, M, MaxIter, ϵ_v, d_max, funcs)

end

main()