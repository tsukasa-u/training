include("M-iLQR.jl")

module ALM_iLQR
    function init()

    end

    function compute_J1(X, U, λ, λf, μf, h, hf, funcs::FUNC)
        return M_iLQR.computeJ(X, U, funcs) + sum(λ'*h+0.5*μ*norm(h)^2 for j in 1:N) + λf'*hf + 0.5*μf*norm(hf)^2
    end

    function compute_h(M, L, funcs::FUNC)
        return [max(0, g(X[i, j, :], U[i, j, :])) for i in 1:M, j in 1:L, g in funcs.g]
    end

    function RunALM_iLQR()

    end
end