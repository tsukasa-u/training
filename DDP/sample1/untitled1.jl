using LinearAlgebra
using ForwardDiff

n = 10

f(x, u) = x
fx(x, u) = ForwardDiff.gradient(x -> f(x, u), x)
fu(x, u) = ForwardDiff.gradient(x -> f(x, u), u)
fxx(x, u) = ForwardDiff.jacobian(x -> fx(x, u), x)
fux(x, u) = ForwardDiff.jacobian(u -> fx(x, u), u)
fuu(x, u) = ForwardDiff.jacobian(u -> fu(x, u), u)

l(x, u) = 1.0
lx(x, u) = ForwardDiff.gradient(x -> l(x, u), x)
lu(x, u) = ForwardDiff.gradient(u -> l(x, u), u)
lxx(x, u) = ForwardDiff.jacobian(x -> lx(x, u), x)
lux(x, u) = ForwardDiff.jacobian(u -> lx(x, u), u)
luu(x, u) = ForwardDiff.jacobian(u -> lu(x, u), u)

lf(x) = 0.0

v(i, x, u) = sum(l(_x, _u) for (_x, _u) in zip(x[i:end], u[i:end])) + lf(x[end])
vx(i, x, u) = ForwardDiff.gradient(x -> v(i, x, u), x)
vxx(i, x, u) = ForwardDiff.jacobian(x -> vx(i, x, u), x)

function compute_Q(i, x, u)

    _lx = lx(x, u)
    _lu = lu(x, u)
    _lxx = lxx(x, u)
    _lux = lux(x, u)
    _luu = luu(x, u)
    _vx = vx(i + 1, x, u)
    _vxx = vxx(i + 1, x, u)
    _fx = fx(x, u)
    _fu = fu(x, u)
    _fxx = fxx(x, u)
    _fux = fux(x, u)
    _fuu = fuu(x, u)

    Qx = _lx + transpose(_fx) * _vx
    Qu = _lu + transpose(_fu) * _vx
    Qxx = _lxx + transpose(_fx) * _vxx * _fx + kron(_vx, _fxx)
    Qux = _lux + transpose(_fu) * _vxx * _fx + kron(_vx, _fux)
    Quu = _luu + transpose(_fu) * _vxx * _fu + kron(_vx, _fuu)

    return Qx, Qu, Qxx, Qux, Quu
end

function normalize_inverse_Quu(Quu)
    return inv(normalize(Quu, 1e-6))
end

function compute_k(Qxx_1, Qu)
    return -Qxx_1 * Qu
end

function compute_K(Quu_1, Qux)
    return -Quu_1 * Qux
end

function BFP(n, x, u)

    k = []
    K = []

    for i in (n-1):-1:1
        Qx, Qu, Qxx, Qux, Quu = compute_Q(i, x, u)
        Quu_1 = normalize_inverse_Quu(Quu)
        push!(k, compute_k(Quu_1, Qu))
        push!(K, compute_K(Quu_1, Qux))
    end

    return k, K
end

function FFP(n, x, u, α)

    _x = copy(x)
    _u = copy(u)
    for i in 1:(n-1)
        _u[i] = u[i] + α*k[i] + K[i] * (_x[i] - x[i])
        _x[i + 1] = f(_x[i], _u[i])
    end
    
    return _x, _u
end

function BFFP(n, x, u, α)

    BFP(n, x, u)
    _x, u = FFP(n, x, u, α)

    ret = max(abs.(_x. - x)) < 1e-6

    x = copy(_x)

    return ret
end

function loop()
    
        x = [0.0 for i in 1:n]
        u = [0.0 for i in 1:n]
    
        for i in 1:100
            BFFP(n, x, u, 0.8)
            println(x)
        end
end

loop()