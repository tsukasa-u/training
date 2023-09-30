module func
    using RobotZoo
    using RobotDynamics
    using ForwardDiff

    export FUNC

    mutable struct FUNC

        f::Function
        fx::Function
        fu::Function
        
        l::Function
        lx::Function
        lu::Function
        lxx::Function
        luu::Function
        lux::Function

        lf::Function
        lfx::Function
        lfxx::Function

        g::Function
        h::Function
        hx::Function
        hu::Function

        a::RobotZoo.Cartpole{Float64}
        _h::Float64
        M::Int64
        N::Int64

        function FUNC()
            _new = new()

            _new.a = RobotZoo.Cartpole()
            _new._h = 0.01

            _new.f = (x, u) -> dynamics_rk4(x,u,_new.a,_new.h)
            _new.fx = (x, u) -> ForwardDiff.jacobian(dx->_new.f(dx,u),x)
            _new.fu = (x, u) -> ForwardDiff.jacobian(du->_new.f(x,du),u)

            _new.l = (x, u) -> @assert false "l is not defined"
            _new.lx = (x, u) -> ForwardDiff.gradient(dx->_new.l(dx,u),x)
            _new.lu = (x, u) -> ForwardDiff.gradient(du->_new.l(x,du),u)
            _new.lxx = (x, u) -> ForwardDiff.hessian(dx->_new.l(dx,u),x)
            _new.luu = (x, u) -> ForwardDiff.hessian(du->_new.l(x,du),u)
            # _new.lux = (x, u) -> ForwardDiff.jacobian(du->ForwardDiff.gradient(dx->_new.l(dx,du),x),u)
            _new.lux = (x, u) -> ForwardDiff.jacobian(dx->_new.lu(dx,u),x)

            _new.lf = (x) -> @assert false "lf is not defined"
            _new.lfx = (x) -> ForwardDiff.gradient(dx->_new.lf(dx),x)
            _new.lfxx = (x) -> ForwardDiff.hessian(dx->_new.lf(dx),x)

            _new.g = (x, u) -> @assert false "g is not defined"
            _new.h = (x, u) -> max.(0.0, _new.g(x, u))
            _new.hx = (x, u) -> ForwardDiff.gradient(dx->_new.h(dx,u),x)
            _new.hu = (x, u) -> ForwardDiff.gradient(du->_new.h(x,du),u)
            

            return _new
        end

    end

    function dynamics_rk4(x,u,a,h)
        #RK4 integration with zero-order hold on u
        f1 = RobotZoo.dynamics(a, x, u)
        f2 = RobotZoo.dynamics(a, x + 0.5*h*f1, u)
        f3 = RobotZoo.dynamics(a, x + 0.5*h*f2, u)
        f4 = RobotZoo.dynamics(a, x + h*f3, u)
        return x + (h/6.0)*(f1 + 2.0*f2 + 2.0*f3 + f4)
    end

    mutable struct Marray
        M::Int64
        N::Int64
        L::Vector{Int64}
        _L::Int64
        n::Vector{Int64}
        a::Array{Float64}
        function Marray(_L, M, N, n)
            _new = new()
            _new.M = M
            _new.N = N
            _new._L = _L
            _new.n = n
            _new.L = [(_L-1)*(i-1) for i in 1:M+1]
            _new.L[M+1] = N
            _new.a = zeros(N, n...)
            return _new
        end
        
        function Marray(_L, M, N, n, a)
            _new = new()
            _new.M = M
            _new.N = N
            _new._L = _L
            _new.n = n
            _new.L = [(_L-1)*(i-1) for i in 1:M+1]
            _new.L[M+1] = N
            _new.a = copy(a)
            return _new
        end
    end

    Base.getindex(a::Marray, i::Int64, j::Int64) = begin
        idx = a.L[i]+j
        if idx <= a.N
            return a.a[idx]
        else
            return NaN64
        end
    end
    Base.getindex(a::Marray, _::Colon, _::Colon) = begin
        return a.a[:]
    end
    Base.getindex(a::Marray, i::Int64, j::Int64, k...) = begin
        idx = a.L[i]+j
        if idx <= a.N
            return a.a[idx, k...]
        else
            return NaN64
        end
    end
    Base.getindex(a::Marray, _::Colon, _::Colon, k...) = begin
        return a.a[:, k...]
    end
    Base.setindex!(a::Marray, v, i::Int64, j::Int64) = begin
        idx = a.L[i]+j
        if idx <= a.N
            a.a[idx] = v
        end
    end
    Base.setindex!(a::Marray, v, _::Colon, _::Colon) = begin
        return a.a[:] = v
    end
    Base.setindex!(a::Marray, v, i::Int64, j::Int64, k...) = begin
        idx = a.L[i]+j
        if idx <= a.N
            a.a[a.L[i]+j, k...] = v
        end
    end
    Base.setindex!(a::Marray, v, _::Colon, _::Colon, k...) = begin
        return a.a[:, k...] = v
    end

    Base.copy(s::Marray) = Marray(s._L, s.M, s.N, s.n, s.a)

    # import Base.:+, Base.:-
    Base.:+(a::Marray, b::Marray) = Marray(a._L, a.M, a.N, a.n, a.a + b.a)
    Base.:-(a::Marray, b::Marray) = Marray(a._L, a.M, a.N, a.n, a.a - b.a)

    Base.copyto!(a::Marray, b::Marray) = (a.a = copy(b.a))
    Base.copyto!(a::Marray, b::Marray, i::Int64) = (a.a[a.L[i]+1:a.L[i+1]] = copy(b.a[b.L[i]+1:b.L[i+1]]))
    Base.copyto!(a::Marray, b::Marray, i::Int64, j::Int64) = (a.a[a.L[i]+j] = copy(b.a[b.L[i]+j]))
    Base.copyto!(a::Marray, b::Marray, i::Int64, j::Int64, k...) = (a.a[a.L[i]+j, k...] = copy(b.a[b.L[i]+j, k...]))
    Base.copyto!(a::Marray, b::Array) = (a.a = copy(b))

    # Base.endof(a::Marray) = a.N

    export sumMarray, endMarray, getEndMarray, setEndMarray, getMarray, setMarray

    sumMarray(f::Function, a::Marray, b::Marray) = sum(f(a.a[i, :], b.a[i, :]) for i in 1:min(a.N, b.N))
    sumMarray(f::Function, a::Marray...) = sum(f([ele.a[i, :] for ele in a]...) for i in 1:min([ele.N for ele in a]...))
    endMarray(f::Function, a::Marray) = f(a.a[a.N, :])
    endMarray(f::Function, a::Marray...) = f([ele.a[ele.N, :] for ele in a]...)
    getMarray(a::Marray, b) = a.a[a.N, b...]
    getEndMarray(a::Marray) = a.a[a.N, a.n...]
    getEndMarray(a::Marray, k) = a.a[a.N, [a.n[i] for i in 1:length(a.n)-size(k)[1]]..., k...]
    setMarray!(a::Marray, b, v) = (a.a[a.N, b...] = v)
    setEndMarray!(a::Marray, v) = (a.a[a.N, a.n...] = v)
    setEndMarray!(a::Marray, v, k) = (a.a[a.N, [a.n[i] for i in 1:length(a.n)-size(k)[1]]..., k...] = v)

    getMarray(a::Marray) = a.a

    broadcast(f::Function, a::Marray...) = begin
        idx = argmin([ele.N for ele in a])
        n = size(f(a.a[idx, [1:n for n in 1:a[idx].n]]))
        return Marray(a[idx]._L, a[idx].M, a[idx].N, n, [f(a.a[i, [1:n for n in 1:a.n]]) for i in 1:a[idx].N])
    end

    # norm(a::Marray) = Marray(a._L, a.M, a.N, [], norm.(a.a, 2))

    # max(a::Marray) = begin
    #     idx = argmin([ele.N for ele in a])
    #     Marray(a[idx]._L, a[idx].M, a[idx].N, a[idx].n, maximum(a.a, dims=1))
    # end

    onesMarray(_L, M, N, n) = Marray(_L, M, N, n, ones(N, n...))
    zerosMarray(_L, M, N, n) = Marray(_L, M, N, n, zeros(N, n...))
    similar(a::Marray) = Marray(a._L, a.M, a.N, a.n, similar(a.a))
end
