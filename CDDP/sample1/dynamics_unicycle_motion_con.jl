using LinearAlgebra
using ForwardDiff

mutable struct FP
    x::Matrix{Float64}
    u::Matrix{Float64}
    y::Matrix{Float64}
    s::Matrix{Float64}
    c::Matrix{Float64}
    q::Matrix{Float64}
    mu::Matrix{Float64}
    horizon::Int64
    filter::Vector{Float64}
    step::Int64
    stepsize::Float64
    failed::Int64
    cost::Float64
    logcost::Float64
    err::Float64

    FP() = new()
end

mutable struct BP
    ku::Array{Float64, 2}
    Ku::Array{Float64, 3}
    ky::Array{Float64, 2}
    Ky::Array{Float64, 3}
    ks::Array{Float64, 2}
    Ks::Array{Float64, 3}

    reg::Float64
    failed::Int64
    recovery::Int64
    opterr::Float64
    dV::Vector{Float64}

    BP() = new()
end

mutable struct FUNC
    f::Function
    fx::Function
    fu::Function
    fxx::Function
    fxu::Function
    fuu::Function
    fxx_::Function
    fxu_::Function
    fuu_::Function

    q::Function
    qx::Function
    qu::Function
    qxx::Function
    qxu::Function
    quu::Function

    p::Function
    px::Function
    pxx::Function
   
    c::Function
    cx::Function
    cu::Function
    cxx::Function
    cxu::Function
    cuu::Function

    dim_c::Int64
    dim_x::Int64
    dim_u::Int64

    FUNC() = new() 
end

# function [funcs, fp, bp] = dynamics_concar(seed)
function dynamics_unicycle_motion_con()
    N=600;
    dim_x=3;
    dim_u=1;   
    dim_c=7; 
    
    # x = sym('x_%d', [dim_x 1], 'real');
    # u = sym('u', [dim_u 1],'real');
    
    # d  = 2.0;
    # h  = 0.03;

    # f = x+[ (d + h*x(4)*cos(u(1)) - sqrt(d^2 - (h*x(4)*sin(u(1)))^2))*[cos(x(3));sin(x(3))];
    #     asin(sin(u(1))*h*x(4)/d);
    #     h*u(2)];
    f(x, u) = begin
        rx, ry, ϕ = x
        u_, = u
        h = 0.1
        v = 0.15
        return [
            h*v*cos(ϕ),
            h*v*sin(ϕ),
            h*u_
        ] + x
    end
    
    

    # fx = cat(2, diff(f,x(1)), diff(f,x(2)), diff(f,x(3)), diff(f,x(4)));
    # fu = cat(2, diff(f,u(1)), diff(f,u(2)));
    # fxx= cat(3, diff(fx,x(1)), diff(fx,x(2)), diff(fx,x(3)), diff(fx,x(4)));
    # fxu= cat(3, diff(fx,u(1)), diff(fx,u(2)));
    # fuu= cat(3, diff(fu,u(1)), diff(fu,u(2)));
    fx(x, u) = ForwardDiff.jacobian(x -> f(x, u), x)
    fu(x, u) = ForwardDiff.jacobian(u -> f(x, u), u)

    function reshape_tensor_3(x, n::Tuple{S, S, S}) where S<:Integer
        return [x[(i-1)*n[1]+v, j] for v in 1:n[1], i in 1:n[2], j in 1:n[3]]
    end

    fxx_(x, u) = ForwardDiff.jacobian(x -> ForwardDiff.jacobian(x -> f(x, u), x), x)
    fxu_(x, u) = ForwardDiff.jacobian(x -> ForwardDiff.jacobian(u -> f(x, u), u), x)
    fuu_(x, u) = ForwardDiff.jacobian(u -> ForwardDiff.jacobian(u -> f(x, u), u), u)
    fxx(x, u) = reshape_tensor_3(ForwardDiff.jacobian(x -> ForwardDiff.jacobian(x -> f(x, u), x), x), (dim_x,dim_x,dim_x))
    fxu(x, u) = reshape_tensor_3(ForwardDiff.jacobian(u -> ForwardDiff.jacobian(x -> f(x, u), x), u), (dim_x,dim_x,dim_u))
    fuu(x, u) = reshape_tensor_3(ForwardDiff.jacobian(u -> ForwardDiff.jacobian(u -> f(x, u), u), u), (dim_x,dim_u,dim_u))

    

    # R=1e-2*diag([1 .01]);
    # Q=1e-3*[1  1]; 
    
    # q=u'*R*u + 1e-3*(sqrt(x(1)^2+Q(1))-Q(1))+1e-3*(sqrt(x(2)^2+Q(2))-Q(2));   #stage cost
    # qx=gradient(q,x);
    # qu=gradient(q,u);
    # qxx=cat(2, diff(qx,x(1)), diff(qx,x(2)), diff(qx,x(3)), diff(qx,x(4)));
    # qxu=cat(2, diff(qx,u(1)), diff(qx,u(2)));
    # quu=cat(2, diff(qu,u(1)), diff(qu,u(2)));
    q(x, u) = begin
        rx, ry, ϕ = x
        u_, = u
        return 0.1*(rx^2+ry^2+ϕ^2+0.1*u_^2)
    end
    qx(x, u) = ForwardDiff.gradient(x -> q(x, u), x)
    qu(x, u) = ForwardDiff.gradient(u -> q(x, u), u)
    qxx(x, u) = ForwardDiff.hessian(x -> q(x, u), x)
    qxu(x, u) = ForwardDiff.jacobian(u -> ForwardDiff.gradient(x -> q(x, u), x), u)
    quu(x, u) = ForwardDiff.hessian(u -> q(x, u), u)
    
    # P=[.01 .01 .01  1];
    # p=0.1*(sqrt(x(1)^2+P(1))-P(1))+0.1*(sqrt(x(2)^2+P(2))-P(2))+...
    #     1*(sqrt(x(3)^2+P(3))-P(3))+0.3*(sqrt(x(4)^2+P(4))-P(4));
    # px=gradient(p,x);
    # pxx=cat(2, diff(px,x(1)), diff(px,x(2)), diff(px,x(3)), diff(px,x(4)));
    p(x) = begin
        rx, ry, ϕ = x
        return 0.1*(rx^2+ry^2+ϕ^2)
    end
    px(x) = ForwardDiff.gradient(x -> p(x), x)
    pxx(x) = ForwardDiff.hessian(x -> p(x), x)
    
    # c=[x(1)-2; -x(1)-2; x(2)-2; -x(2)-2; -u(1)-.5; -u(1)-.5; u(2)-2; -u(2)-2;]; #constraints
    # # %c=[u-0.1; -u-0.1; 0.5+x(2); 0.5-x(2)]; %constraints
    # cx= cat(2, diff(c,x(1)), diff(c,x(2)), diff(c,x(3)), diff(c,x(4)));
    # cu= cat(2, diff(c,u(1)), diff(c,u(2)));
    c(x, u) = begin
        rx, ry, ϕ = x
        u_, = u
        return [
            # u_ - 1.5,
            # -u_ - 1.5,
            # ry - 1.0,
            # -ry - 1.0,
            # 1.0^2 - ((rx + 5.0)^2 + (ry + 1.0)^2),
            # 0.5^2 - ((rx + 8.0)^2 + (ry - 0.2)^2),
            # 1.5^2 - ((rx + 2.5)^2 + (ry - 1.0)^2),
            0,0,0,0,0,0,0
        ]
    end
    cx(x, u) = ForwardDiff.jacobian(x -> c(x, u), x)
    cu(x, u) = ForwardDiff.jacobian(u -> c(x, u), u)
    cxx(x, u) = reshape_tensor_3(ForwardDiff.jacobian(x -> ForwardDiff.jacobian(x -> c(x, u), x), x), (dim_c,dim_x,dim_x))
    cxu(x, u) = reshape_tensor_3(ForwardDiff.jacobian(u -> ForwardDiff.jacobian(x -> c(x, u), x), u), (dim_c,dim_x,dim_u))
    cuu(x, u) = reshape_tensor_3(ForwardDiff.jacobian(u -> ForwardDiff.jacobian(u -> c(x, u), u), u), (dim_c,dim_u,dim_u))
    
    
    # vars=[x;u];
    
    # funcs.f=matlabFunction(f, 'vars', {vars});
    # funcs.fx=matlabFunction(fx, 'vars', {vars});
    # funcs.fu=matlabFunction(fu, 'vars', {vars});
    # funcs.fxx=matlabFunction(fxx, 'vars', {vars});
    # funcs.fxu=matlabFunction(fxu, 'vars', {vars});
    # funcs.fuu=matlabFunction(fuu, 'vars', {vars});

    # funcs.q=matlabFunction(q, 'vars', {vars});
    # funcs.qx=matlabFunction(qx, 'vars', {vars});
    # funcs.qu=matlabFunction(qu, 'vars', {vars});
    # funcs.qxx=matlabFunction(qxx, 'vars', {vars});
    # funcs.qxu=matlabFunction(qxu, 'vars', {vars});
    # funcs.quu=matlabFunction(quu, 'vars', {vars});

    # funcs.p=matlabFunction(p, 'vars', {vars});
    # funcs.px=matlabFunction(px, 'vars', {vars});
    # funcs.pxx=matlabFunction(pxx, 'vars', {vars});
   
    # if dim_c==0
    #   funcs.c=[];
    #   funcs.cx=[];
    #   funcs.cu=[];
    # else
    #   funcs.c=matlabFunction(c, 'vars', {vars});
    #   funcs.cx=matlabFunction(cx, 'vars', {vars});
    #   funcs.cu=matlabFunction(cu, 'vars', {vars});
    # end

    funcs = FUNC()
    funcs.f = f
    funcs.fx = fx
    funcs.fu = fu
    funcs.fxx = fxx
    funcs.fxu = fxu
    funcs.fuu = fuu
    funcs.fxx_ = fxx_
    funcs.fxu_ = fxu_
    funcs.fuu_ = fuu_

    funcs.q = q
    funcs.qx = qx
    funcs.qu = qu
    funcs.qxx = qxx
    funcs.qxu = qxu
    funcs.quu = quu

    funcs.p = p
    funcs.px = px
    funcs.pxx = pxx

    funcs.c = c
    funcs.cx = cx
    funcs.cu = cu
    funcs.cxx = cxx
    funcs.cxu = cxu
    funcs.cuu = cuu
    
    funcs.dim_c=dim_c;
    funcs.dim_x=dim_x;
    funcs.dim_u=dim_u;

    fp = FP()
    
    fp.x=zeros(dim_x, N+1);
    fp.x[:,1]=[-10;0;3.15];

    fp.u=0.02*zeros(dim_u, N);
    fp.y=0.01*ones(dim_c, N);
    fp.s=0.1*ones(dim_c, N);
    fp.c=zeros(dim_c, N);
    fp.mu=fp.y.*fp.s;
    fp.horizon=N;
    fp.filter=[Inf;0];
    fp.q=zeros(1,N);

    bp = BP()
    bp.ku=zeros(dim_u,N);
    bp.Ku=zeros(dim_u,dim_x,N);
    bp.ky=zeros(dim_c,N);
    bp.Ky=zeros(dim_c,dim_x,N);
    bp.ks=zeros(dim_c,N);
    bp.Ks=zeros(dim_c,dim_x,N);

    return funcs, fp, bp
end