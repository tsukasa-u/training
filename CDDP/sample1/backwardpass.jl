function tensdot(a, b)
  dim_a = size(a)
  dim_b = size(b)
  # println(dim_a, " ", dim_b)
  if length(dim_a) == 3
    return sum(a[i,:,:]*b[i] for i in 1:dim_a[3])
  elseif length(dim_b) == 3
    return sum(a[i]*b[i,:,:] for i in 1:dim_b[3])
  else
    return a * b
  end
end

# function bp = backwardpass(alg, funcs, fp, bp)
function backwardpass(alg, funcs, fp, bp)
    # tensdot = @(a,b) permute(sum(bsxfun(@times,a,b),1), [2 3 1]);
    N=fp.horizon;
    dim_x=length(fp.x[:,1]);
    dim_u=length(fp.u[:,1]);
    dim_c=length(fp.c[:,1]);
  
    dV=[0;0];
    c_err=0;
    mu_err=0;
    Qu_err=0;
  
    if fp.failed>0 || bp.failed>0
       bp.reg=bp.reg+1;
    elseif fp.step==1
       bp.reg=bp.reg-1;
    elseif fp.step<=4
       bp.reg=bp.reg;
    else
       bp.reg=bp.reg+1;
    end
    
    if bp.reg<0
      bp.reg=0;
    elseif bp.reg>24
      bp.reg=24;
    end
  
    V=funcs.p(fp.x[:,N+1]);
    Vx=funcs.px(fp.x[:,N+1]);
    Vxx=funcs.pxx(fp.x[:,N+1]);
    for i=N:-1:1
      x=fp.x[:,i];
      u=fp.u[:,i];
      s=fp.s[:,i];
      y=fp.y[:,i];
      q=fp.q[:,i];
      c=fp.c[:,i]; 
  
      fx=funcs.fx(x,u); 
      fu=funcs.fu(x,u);
  
      Qsx=funcs.cx(x,u); 
      Qsu=funcs.cu(x,u);
  
      # println(size(funcs.fxx_(x,u)))
      # println(size(funcs.fxu_(x,u)))
      # println(size(funcs.fuu_(x,u)))
      fxx=funcs.fxx(x,u); 
      fxu=funcs.fxu(x,u); 
      fuu=funcs.fuu(x,u);
  
      qx=funcs.qx(x,u); 
      qu=funcs.qu(x,u);
  
      Qx=qx+transpose(Qsx)*s+transpose(fx)*Vx;
      Qu=qu+transpose(Qsu)*s+transpose(fu)*Vx; 
  
      Qxx=funcs.qxx(x,u)+transpose(fx)*Vxx*fx+tensdot(Vx,fxx);
      # println(size(funcs.qxu(x,u)), " ", size(transpose(fx)*Vxx*fu), " ", size(tensdot(Vx,fxu)))
      Qxu=funcs.qxu(x,u)+transpose(fx)*Vxx*fu+tensdot(Vx,fxu);
      quu=funcs.quu(x,u);
      Quu=quu+transpose(fu)*Vxx*fu+tensdot(Vx,fuu);
  
    #   if isfield(funcs, 'cxx')
        Qxx=Qxx+tensdot(funcs.cxx(x,u), s);
    #   end
    #   if isfield(funcs, 'cxu')
        Qxu=Qxu+tensdot(funcs.cxu(x,u), s);
    #   end
    #   if isfield(funcs, 'cuu')
        Quu=Quu+tensdot(funcs.cuu(x,u), s);
    #   end
      
      S=diagm(s);
    
      Quu_reg=Quu+quu*(1.6^bp.reg-1);
  
      if alg.infeas
          r=s.*y-alg.mu;
          rhat=s.*(c+y)-r;
          yinv=1.0./y;
          SYinv=diagm(s.*yinv);
          
        #   [R, failed]=chol(Quu_reg+Qsu'*SYinv*Qsu);
        #   if failed
        #       bp.failed=1;
        #       break;
        #   end
        #   kK=-R\(R'\[Qu+Qsu'*(yinv.*rhat), Qxu'+Qsu'*SYinv*Qsx]);
          R = Quu_reg+transpose(Qsu)*SYinv*Qsu
          kK=-R\[Qu+transpose(Qsu)*(yinv.*rhat), transpose(Qxu)+transpose(Qsu)*SYinv*Qsx]

          ku= kK[:,1];
          Ku= kK[:,2:end];
          ks= yinv.*(rhat+S*Qsu*ku);
          Ks= SYinv*(Qsx+Qsu*Ku);
          ky= -(c+y)-Qsu*ku;
          Ky= -Qsx-Qsu*Ku;
        
          Quu=Quu+transpose(Qsu)*SYinv*Qsu;
          Qxu=Qxu+transpose(Qsx)*SYinv*Qsu;
          Qxx=Qxx+transpose(Qsx)*SYinv*Qsx;
  
          Qu=Qu+transpose(Qsu)*(yinv.*rhat);
          Qx=Qx+transpose(Qsx)*(yinv.*rhat);
      else
          r=S*c.+alg.mu;
          cinv=1.0./c;
          SCinv=diagm(s.*cinv);
          
        #   [R, failed]=chol(Quu_reg-Qsu'*SCinv*Qsu);
        #   if failed
        #       bp.failed=1;
        #       break;
        #   end
        #   kK=-R\(R'\[Qu-Qsu'*(cinv.*r), Qxu'-Qsu'*SCinv*Qsx]);

          R = Quu_reg-Qsu'*SCinv*Qsu
          kK=-R\[Qu-transpose(Qsu)*(cinv.*r) transpose(Qxu)-transpose(Qsu)*SCinv*Qsx]

          ku= kK[:,1];
          Ku= kK[:,2:end];
          # println(size(cinv), " ", size(r), " ", size(S), " ", size(Qsu), " ", size(ku))
          ks= -cinv.*(r+S*Qsu*ku);
          Ks= -SCinv*(Qsx+Qsu*Ku);
          ky= zeros(dim_c,1);
          Ky= zeros(dim_c,dim_x);
          
          Quu=Quu-transpose(Qsu)*SCinv*Qsu;
          Qxu=Qxu-transpose(Qsx)*SCinv*Qsu;
          Qxx=Qxx-transpose(Qsx)*SCinv*Qsx;
  
          Qu=Qu-transpose(Qsu)*(cinv.*r);
          Qx=Qx-transpose(Qsx)*(cinv.*r);
      end
      dV=dV+ [transpose(ku)*Qu; 0.5*transpose(ku)*Quu*ku];
      Vx=Qx+transpose(Ku)*Qu+transpose(Ku)*Quu*ku + Qxu*ku;
      Vxx=Qxx+transpose(Ku)*transpose(Qxu)+Qxu*Ku+ transpose(Ku)*Quu*Ku;
     
      bp.ku[:,i]=ku;
      bp.ky[:,i]=ky;
      bp.ks[:,i]=ks;
  
      bp.Ku[:,:,i]=Ku;
      bp.Ky[:,:,i]=Ky;
      bp.Ks[:,:,i]=Ks;
  
    #   %% Optimality error
      Qu_err=max(Qu_err,norm(Qu, Inf));
      mu_err=max(mu_err,norm(r, Inf));
      if alg.infeas
        c_err=max(c_err, norm(fp.c[:,i]+fp.y[:,i],Inf));
      end
    end
    bp.failed=0;
    bp.opterr=max(Qu_err, c_err, mu_err);
    bp.dV=dV;
    return bp
  end