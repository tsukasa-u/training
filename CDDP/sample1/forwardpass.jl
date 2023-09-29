using InvertedIndices

# function fp = forwardpass(alg, funcs, fp, bp)
function forwardpass(alg, funcs, fp, bp)
  failed=0
  cost=0
  logcost=0
  err=0
  stepsize=0
    N=fp.horizon;
    dim_x=length(fp.x[:,1]);
    dim_u=length(fp.u[:,1]);
    dim_c=length(fp.c[:,1]);
    xold=fp.x;
    uold=fp.u;
    yold=fp.y;
    sold=fp.s;
    cold=fp.c;
    tau=max(0.99, 1-alg.mu);
    # steplist=2.0.^linspace(0,-10, 11);
    steplist=2.0.^[0,-1,-2,-3,-4,-5,-6,-7,-8,-9,-10]

    xnew=zeros(dim_x,N+1);
    unew=zeros(dim_u,N);
    ynew=zeros(dim_c,N);
    snew=zeros(dim_c,N);
    cnew=zeros(dim_c,N);
    qnew=zeros(1,N);

    step=0
    for step in eachindex(steplist)
      
      xnew=zeros(dim_x,N+1);
      unew=zeros(dim_u,N);
      ynew=zeros(dim_c,N);
      snew=zeros(dim_c,N);
      cnew=zeros(dim_c,N);
      qnew=zeros(1,N);
      
      failed=0;
      stepsize=steplist[step];
      xnew[:,1]=xold[:,1];
      if alg.infeas
        for i in 1:N
          ynew[:,i] = yold[:,i] + stepsize*bp.ky[:,i]+bp.Ky[:,:,i]*(xnew[:,i]-xold[:,i]);
          snew[:,i] = sold[:,i] + stepsize*bp.ks[:,i]+bp.Ks[:,:,i]*(xnew[:,i]-xold[:,i]);
          if or(any(ynew[:,i]<(1-tau)*yold[:,i]), any(snew[:,i]<(1-tau)*sold[:,i]))
            failed=1;
            break;
          end
          unew[:,i] = uold[:,i] + stepsize*bp.ku[:,i]+bp.Ku[:,:,i]*(xnew[:,i]-xold[:,i]);
          xnew[:,i+1] = funcs.f([xnew[:,i]; unew[:,i]]);
        end
      else
        for i in 1:N 
          snew[:,i] = sold[:,i] + stepsize*bp.ks[:,i]+bp.Ks[:,:,i]*(xnew[:,i]-xold[:,i]);
          unew[:,i] = uold[:,i] + stepsize*bp.ku[:,i]+bp.Ku[:,:,i]*(xnew[:,i]-xold[:,i]);
          cnew[:,i] = funcs.c(xnew[:,i], unew[:,i]);
          if any(cnew[:,i]>(1-tau)*cold[:,i]) || any(cnew[:,i].>0)|| any(snew[:,i]<(1-tau)*sold[:,i])
            failed=1;
            println("c ", i)
            break;
          end
          xnew[:,i+1] = funcs.f(xnew[:,i], unew[:,i]);
        end
      end
      # println(cnew)
      
      if failed>0
        continue;
      else

        plot(
          fp.x[1,:],
          fp.x[2,:]
        )
        savefig("cnt.png")

        for i=1:N
          qnew[1,i] = funcs.q(xnew[:,i], unew[:,i]);
        end
        cost=sum(qnew[1,i] for i in 1:N)+funcs.p(xnew[:,N+1]);
        
        if alg.infeas
          logcost=cost-alg.mu*sum(log(reshape(ynew,1,[])));   
          for i=1:N
            cnew[:,i] = funcs.c([xnew[:,i]; unew[:,i]]);
          end
          err=max(alg.tol, norm(reshape(cnew+ynew,1,[]),1));
        else
          for i in 1:N
            for j in 1:dim_c
              if cnew[j,i]>0
                println("cnew ", i, " ", j, " ", cnew[j,i], " ", xnew[:,i], " ", unew[:,i])
              end
            end
          end
          logcost=sum(qnew[1,i] for i in 1:N)+funcs.p(xnew[:,N+1])-alg.mu*sum(log.(reshape(-cnew,1,:)[1,:]));
          err=0;
        end
        
        candidate=[logcost;err];
        if any(all(candidate>=fp.filter))
          failed=2;
          continue;
        else
          idx=all(candidate<=fp.filter);
          for (i, ele) in enumerate(idx)
            if ele
              # deleteat!(fp.filter, [:,idx])
              # fp.filter[:,idx]=[]
              # println(fp.filter[Not(i)], " ", fp.filter)
              fp.filter=fp.filter[Not(i)]
            end
          end
          fp.filter=[fp.filter; candidate];
          break;
        end      
      end
    end
    
    if failed==1
        fp.failed=failed;
        fp.stepsize=0;
    else
        fp.cost=cost;
        fp.logcost=logcost;
        fp.x=xnew;
        fp.u=unew;
        fp.y=ynew;
        fp.s=snew;
        fp.c=cnew;
        fp.q=qnew;
        fp.err=err;
        fp.stepsize=stepsize;
        fp.step=step;
        fp.failed=0;
    end
    return fp
  end