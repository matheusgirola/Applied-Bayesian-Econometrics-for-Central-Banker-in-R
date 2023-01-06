function  [yhatx,yhat,hhat]  = getforecast( varcoef,horizon,L,LH,N,mumat,Fmat,Qmat,iamat,hlast,y)
LL=max(L,LH);
    hhat=zeros(horizon+LL,N);
    hhat(LL:LL,:)=log(hlast(end,1:N));
    vhat=randn(horizon+LL,N);
    ehat=randn(horizon+LL,N);
    yhat=zeros(horizon+LL,N);
    yhat((LL-L)+1:(LL-L)+L,:)=y(end-L+1:end,:);
    for m=LL+1:horizon+LL
        hhat(m,:)=mumat(1:N)'+hhat(m-1,:)*Fmat(1:N,1:N)+vhat(m,:)*chol(Qmat(1:N,1:N));
        xhat=[];
        for i=1:L
            xhat=[xhat yhat(m-i,:)];
        end
        xhat=[xhat hhat(m,:)];
        for i=1:LH
            xhat=[xhat hhat(m-i,:)];
        end
        xhat=[xhat 1];
        smat=iamat*diag(exp(hhat(m,:)))*iamat';
        
        yhat(m,:)=xhat*varcoef+ehat(m,:)*cholx(smat);
    end
    yhatx=yhat(LL+1:end,:);
end

