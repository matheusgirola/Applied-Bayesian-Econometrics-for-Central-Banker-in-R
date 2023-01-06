
function[dens,hnew,rnew,pmatnew]= getdens(b0,parth,...
    b0r,partr,...
    Sbig,b0p,partp,N,L,EX,ns,y,MU,F,H,MUr,Fr,Qr,...
    nfact,varcoef,iAA,bload,rho)
%draw states
        rtemp=MUr+b0r.*Fr+partr.*sqrt(Qr);
        rtemplag=b0r; %lagged state
        rnew=rtemp;
        
        htemp=MU'+b0*F'+parth;
        hnew=htemp;
        
        %volatilities
        
        sigma=iAA*diag(Sbig.*exp(htemp(:,1)))*iAA';
        csigma=zeros(ns,ns);
        csigma(1:nfact,1:nfact)=cholx(sigma);
        xtemp=zeros(1,N*L+EX);
        xtemp(:,(N*L)+2:end-1)=b0(:,2:end);
        xtemp(:,(N*L)+1:(N*L)+1)=htemp(:,1);
        xtemp(:,end)=1;
        exo=xtemp(:,N*L+1:end)*varcoef(N*L+1:end,:);
        mup=zeros(1,ns);
        mup(:,1:nfact)=exo;
        fp=comp(vec(varcoef),N,L,EX);
        
      
        pmattemp=mup+b0p*fp'+partp*csigma;
        
        pmatnew=pmattemp;
        
        %likelihood
        exor=bload.*(rtemp-rtemplag.*rho);
        yhat=(pmattemp*H');
       
        res=y-yhat-exor;
        RR=diag(exp(rtemp));
        dR=logdet(RR);
        iRR=diag(1./exp(rtemp));
        
            resx=res*iRR*res';
%             densx=(1/sqrt(dR))*exp(-0.5*resx);
denslog=-0.5*dR+(-0.5*resx)    ;

tempdens=(denslog);
    if isnan(tempdens) || isinf(tempdens) || ~isreal(tempdens)
        dens=exp(-100);
    else
        dens=tempdens;
    end